"""Convenience wrapper for looping and clipping identically gridded in inputs

Everything here can be customized using the classes in the subset package

"""
import sys
import argparse
import os
from pathlib import Path
import logging
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from pfsubset.subset.clipper import MaskClipper, BoxClipper
from pfsubset.subset.utils.arguments import is_valid_file, is_valid_path
from pfsubset.subset import TIF_NO_DATA_VALUE_OUT as NO_DATA
from pfsubset.subset.mask import SubsetMask
import pfsubset.subset.utils.io as file_io_tools


def parse_args(args) -> argparse.Namespace:
    """Parse the command line arguments

    Parameters
    ----------
    args : list
        list of arguments from sys.argv

    Returns
    -------
    Namespace
        populated Namespace object from the parsed command line options
    """
    parser = argparse.ArgumentParser('Clip a list of identically gridded files and extract the data within the mask')

    exclusive_group = parser.add_mutually_exclusive_group(required=True)

    exclusive_group.add_argument("--maskfile", "-m", dest="mask_file", required=False,
                        type=lambda x: is_valid_file(parser, x),
                        help="gridded full_dim_mask file to full extent of files to be clipped")

    exclusive_group.add_argument("--bboxfile", "-b", dest="bbox_file", required=False,
                        type=lambda x: is_valid_file(parser, x),
                        help="a tab separated text file indicating x1,y1,nx,ny of the files to be clipped")

    exclusive_group.add_argument("--inline-bbox", "-l", dest="bbox_def", nargs=4, metavar=('X1', 'Y1', 'NX', 'NY'),
                                 required=False, type=int, help="bbox defined by x1 y1 nx ny")

    exclusive_group2 = parser.add_mutually_exclusive_group(required=True)

    exclusive_group2.add_argument("--datafiles", "-d", dest="data_files", required=False,
                                nargs='+',
                                type=lambda x: is_valid_file(parser, x),
                                help="the list of gridded data files (.pfb or .tif) to clip from")

    exclusive_group2.add_argument("--glob", "-g", dest="glob_pattern", required=False, type=str,
                                  help="a filter string for filtering files to clip")

    parser.add_argument("--input_path", "-i", dest="input_path", required=False, default='.',
                        type=lambda x: is_valid_path(parser, x),
                        help="the folder to look for input files")

    parser.add_argument("--ref_file", "-r", dest="ref_file", required=False,
                        type=lambda x: is_valid_file(parser, x),
                        help="the reference tif, if writing TIF outputs")

    parser.add_argument("--out_dir", "-o", dest="out_dir", required=False,
                        default='.',
                        help="the directory to write outputs to",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("--pfb_outs", "-p", dest="write_pfbs", required=False,
                        action='store_true', default=True, help="write pfb output files")

    parser.add_argument("--tif_outs", "-t", dest="write_tifs", required=False,
                        action='store_true', help="write tif output files")

    return parser.parse_args(args)


def mask_clip(mask_file, data_files, out_dir='.', pfb_outs=1, tif_outs=0) -> None:
    """clip a list of files using a full_dim_mask and a domain reference tif

    Parameters
    ----------
    mask_file : str
        full_dim_mask file generated from shapefile to mask utility no_data,0's=bbox,1's=mask
    data_files : list
        list of data files (tif, pfb) to clip from
    out_dir : str, optional
        output directory (optional) (Default value = '.')
    pfb_outs : int, optional
        write pfb files as outputs (optional) (Default value = 1)
    tif_outs : int, optional
        write tif files as outputs (optional) (Default value = 0)

    Returns
    -------
    None
    """
    # read the full_dim_mask file
    mask = SubsetMask(mask_file)

    # create clipper with full_dim_mask
    clipper = MaskClipper(mask_file, no_data_threshold=-1)
    # clip all inputs and write outputs
    clip_inputs(clipper, input_list=data_files, out_dir=out_dir, pfb_outs=pfb_outs,
                tif_outs=tif_outs)


def box_clip(bbox, data_files, out_dir='.', pfb_outs=1, tif_outs=0) -> None:
    """clip a list of files using a bounding box

    Parameters
    ----------
    bbox : tuple
        tuple of x, y, nx, ny to specifies the bounding region
    data_files : list
        list of data files (tif, pfb) to clip from
    out_dir : str, optional
        output directory (optional) (Default value = '.')
    pfb_outs : int, optional
        write pfb files as outputs (optional) (Default value = 1)
    tif_outs : int, optional
        write tif files as outputs (optional) (Default value = 0)

    Returns
    -------
    None
    
    """

    # create clipper with bbox
    data_file = ref_array=file_io_tools.read_file(data_files[0]),
    clipper = BoxClipper( x=bbox[0], y=bbox[1], nx=bbox[2], ny=bbox[3])
    # clip all inputs and write outputs
    clip_inputs(clipper, input_list=data_files, out_dir=out_dir, pfb_outs=pfb_outs,
                tif_outs=tif_outs)


def locate_tifs(file_list) -> list:
    """identify the .tif files in a list of files

    Parameters
    ----------
    file_list : list
        list of files to parse

    Returns
    -------
    list
        list of files where the extension is .tif

    """
    return list([s for s in file_list if '.tif' in s.lower()])


def _clip(clipper, data_file, out_dir, pfb_outs, tif_outs, output_suffix, ref_proj, no_data):
    # A top-level function that is capable of being serialized and executed across cores.
    # The arguments and semantics of the inputs are identical to the public `clip_inputs` function.
    filename = Path(data_file).stem
    #return_arr, new_geom, _, _ = clipper.subset(file_io_tools.read_file(data_file))
    return_arr, new_geom, _, _ = clipper.subset(data_file)
    if pfb_outs:
        file_io_tools.write_pfb(return_arr, os.path.join(out_dir, f'{filename}{output_suffix}.pfb'))
    if tif_outs and new_geom is not None and ref_proj is not None:
        file_io_tools.write_array_to_geotiff(os.path.join(out_dir, f'{filename}{output_suffix}.tif'),
                                             return_arr, new_geom, ref_proj, no_data=no_data)


def clip_inputs(clipper, input_list, out_dir='.', pfb_outs=1, tif_outs=0, no_data=NO_DATA, output_suffix='_clip',
                n_workers=1) -> None:
    """clip a list of files using a clipper object

    Parameters
    ----------
    clipper : Clipper
        clipper object prepared with full_dim_mask and reference dataset
    input_list : list
        list of data files (tif, pfb) to clip from
    out_dir : str, optional
        output directory (optional) (Default value = '.')
    pfb_outs : int, optional
        write pfb files as outputs (optional) (Default value = 1)
    tif_outs : int, optional
        write tif files as outputs (optional) (Default value = 0)
    no_data : int, optional
        no_data value for tifs (optional) (Default value = NO_DATA)
    output_suffix : str, optional
        filename suffix to add to all output pfb/tif files
    n_workers: int, optional
        No. of parallel workers to launch for clipping files (Default value = 1)
    Returns
    -------
    None
    """
    ref_proj = None
    if tif_outs:
        # identify projection
        ref_proj = clipper.subset_mask.mask_tif.GetProjection()

    with ProcessPoolExecutor(max_workers=min(n_workers, len(input_list))) as executor:
        to_do = []
        for data_file in input_list:
            future = executor.submit(_clip, clipper, data_file, out_dir, pfb_outs, tif_outs, output_suffix, ref_proj, no_data)
            to_do.append(future)

        for future in as_completed(to_do):
           # Even though we don't need the result, trying to obtain future.result() will re-raise Exceptions, if any.
           future.result()


def get_file_list(input_dir, glob_pattern=None, files=None) -> list:
    """get a list of proper paths for files either in the list or matching the glob pattern

    Parameters
    ----------
    input_dir : Path
        directory where input files are located
    glob_pattern : str
        filename pattern to match
    files : list
        list of input filenames to use
    Returns
    -------
    file_list : list of Paths
        the assembled list of filenames with paths
    """
    file_list = []
    if glob_pattern is not None:
        file_list = input_dir.glob(glob_pattern)
    elif files is not None:
        file_list = [input_dir / filename for filename in files]
    return list(file_list)


def main():
    # setup logging
    logging.basicConfig(filename='bulk_clipper.log', filemode='w', level=logging.INFO)
    start_date = datetime.utcnow()
    logging.info(f'start process at {start_date} from command {" ".join(sys.argv[:])}')
    args = parse_args(sys.argv[1:])
    # get the input path and data file list
    input_path = Path(args.input_path)
    data_files = get_file_list(input_dir=input_path, files=args.data_files, glob_pattern=args.glob_pattern)
    # If tif out specified, look for a reference tif
    if args.write_tifs and not args.ref_file:
        if 'tif' not in args.mask_file.lower():
            input_tifs = locate_tifs(data_files)
            if len(input_tifs) < 1:
                raise Exception('Must include at least one geotif input or a ref_file when tif_outs is selected')
    if args.mask_file:
        mask_clip(args.mask_file, data_files, args.out_dir, args.write_pfbs, args.write_tifs)
    elif args.bbox_file:
        box_clip(file_io_tools.read_bbox(args.bbox_file), data_files, args.out_dir, args.write_pfbs,
                 args.write_tifs)
    elif args.bbox_def:
        box_clip(args.bbox_def, data_files, args.out_dir, args.write_pfbs, args.write_tifs)
    end_date = datetime.utcnow()
    logging.info(f'completed process at {end_date} for a runtime of {end_date-start_date}')


if __name__ == '__main__':
    main()
