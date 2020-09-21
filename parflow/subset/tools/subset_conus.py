"""Convenience wrapper for clipping inputs for the CONUS1 or CONUS2 models

Everything here can be customized using the classes in the subset package
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from parflow.subset.utils.arguments import is_valid_path, is_positive_integer, is_valid_file
from parflow.subset.clipper import MaskClipper
from parflow.subset.domain import Conus
from parflow.subset.rasterizer import ShapefileRasterizer
from datetime import datetime
import parflow.subset.tools.bulk_clipper as bulk_clipper
import parflow.subset.builders.solidfile as solidfile_generator
from parflow.subset.builders.tcl import build_tcl
from parflow.subset.clipper import ClmClipper
from parflow.subset.data import parkinglot_template, conus_manifest


def parse_args(args):
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
    parser = argparse.ArgumentParser('Subset a ParFlow CONUS domain')

    """
    Required Arguments
    """
    parser.add_argument("--input_path", "-i", dest="input_path", required=True,
                        type=lambda x: is_valid_path(parser, x),
                        help="the input path to the shapefile file set")

    parser.add_argument("--shapefile", "-s", dest="shapefile", required=True,
                        help="the name of the shapefile file set")

    parser.add_argument("--conus_files", "-f", dest="conus_files", required=True,
                        help="local path to the CONUS inputs to subset",
                        type=lambda x: is_valid_path(parser, x))
    """
    Optional Arguments
    """
    parser.add_argument("--manifest", "-m", dest="manifest_file", required=False,
                        default=conus_manifest,
                        type=lambda x: is_valid_file(parser, x),
                        help="the manifest of CONUS filenames to build the domain from")

    parser.add_argument("--version", "-v", dest="conus_version", required=False,
                        default=1, choices=[1, 2],
                        type=lambda x: is_positive_integer(parser, x),
                        help="the version of CONUS to subset from")

    parser.add_argument("--out_dir", "-o", dest="out_dir", required=False,
                        default='.',
                        help="the directory to write outputs to",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("--out_name", "-n", dest="out_name", required=False,
                        default=None,
                        help="the name to give the outputs")

    parser.add_argument("--clip_clm", "-c", dest="clip_clm", required=False,
                        action='store_true',
                        help="also clip inputs for CLM")

    parser.add_argument("--write_tcl", "-w", dest="write_tcl", required=False,
                        action='store_true',
                        help="generate the .tcl script for this subset")

    parser.add_argument("--padding", "-p", dest="padding", nargs=4, metavar=('Top', 'Right', 'Bottom', 'Left'),
                        required=False, default=(0, 0, 0, 0), type=int,
                        help="integer padding value for bounding box on x-sides"
                        )

    parser.add_argument("--attribute_ids", "-a", dest="attribute_ids", required=False,
                        default=[1], nargs='+',
                        help="list of attribute ID's to clip",
                        type=lambda x: is_positive_integer(parser, x))

    parser.add_argument("--attribute_name", "-e", dest="attribute_name", required=False,
                        default="OBJECTID",
                        help="name of the attribute field to query for attribute ids",
                        type=str)

    parser.add_argument("--tif_outs", "-t", dest="write_tifs", required=False,
                        action='store_true', help="write tif output files")

    return parser.parse_args(args)


def subset_conus(input_path, shapefile, conus_version=1, conus_files='.', out_dir='.', out_name=None, clip_clm=False,
                 write_tcl=False, padding=(0, 0, 0, 0), attribute_name='OBJECTID', attribute_ids=None, write_tifs=False,
                 manifest_file=conus_manifest):
    """subset a conus domain inputs for running a regional model

    Parameters
    ----------
    input_path : str
        path to input shapefile parts to use as mask
    shapefile : str
        name of shapefile to use as mask
    conus_version : int, optional
        version of the CONUS domain to use (1 or 2) (Default value = 1)
    conus_files : str, optional
        path to the CONUS source input files listed in conus_manifest.yaml (Default value = '.')
    out_dir : str, optional
        directory to write the outputs (default .)
    out_name : str, optional
        name to give the outputs (default shapefile name)
    clip_clm : int, optional
        whether or not to clip the CLM input files too (default no)
    write_tcl : int, optional
        whether or not to write a TCL file for the subset (default no)
    padding : tuple, optional
        grid cells of no_data to add around domain mask. CSS Style (top, right, bottom, left) default 0
    attribute_name : str, optional
        attribute name defined in shapefile to select as mask default 'OBJECTID'
    attribute_ids : list, optional
        list of attribute ID's defined in shapefile to use as mask input. default [1]
    write_tifs : int, optional
        whether or not to write outputs as TIF's in addition to PFB's. (default no)

    Returns
    -------
    None
    """
    if out_name is None:
        out_name = shapefile
    conus = Conus(version=conus_version, local_path=conus_files, manifest_file=manifest_file)

    if attribute_ids is None:
        attribute_ids = [1]
    # Step 1, rasterize shapefile

    rasterizer = ShapefileRasterizer(input_path, shapefile, reference_dataset=conus.get_domain_tif(),
                                     no_data=-999, output_path=out_dir, )
    rasterizer.rasterize_shapefile_to_disk(out_name=f'{out_name}_raster_from_shapefile.tif',
                                           padding=padding,
                                           attribute_name=attribute_name,
                                           attribute_ids=attribute_ids)

    subset_mask = rasterizer.subset_mask

    # Step 2, Generate solid file
    clip = MaskClipper(subset_mask, no_data_threshold=-1)
    batches = solidfile_generator.make_solid_file(clipped_mask=clip.clipped_mask,
                                                  out_name=os.path.join(out_dir, out_name))
    if len(batches) == 0:
        raise Exception("Did not make solid file correctly")

    # Step 3. Clip all the domain data inputs
    bulk_clipper.clip_inputs(clip,
                             [os.path.join(conus.local_path, value) for key, value in conus.required_files.items()
                              if key not in ['DOMAIN_MASK', 'CHANNELS']],
                             out_dir=out_dir, tif_outs=write_tifs)

    # Step 4. Clip CLM inputs
    if clip_clm == 1:
        clm_clipper = ClmClipper(subset_mask)
        latlon_formatted, latlon_data = clm_clipper.clip_latlon(os.path.join(conus.local_path,
                                                                             conus.optional_files.get('LAT_LON')))

        clm_clipper.write_lat_lon(latlon_formatted, os.path.join(out_dir, f'{out_name}_latlon.sa'),
                                  x=latlon_data.shape[2], y=latlon_data.shape[1], z=latlon_data.shape[0])

        land_cover_data, vegm_data = clm_clipper.clip_land_cover(lat_lon_array=latlon_formatted,
                                                                 land_cover_file=os.path.join(conus.local_path,
                                                                                              conus.optional_files.get(
                                                                                                  'LAND_COVER')))

        clm_clipper.write_land_cover(vegm_data, os.path.join(out_dir, f'{out_name}_vegm.dat'))

    # Step 5. Generate TCL File
    if write_tcl == 1:
        build_tcl(os.path.join(out_dir, f'{out_name}.tcl'),
                  parkinglot_template,
                  out_name,
                  os.path.join(out_dir, f'{Path(conus.required_files.get("SLOPE_X")).stem}_clip.pfb'),
                  os.path.join(out_dir, f'{out_name}.pfsol'),
                  os.path.join(out_dir, 'pme.pfb'), end_time=10, batches=batches,
                  p=2, q=1, r=1, timestep=1, constant=1)


def main():
    # setup logging
    start_date = datetime.utcnow()
    # parse the command line arguments
    cmd_line_args = sys.argv
    args = parse_args(cmd_line_args[1:])
    logging.basicConfig(filename=os.path.join(args.out_dir, 'subset_conus.log'), filemode='w', level=logging.INFO)
    logging.info(f'start process at {start_date} from command {" ".join(cmd_line_args[:])}')
    subset_conus(input_path=args.input_path, shapefile=args.shapefile, conus_version=args.conus_version,
                 conus_files=args.conus_files, out_dir=args.out_dir, out_name=args.out_name, clip_clm=args.clip_clm,
                 write_tcl=args.write_tcl, padding=args.padding, attribute_ids=args.attribute_ids,
                 attribute_name=args.attribute_name, write_tifs=args.write_tifs, manifest_file=args.manifest_file)

    end_date = datetime.utcnow()
    logging.info(f'completed process at {end_date} for a runtime of {end_date - start_date}')


if __name__ == '__main__':
    main()
