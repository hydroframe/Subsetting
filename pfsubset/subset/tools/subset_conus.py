"""Convenience wrapper for clipping inputs for the CONUS1 or CONUS2 models

Everything here can be customized using the classes in the subset package
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from datetime import datetime
import numpy as np
import shapefile
from pfsubset.subset.utils.arguments import is_valid_path, is_positive_integer, is_valid_file
from pfsubset.subset.clipper import MaskClipper
from pfsubset.subset.domain import Conus
from pfsubset.subset.rasterizer import ShapefileRasterizer
import pfsubset.subset.tools.bulk_clipper as bulk_clipper
from pfsubset.subset.clipper import ClmClipper
from pfsubset.subset.data import parkinglot_template, conus_manifest
from parflow.tools.builders import SolidFileBuilder
from pfsubset.subset import TIF_NO_DATA_VALUE_OUT
from pfsubset.subset.data import parking_lot_template
from pfsubset.subset.utils import huc2shape
from pfsubset.subset.mask import SubsetMask

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
    exclusive_group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument("--input_path", "-i", dest="input_path", required=True,
                        type=lambda x: is_valid_path(parser, x),
                        help="the input path to the shapefile file set")

    exclusive_group.add_argument("--shapefile", "-s", dest="shapefile", required=False, default=None,
                                help="the name of the shapefile file set")

    exclusive_group.add_argument("--watershed", "-w", dest="hucs", required=False, nargs="+", default=None,
                                 help="one or more HUC 12 values that identify a contiguous watershed")

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

    parser.add_argument("--run_script", "-r", dest="run_script", required=False,
                        action='store_true',
                        help="generate the run script for this subset")

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


def subset_conus(input_path=None, shapefile=None, subset_tif=None, mask_value=None, conus_version=1, conus_files='.',
                 out_dir='.', out_name=None, clip_clm=False, run_script=False, padding=(0, 0, 0, 0),
                 attribute_name='OBJECTID', attribute_ids=None, write_tifs=False, manifest_file=conus_manifest):
    """subset a conus domain inputs for running a regional model

    Parameters
    ----------
    input_path : str, optional
        path to input shapefile parts to use as mask. Only applicable if shapefile is also provided.
    shapefile : str, optional
        name of shapefile to use as mask. Required if subset_tif is not provided.
    subset_tif : str, optional
        path to tiff file containing mask. Required if shapefile is not provided.
    mask_value : int, optional
        integer value specifying the data value in the tiff file to consider as the masking value.
        If None, all +ve values are considered as the masking value
        Only applicable if subset_tif is provided.
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
    run_script : int, optional
        whether or not to build and return a Run object for the subset (default no)
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
    run_script : pfsubset.tools.Run
        The Run object which can be used to execute the ParFlow model subset that was created by subset_conus
    """
    assert any((shapefile, subset_tif)) and not all((shapefile, subset_tif)), \
        'Specify either a shapefile or a subset_tif file.'

    if out_name is None:
        out_name = shapefile or os.path.splitext(os.path.basename(subset_tif))[0]
    conus = Conus(version=conus_version, local_path=conus_files, manifest_file=manifest_file)

    if attribute_ids is None:
        attribute_ids = [1]

    # Step 1, rasterize shapefile
    if shapefile is not None:
        rasterizer = ShapefileRasterizer(input_path, shapefile, reference_dataset=conus.get_domain_tif(),
                                         no_data=TIF_NO_DATA_VALUE_OUT, output_path=out_dir, )
        mask_array = rasterizer.rasterize_shapefile_to_disk(out_name=f'{out_name}_raster_from_shapefile.tif',
                                               padding=padding,
                                               attribute_name=attribute_name,
                                               attribute_ids=attribute_ids)

        subset_mask = rasterizer.subset_mask
    else:
        subset_mask = SubsetMask(tif_file=subset_tif, mask_value=mask_value)
        mask_array = subset_mask.mask_array

    # Step 2, Generate solid file
    clip = MaskClipper(subset_mask, no_data_threshold=-1)

    sfb = SolidFileBuilder(top=3, bottom=6, side=1).mask(clip.subset(mask_array, crop_inner=0)[0][0, :, :])
    # identify the unique patch ID's assigned to the solid file

    #TODO: get patch defs from a class
    sfb.top_ids(clip.subset(conus.get_domain_mask())[0][0, :, :])
    sfb.side_ids(clip.subset(conus.get_border_mask())[0][0, :, :])

    top_patchIDs = np.unique(clip.subset(conus.get_domain_mask())[0][0, :, :])
    side_patchIDs = np.unique(clip.subset(conus.get_border_mask())[0][0, :, :])
    side_patchIDs[side_patchIDs == 0] = 2
    botom_patchIDs = [6]
    patch_ids = np.unique(np.concatenate((top_patchIDs, side_patchIDs, botom_patchIDs)))

    sfb = sfb.write(os.path.join(out_dir, f'{out_name}.pfsol'), cellsize=1000, vtk=True)

    # Step 3. Clip all the domain data inputs
    bulk_clipper.clip_inputs(clip,
                             [os.path.join(conus.local_path, value) for key, value in conus.required_files.items()
                              if key not in ['DOMAIN_MASK', 'CHANNELS']],
                             out_dir=out_dir, tif_outs=write_tifs)

    # Step 4. Clip CLM inputs
    if clip_clm == 1:
        clm_clipper = ClmClipper(subset_mask.get_bbox())
        latlon_formatted, latlon_data = clm_clipper.clip_latlon(os.path.join(conus.local_path,
                                                                             conus.optional_files.get('LAT_LON')))

        clm_clipper.write_lat_lon(latlon_formatted, os.path.join(out_dir, f'{out_name}_latlon.sa'),
                                  x=latlon_data.shape[2], y=latlon_data.shape[1], z=latlon_data.shape[0])

        land_cover_data, vegm_data = clm_clipper.clip_land_cover(lat_lon_array=latlon_formatted,
                                                                 land_cover_file=os.path.join(conus.local_path,
                                                                                              conus.optional_files.get(
                                                                                                  'LAND_COVER')))

        clm_clipper.write_land_cover(vegm_data, os.path.join(out_dir, f'{out_name}_vegm.dat'))

    # Step 5. Generate Run Script

    if run_script == 1:
        # the Run object reads sys.argv, and this is problematic because they share a common flag -r
        sys.argv = ['Run']
        slopex_file = os.path.join(out_dir, f'{Path(conus.required_files.get("SLOPE_X")).stem}_clip.pfb')
        slopey_file = os.path.join(out_dir, f'{Path(conus.required_files.get("SLOPE_Y")).stem}_clip.pfb')
        solid_file = os.path.join(out_dir, f'{out_name}.pfsol')
        bbox = subset_mask.get_bbox()
        extents = bbox.get_padded_extents()

        NX = int(extents[3] - extents[2])
        NY = int(extents[1] - extents[0])

        out_name = f'{out_name}.conus{conus_version}.parking_lot'
        # TODO: associate model templates with models and versions, provide method to override boundary conditions
        run_script = parking_lot_template.get_parking_lot_model(out_name, slopex_file, slopey_file, solid_file, NX, NY)
        patch_names = [conus.get_patch_name(patch_id) for patch_id in patch_ids[patch_ids>TIF_NO_DATA_VALUE_OUT]]
        run_script.Geom.domain.Patches = ' '.join(patch_names)

        # convert patch ID's to patch names for run script
        for patch in patch_names:
            bc = parking_lot_template.get_parking_lot_model_boundary(patch)
            for k, v in bc.items():
                # assign patch Boundary Conditions as defined by CONUS Model 1 or 2
                run_script.Patch.pfset(key=f'{patch}.BCPressure.{k}', value=v)

        if conus_version == 1:
            # CONUS1 doesn't seem to work well with OverlandKinematic
            run_script.Patch.top.BCPressure.Type = 'OverlandFlow'

        run_script.validate()
        run_script.write(file_name=os.path.join(out_dir, out_name), file_format='pfidb')
        run_script.write(file_name=os.path.join(out_dir, out_name), file_format='yaml')
        run_script.write(file_name=os.path.join(out_dir, out_name), file_format='json')
        return run_script


def main():
    # setup logging
    start_date = datetime.utcnow()
    # parse the command line arguments
    cmd_line_args = sys.argv
    args = parse_args(cmd_line_args[1:])
    logging.basicConfig(filename=os.path.join(args.out_dir, 'subset_conus.log'), filemode='w', level=logging.INFO)
    logging.info(f'start process at {start_date} from command {" ".join(cmd_line_args[:])}')
    shape = args.shapefile
    ids = args.attribute_ids
    out_name = args.out_name
    attribute_name = args.attribute_name
    if args.hucs is not None:
        if args.out_name is None:
            out_name = 'watershedClip'

        # query the geometries and build the watershed shape object
        watershed = huc2shape.build_shape_object(args.hucs)

        shape = out_name
        attribute_name = 'ID'
        # write the watershed shape object to *.shp
 
        watershed.write_shapefile((Path(args.input_path) / (shape  + '.shp')).as_posix())

        # collect all the ID values in the generated shapefile
        ids = []
        with shapefile.Reader((Path(args.input_path) / (shape  + '.shp')).as_posix()) as shp:
            for record in shp.records():
                ids.append(str(record[0]))

    subset_conus(input_path=args.input_path, shapefile=shape, conus_version=args.conus_version,
                 conus_files=args.conus_files, out_dir=args.out_dir, out_name=out_name, clip_clm=args.clip_clm,
                 run_script=args.run_script, padding=args.padding, attribute_ids=ids,
                 attribute_name=attribute_name, write_tifs=args.write_tifs, manifest_file=args.manifest_file)

    end_date = datetime.utcnow()
    logging.info(f'completed process at {end_date} for a runtime of {end_date - start_date}')


if __name__ == '__main__':
    main()
