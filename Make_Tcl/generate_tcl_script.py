#!/usr/bin/env python3

import numpy as np
import os
import sys
import argparse
from os.path import join


def generate_tcl(pfsol_file, slope_name, temp_file, run_name,
                 asc_path, topo_path):

    # Check if pfsol, pfb file exists

    top_asc_file = join(asc_path, "Top_Border.asc")
    bottom_asc_file = join(asc_path, "Bottom_Border.asc")
    front_asc_file = join(asc_path, "Front_Border.asc")
    right_asc_file = join(asc_path, "Right_Border.asc")
    left_asc_file = join(asc_path, "Left_Border.asc")
    back_asc_file = join(asc_path, "Back_Border.asc")

    pfb_slopex = f'{slope_name}x.pfb'
    pfb_slopey = f'{slope_name}y.pfb'

    for req_file in [pfsol_file, pfb_slopex, pfb_slopey, top_asc_file,
                     bottom_asc_file, front_asc_file, right_asc_file,
                     left_asc_file, back_asc_file]:
        if not os.path.isfile(req_file):
            print('missing '+req_file+'...please check in '
                  'Domain or Topography generation')
            sys.exit()

    # Check if template file exits
    if not os.path.isfile(temp_file):
        print('template file is missing')
        sys.exit()

    # If all the file exits, copy them to current folder
    for req_file in [pfsol_file, pfb_slopex, pfb_slopey, top_asc_file,
                     bottom_asc_file, front_asc_file, right_asc_file,
                     left_asc_file, back_asc_file]:
        os.system('cp '+req_file+' .')

    # Get dimensions of the files
    with open(top_asc_file, 'r') as fi:
        head = [next(fi) for x in range(6)]

    for line in head:
        line = line.strip()
        if 'ncols' in line:
            nx = int(line.split()[1])
        if 'nrows' in line:
            ny = int(line.split()[1])

    patches_array = np.array([]).reshape(nx*ny, 0)

    # Get patches values
    for patch_file in [top_asc_file, bottom_asc_file, front_asc_file,
                       right_asc_file, left_asc_file, back_asc_file]:
        temp_arr = np.loadtxt(patch_file, skiprows=6)
        patches_array = np.hstack([patches_array, temp_arr.reshape(-1, 1)])

    list_patches = np.unique(patches_array).tolist()
    patches_dictionary = {0.: 'land',
                          1.: 'ocean',
                          3.: 'top',
                          4.: 'lake',
                          5.: 'sink',
                          6.: 'bottom'}

    patches_string = ''
    for p in list_patches:
        patches_string += patches_dictionary[p]+' '

    # read current template file
    with open(temp_file, 'r') as fi:
        tcl_content = fi.read()

    tcl_content = tcl_content.split('\n')

    new_content = []
    for line in tcl_content:
        # change runname
        if 'set runname' in line:
            line = line.replace('CONUS2_Parkinglot', run_name)
        # set computational grid nx
        if 'pfset ComputationalGrid.NX' in line:
            line = line.replace('4442', str(nx))
        # set computational grid ny
        if 'pfset ComputationalGrid.NY' in line:
            line = line.replace('3256', str(ny))
        # set name of the solid file
        if 'pfset GeomInput.domaininput.FileName' in line:
            line = line.replace('CONUS2.0_test1.pfsol',
                                os.path.basename(pfsol_file))
        # set available patches
        if 'pfset Geom.domain.Patches' in line:
            line = 'pfset Geom.domain.Patches "'+patches_string+'"'
        # set patch name for BC
        if 'pfset BCPressure.PatchNames' in line:
            line = 'pfset BCPressure.PatchNames "'+patches_string+'"'
        # set slope .pfb file as input
        if 'pfset TopoSlopesY.Type' in line:
            line = 'pfset TopoSlopesY.Type PFBFile'
        # set slope .pfb file as input
        if 'pfset TopoSlopesX.Type' in line:
            line = 'pfset TopoSlopesX.Type PFBFile'
        # set name of the pfb file
        if 'pfset TopoSlopesX.FileName' in line:
            line = 'pfset TopoSlopesX.FileName '+os.path.basename(pfb_slopex)
        # set name of the pfb file
        if 'pfset TopoSlopesY.FileName' in line:
            line = 'pfset TopoSlopesY.FileName '+os.path.basename(pfb_slopey)
        if 'pfdist CONUS.slopex.pfb' in line:
            line = line.replace('CONUS.slopex.pfb',
                                os.path.basename(pfb_slopex))
        if 'pfdist CONUS.slopey.pfb' in line:
            line = line.replace('CONUS.slopey.pfb',
                                os.path.basename(pfb_slopey))
        new_content.append(line+'\n')

    new_content.append('pfundist '+os.path.basename(pfb_slopex)+'\n')
    new_content.append('pfundist '+os.path.basename(pfb_slopey)+'\n')

    with open(run_name+'.tcl', 'w') as fo:
        for line in new_content:
            fo.write(line)

    # Remove all the previous simulation
    os.system('rm -rf '+run_name+'.out.*')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate TCL script for '
                                                 'executing parFlow '
                                                 'simulations')

    parser.add_argument('-pfsol', '--pfsol_file', type=str, required=True,
                        help='ParFlow solid file, e.g. subset.pfsol')
    parser.add_argument('-s', '--slope', type=str, required=True,
                        help='ParFlow slope file without extension or direction, '
                             'e.g. output_slope instead of output_slopex.pfb')
    parser.add_argument('-t', '--template', type=str, required=True,
                        help='ParFlow template file')
    parser.add_argument('-r', '--run_name', type=str, required=False,
                        default="",
                        help='ParFlow run name')
    parser.add_argument('-a', '--asc_path', type=str, required=False,
                        default="../Domain",
                        help='Path to asc files, e.g. Back_Border.asc')
    parser.add_argument('-topo', '--topography_path', type=str, required=False,
                        default="../Topography",
                        help='Path to topography files, e.g. *.pfb')

    args = parser.parse_args()

    # run default run name
    if args.run_name == "":
        args.run_name = f'{args.mask}_{args.slope_name}'

    generate_tcl(args.pfsol_file,
                 args.slope,
                 args.template,
                 args.run_name,
                 args.asc_path,
                 args.topography_path)
