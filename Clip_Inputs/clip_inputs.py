#!/usr/bin/env python3

import gdal
import ogr
import numpy as np
import pandas as pd
import argparse
import os
import sys
import pfio
from Define_Watershed import DelinWatershed


def latlon2pixzone(xul, yul, dx, dy, lat0, lon0):
    rl = abs((yul - lat0)/dy)
    cu = abs((lon0 - xul)/dx)
    return int(round(rl)), int(round(cu))


def rasterize(out_raster, in_shape, ds_ref,
              dtype=gdal.GDT_Int32,
              ndata=-99,
              attribute_name='OBJECTID'):

    # target raster file
    geom_ref = ds_ref.GetGeoTransform()
    target_ds = gdal.GetDriverByName('GTiff').Create(out_raster,
                                                     ds_ref.RasterXSize,
                                                     ds_ref.RasterYSize,
                                                     1, dtype)
    target_ds.SetProjection(ds_ref.GetProjection())
    target_ds.SetGeoTransform(geom_ref)
    target_ds.GetRasterBand(1).SetNoDataValue(ndata)

    # shapefile
    shp_source = ogr.Open(in_shape)
    shp_layer = shp_source.GetLayer()

    # Rasterize layer
    if gdal.RasterizeLayer(target_ds, [1],
                           shp_layer,
                           options=[f"ATTRIBUTE={attribute_name}"]) != 0:
        raise Exception("error rasterizing layer: %s" % shp_layer)
    else:
        target_ds.FlushCache()
        return 0


def read_from_file(infile):
    # get extension
    ext = os.path.splitext(os.path.basename(infile))[1]
    if ext in ['.tif', '.tiff']:
        res_arr = gdal.Open(infile).ReadAsArray()
    # parflow ascii file
    elif ext == '.sa':
        with open(infile, 'r') as fi:
            header = fi.readline()
        nx, ny, nz = [int(x) for x in header.strip().split(' ')]
        arr = pd.read_csv(infile, skiprows=1, header=None).values
        res_arr = np.reshape(arr, (nz, ny, nx))[:, ::-1, :]
    # parflow binary file
    elif ext == '.pfb':
        res_arr = pfio.pfread(infile)
    else:
        print('can not read file type '+ext)
        sys.exit()
    return res_arr


def subset(arr, mask_arr, ds_ref, crop_to_domain, ndata=0):

    arr1 = arr.copy()

    # create new geom
    old_geom = ds_ref.GetGeoTransform()

    # find new up left index
    yy, xx = np.where(mask_arr == 1)
    len_y = max(yy)-min(yy)+1
    len_x = max(xx)-min(xx)+1

    # add grid cell to make dimensions as multiple of 32 (nicer PxQxR)
    new_len_y = ((len_y//32)+1)*32
    n1 = (new_len_y-len_y)//2
    n2 = new_len_y-len_y-n1
    new_len_x = ((len_x//32)+1)*32
    n3 = (new_len_x-len_x)//2
    n4 = new_len_x-len_x-n3

    # create new geom
    new_x = old_geom[0] + old_geom[1]*(min(xx)+1)
    new_y = old_geom[3] + old_geom[5]*(min(yy)+1)
    new_geom = (new_x, old_geom[1], old_geom[2], new_y,
                old_geom[4], old_geom[5])

    # start subsetting
    if len(arr.shape) == 2:
        if crop_to_domain:
            arr1[mask_arr != 1] = ndata
            return_arr = np.zeros((new_len_y, new_len_x))
            return_arr[n1:-n2, n3:-n4] = arr1[min(yy):max(yy) + 1,
                                              min(xx):max(xx)+1]
        else:
            return_arr = arr1[min(yy)-n1:max(yy)+n2+1,
                              min(xx)-n3:max(xx)+n4+1]
        return_arr = return_arr[np.newaxis, ...]
    else:
        if crop_to_domain:
            arr1[:, mask_arr != 1] = ndata
            return_arr = np.zeros((arr1.shape[0], new_len_y, new_len_x))
            return_arr[:, n1:-n2, n3:-n4] = arr1[:, min(yy):max(yy)+1,
                                                 min(xx):max(xx)+1]
        else:
            return_arr = arr1[:, min(yy)-n1:max(yy)+n2+1,
                              min(xx)-n3:max(xx)+n4+1]
    return return_arr, new_geom


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a clipped input '
                                                 'ParFlow binary files')
    parser.add_argument('-i', '--input_file', type=str,
                        help='input file for subsetting')
    parser.add_argument('--crop_to_domain', type=int,
                        help='crop to domain (i.e. value outside of the '
                             'domain will be assign as nodata -- optional). '
                             'Default is 1')
    parser.add_argument('-x0', type=int, help='x0 (optional).Default is 0')
    parser.add_argument('-y0', type=int, help='y0 (optional).Default is 0')
    parser.add_argument('-z0', type=int, help='z0 (optional).Default is 0')
    subparsers = parser.add_subparsers(dest='type',
                                       help='subset using three options:')

    # group 1: using shapefile
    parser_a = subparsers.add_parser('shapefile',
                                     help='subset using shapefile and the '
                                          'selected id of watershed')
    parser_a.add_argument('-shp_file', type=str, help='input shapefile')
    parser_a.add_argument('-id', type=int, help='id of the selected watershed')
    parser_a.add_argument('-out_name', type=str,
                          help='Name of output solidfile (optional)')
    parser_a.add_argument('-dx', type=int,
                          help='Spatial resolution of solidfile (optional). '
                               'Default is 1000')
    parser_a.add_argument('-dz', type=int,
                          help='Lateral resolution of solidfile (optional). '
                               'Default is 1000')
    parser_a.add_argument('-printmask', type=int,
                          help='Print mask (optional). Default is 0')
    parser_a.add_argument('-att', type=str, required=False, default='OBJECTID',
                          help='Column name of the shape attribute to use '
                               'during rasterize')
    # parser_a.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
    # parser_a.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

    # group 2: using mask file
    parser_b = subparsers.add_parser('mask', help='subset using a mask file')
    parser_b.add_argument('-mask_file', type=str, help='input mask file')
    parser_b.add_argument('-out_name', type=str,
                          help='Name of output solidfile (optional)')
    parser_b.add_argument('-dx', type=int,
                          help='Spatial resolution of solidfile (optional). '
                               'Default is 1000')
    parser_b.add_argument('-dz', type=int,
                          help='Lateral resolution of solidfile (optional). '
                               'Default is 1000')
    parser_b.add_argument('-printmask', type=int,
                          help='print mask (optional). Default is 0')
    # parser_b.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
    # parser_b.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

    # group 3: using custom watershed
    parser_c = subparsers.add_parser('define_watershed',
                                     help='subset using a newly created '
                                          'watershed')
    parser_c.add_argument('-dir_file', type=str, help='input direction file',)
    parser_c.add_argument('-outlet_file', type=str,
                          help='file contains coordinates of outlet points')
    parser_c.add_argument('-out_name', type=str,
                          help='name of output solidfile (required)')
    parser_c.add_argument('-dx', type=int,
                          help='Spatial resolution of solidfile (optional). '
                               'Default is 1000')
    parser_c.add_argument('-dz', type=int,
                          help='Lateral resolution of solidfile (optional). '
                               'Default is 1000')
    parser_c.add_argument('-printmask', type=int,
                          help='print mask (optional). Default is 0')
    # parser_c.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
    # parser_c.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

    # required raster files
    conus_pf_1k_mask = 'conus_1km_PFmask2.tif'

    avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/domain/'

    # Check if file exits, if not we need to login to avra and download.
    # This part requires icommand authorization
    if not os.path.isfile(conus_pf_1k_mask):
        print(conus_pf_1k_mask+' does not exits...downloading from avra')
        auth = os.system('iinit')
        if auth != 0:
            print('Authentication failed...exit')
            sys.exit()

        os.system('iget -K '+avra_path_tif+conus_pf_1k_mask+' .')

    # read domain raster
    ds_ref = gdal.Open(conus_pf_1k_mask)
    arr_ref = ds_ref.ReadAsArray()

    # parsing arguments
    args = parser.parse_args()

    # read input file
    try:
        infile = args.i
    except Exception:
        infile = args.input_file

    if not os.path.isfile(infile):
        print(infile+' does not exits...exitting')
        sys.exit()

    file_ext = os.path.splitext(os.path.basename(infile))[1]
    if file_ext == '.tif':
        ds_in = gdal.Open(infile)

        # check if direction file has the same projection and
        # extent with the domain mask file
        if any([ds_ref.GetProjection() != ds_in.GetProjection(),
                sorted(ds_ref.GetGeoTransform()) !=
                sorted(ds_in.GetGeoTransform())]):
            print('input and domain do not match...exit')
            sys.exit()

    arr_in = read_from_file(infile)

    # deal with optional arguments
    if args.crop_to_domain is None:
        crop_to_domain = 1
    else:
        crop_to_domain = args.crop_to_domain

    if args.dx is None:
        dx = 1000.
    else:
        dx = args.dx

    if args.dz is None:
        dz = 1000.
    else:
        dz = args.dz

    if args.printmask is None:
        printmask = 0
    else:
        printmask = 1

    if args.x0 is None:
        x0 = 0.
    else:
        x0 = args.x0

    if args.y0 is None:
        y0 = 0.
    else:
        y0 = args.y0

    if args.z0 is None:
        z0 = 0.
    else:
        z0 = args.z0

    # main arguments
    if args.type == 'shapefile':
        basin_id = args.id
        region_shp = args.shp_file

        # check if shapefile exits locally
        avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/' + \
                        'SteadyState_Final/Subdomain_Extraction/' + \
                        'Shape_Files/Regions_shp/'

        region_shps = [region_shp.split('.')[0]+x for x in
                       ['.shp', '.dbf', '.prj', '.shx', '.sbx', '.sbn']]

        region_raster = 'Regions.tif'

        if not os.path.isfile(region_shp):
            print(region_shp+' does not exits...downloading from avra')
            auth = os.system('iinit')
            if auth != 0:
                print('Authentication failed...exit')
                sys.exit()

            for shp_component_file in region_shps:
                os.system('iget -K '+avra_path_shp+shp_component_file+' .')

        # rasterize region shapefile
        if os.path.isfile(region_raster):
            os.remove(region_raster)

        rasterize(region_raster, region_shp, ds_ref, attribute_name=args.att)

        shp_raster_arr = gdal.Open(region_raster).ReadAsArray()

        # mask array
        mask_arr = (shp_raster_arr == basin_id).astype(np.int)

    elif args.type == 'mask':
        mask_file = args.mask_file
        if not os.path.isfile(mask_file):
            print(mask_file+' does not exits...please create one')
            sys.exit()

        file_ext = os.path.splitext(os.path.basename(mask_file))[1]
        if file_ext == '.tif':
            ds_mask = gdal.Open(mask_file)

            # check if mask file has the same projection and
            # extent with the domain mask file
            if any([ds_ref.GetProjection() != ds_mask.GetProjection(),
                    sorted(ds_ref.GetGeoTransform()) !=
                    sorted(ds_mask.GetGeoTransform())]):
                print('mask and domain do not match...exit')
                sys.exit()

        mask_arr = read_from_file(mask_file)

    elif args.type == 'define_watershed':
        dir_file = args.dir_file
        if not os.path.isfile(dir_file):
            print(dir_file+' does not exits...downloading from avra')
            auth = os.system('iinit')
            if auth != 0:
                print('Authentication failed...exit')
                sys.exit()

            avra_path_direction = '/iplant/home/shared/avra/CONUS2.0/' + \
                                  'Inputs/Topography/Str5Ep0/'
            os.system('iget -K '+avra_path_direction+dir_file+' .')

        file_ext = os.path.splitext(os.path.basename(dir_file))[1]
        if file_ext == '.tif':
            ds_dir = gdal.Open(dir_file)

            # check if direction file has the same projection and
            # extent with the domain mask file
            if any([ds_ref.GetProjection() != ds_dir.GetProjection(),
                    sorted(ds_ref.GetGeoTransform()) !=
                    sorted(ds_dir.GetGeoTransform())]):
                print('direction and domain do not match...exit')
                sys.exit()

        outlet_file = args.outlet_file
        if not os.path.isfile(outlet_file):
            print(outlet_file+' does not exits...please create one')
            sys.exit()

        dir_arr = read_from_file(dir_file)
        queue = np.loadtxt(outlet_file)
        queue = queue.reshape(-1, 2)

        # get the mask array from DelinWatershed function
        mask_arr = DelinWatershed(queue, dir_arr, printflag=True)

    if printmask:
        import matplotlib.pyplot as plt
        plt.imshow(arr_ref+mask_arr*2)
        plt.savefig('mask.png')

    # crop to get a tighter mask
    clip_arr, new_geom = subset(arr_in, mask_arr, ds_ref, crop_to_domain)

    # create clipped outputs
    if args.out_name:
        out_name = args.out_name
        out_name = os.path.splitext(out_name)[0]
        out_pfb = out_name+'.pfb'
    else:
        if args.type == 'shapefile':
            out_name = str(basin_id)
        elif args.type == 'mask':
            out_name = os.path.splitext(os.path.basename(args.mask_file))[0]
        elif args.type == 'define_watershed':
            print('need to specified name of the solid file')
            sys.exit()
        out_pfb = out_name+'.pfb'

    if os.path.isfile(out_pfb):
        os.remove(out_pfb)

    pfio.pfwrite(clip_arr, out_pfb,
                 float(x0), float(y0),
                 float(z0), float(dx),
                 float(dx), float(dz))
