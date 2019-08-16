#!/usr/local/bin/python3

import gdal
import ogr
import osr
import numpy as np
import os
import sys
import argparse
from pyproj import Proj, transform
from os.path import join, abspath, dirname


def pixzone2latlon(xul, yul, dx, dy, x0, y0):
    lat = yul - dy*y0
    lon = xul + dx*x0
    return lat, lon


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
                           options=[f"ATTRIBUTE={attribute_name}",
                                    "ALL_TOUCHED=TRUE"]) != 0:
        raise Exception("error rasterizing layer: %s" % shp_layer)
    else:
        target_ds.FlushCache()
        return 0


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-shp_file', type=str, required=True,
                        help='input shapefile')
    parser.add_argument('-id', type=int, nargs='+', required=True,
                        help='id of the selected watershed')
    parser.add_argument('-att', type=str, required=False,
                        default='OBJECTID',
                        help='Column name of the shape attribute to use '
                             'during rasterize')
    parser.add_argument('-pfmask', '--pf_conus_mask_1km', type=str,
                        default='conus_1km_PFmask2.tif',
                        help='ParFlow CONUS 1km mask')
    parser.add_argument('-out_name', type=str, default='UC.latlon.txt',
                        help='Name of output coordinate file')

    args = parser.parse_args()

    basin_ids = args.id
    out_file = args.out_name
    conus_pf_1k_mask = args.pf_conus_mask_1km
    region_shp = args.shp_file
    conus_pf_1k_tifs = [conus_pf_1k_mask]

    outdir = join(abspath(dirname(out_file)))

    region_raster = f'{outdir}/Regions.tif'

    region_shps = [region_shp.split('.')[0]+x for x in ['.shp',
                                                        '.dbf',
                                                        '.prj',
                                                        '.shx',
                                                        '.sbx',
                                                        '.sbn']]

    # check if file exits, if not we need to login to avra and download.
    # This part requires icommand authorization
    if any([not os.path.isfile(x) for x in conus_pf_1k_tifs]):
        avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/domain/'
        avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/' + \
                        'SteadyState_Final/Subdomain_Extraction/' + \
                        'Shape_Files/Regions_shp/'
        print(conus_pf_1k_mask+' does not exits...downloading from avra')
        auth = os.system('iinit')
        if auth != 0:
            print('Authentication failed...exit')
            sys.exit()

        for tif_file in conus_pf_1k_tifs:
            os.system('iget -K '+avra_path_tif+tif_file+' .')

    if not os.path.isfile(region_shp):
        print(region_shp+' does not exits...downloading from avra')
        auth = os.system('iinit')
        if auth != 0:
            print('Authentication failed...exit')
            sys.exit()

        for shp_component_file in region_shps:
            os.system('iget -K '+avra_path_shp+shp_component_file+' .')

    # read domain raster
    ds_ref = gdal.Open(conus_pf_1k_mask)
    arr_ref = ds_ref.ReadAsArray()
    geom_ref = ds_ref.GetGeoTransform()

    # makes an empty spatial ref object
    inSRS_converter = osr.SpatialReference()

    # populates the spatial ref object with our WKT SRS
    inSRS_converter.ImportFromWkt(ds_ref.GetProjection())

    # Exports an SRS ref as a Proj4 string usable by PyProj
    inSRS_forPyProj = inSRS_converter.ExportToProj4()
    inProj = Proj(inSRS_forPyProj)
    outProj = Proj(init='epsg:4326')

    # rasterize region shapefile
    if os.path.isfile(region_raster):
        os.remove(region_raster)

    rasterize(region_raster, region_shp, ds_ref, attribute_name=args.att)

    shp_raster_arr = gdal.Open(region_raster).ReadAsArray()

    yy, xx = np.where(np.isin(shp_raster_arr, basin_ids))

    # find the extends
    new_arr = shp_raster_arr[min(yy):max(yy)+1, min(xx):max(xx)+1]

    # add n extra grid cells to every direction
    n = int(max(new_arr.shape)*0.02)

    min_x = max(min(xx)-n, 0)
    min_y = max(min(yy)-n, 0)
    max_x = max(xx)+n+1
    max_y = max(yy)+n+1

    min_lat, min_lon = pixzone2latlon(np.round(geom_ref[0], 5),
                                      np.round(geom_ref[3], 5),
                                      geom_ref[1],
                                      geom_ref[1],
                                      min_x,
                                      max_y)

    max_lat, max_lon = pixzone2latlon(np.round(geom_ref[0], 5),
                                      np.round(geom_ref[3], 5),
                                      geom_ref[1],
                                      geom_ref[1],
                                      max_x,
                                      min_y)

    lats = np.arange(min_lat, max_lat, geom_ref[1])
    lons = np.arange(min_lon, max_lon, geom_ref[1])

    # write lat long file
    with open(out_file, 'w') as fo:
        for lat in lats:
            for lon in lons:
                results = ' '.join([str(np.round(x, 5)) for x in
                                    transform(inProj,
                                              outProj,
                                              lon,
                                              lat)][::-1])
                fo.write(results+'\n')

