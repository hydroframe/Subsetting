#!/usr/bin/env python3

import numpy as np
import gdal
import osr
from pyproj import Proj, transform
import argparse


def create_vegm_latlon(igbp_file,
                       latlon_file,
                       out_file='drv_vegm.UC.dat',
                       nx=637,
                       ny=903,
                       sand=0.16,
                       clay=0.26,
                       color=2):

    # open land cover file
    lc_ds = gdal.Open(igbp_file)

    # get land cover file projection as outProj
    # makes an empty spatial ref object
    outSRS_converter = osr.SpatialReference()

    # populates the spatial ref object with our WKT SRS
    outSRS_converter.ImportFromWkt(lc_ds.GetProjection())

    # Exports an SRS ref as a Proj4 string usable by PyProj
    outSRS_forPyProj = outSRS_converter.ExportToProj4()
    outProj = Proj(outSRS_forPyProj)

    # get land cover array
    lc_arr = lc_ds.ReadAsArray()

    # get land cover geo transform
    xul, xres, _, yul, _, yres = lc_ds.GetGeoTransform()

    # open lat long file
    latlon_arr = np.loadtxt(latlon_file)
    inProj = Proj(init='epsg:4326')

    # get value of land cover for each coordinate
    npoints = latlon_arr.shape[0]

    # make output matrix
    output = np.zeros((npoints, 25))
    output[:, 4] = sand
    output[:, 5] = clay
    output[:, 6] = color

    for ii in range(npoints):
        in_lat = latlon_arr[ii, 0]
        in_lon = latlon_arr[ii, 1]
        out_lon, out_lat = transform(inProj, outProj, in_lon, in_lat)
        lc_val = lc_arr[int(float(out_lat - yul)/yres),
                        int(float(out_lon - xul)/xres)]
        if (ii+1) % nx == 0:
            x = nx
        else:
            x = (ii+1) % nx
        y = ii // nx + 1

        # start writing to output matrix
        output[ii, 0:4] = x, y, in_lat, in_lon
        output[ii, lc_val+6] = 1

    heading = "x y lat lon sand clay color fractional coverage of grid, " + \
              "by vegetation class (Must/Should Add to 1.0)"

    col_names = ['',  '',  '(Deg)', '(Deg)', '(%/100)', '', 'index', '1', '2',
                 '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
                 '14', '15', '16', '17', '18']

    # write to out file
    with open(out_file, 'w') as fo:
        fo.write(heading+'\n')
        fo.write(' '.join(col_names)+'\n')
        np.savetxt(fo, output, fmt=['%d']*2+['%.5f']*2+['%.2f']*2+['%d']*19)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create vegm Lat Lon script')
    parser.add_argument('-input_igbp', type=str, required=True,
                        help='input igbp file')
    parser.add_argument('-input_latlon', type=str, required=True,
                        help='input latlon file')
    parser.add_argument('-out_name', type=str, default='drv_vegm.UC.dat',
                        help='input igbp file')
    parser.add_argument('-nx', type=int, default=637,
                        help='nx')
    parser.add_argument('-ny', type=int, default=903,
                        help='ny')
    parser.add_argument('-sand', type=float, default=0.16,
                        help='sand param')
    parser.add_argument('-clay', type=float, default=0.26,
                        help='clay param')
    parser.add_argument('-color', type=int, default=2,
                        help='color param')

    # parse arguments
    args = parser.parse_args()

    create_vegm_latlon(args.input_igbp,
                       args.input_latlon,
                       args.out_name,
                       args.nx,
                       args.ny,
                       args.sand,
                       args.clay,
                       args.color)

