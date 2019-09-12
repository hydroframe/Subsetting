#!/usr/local/bin/python3

import numpy as np
import gdal, ogr, osr
from pyproj import Proj, transform
import pfio
import matplotlib.pyplot as plt
import os
import gzip
import sys

def project_array(grid_coor, inProj, outProj):
    fx, fy = transform(inProj, outProj, grid_coor[:,1], grid_coor[:,0])
    # Re-create (n,2) coordinates
    return fx, fy

def createTif(fname, res_arr, geom, prj, dtype=gdal.GDT_Float32, ndata=-99):
	driver = gdal.GetDriverByName('GTiff')
	dataset = driver.Create(
			fname,
			res_arr.shape[2],
			res_arr.shape[1],
			res_arr.shape[0],
			dtype,
			options = ['COMPRESS=LZW'])
	dataset.SetGeoTransform(geom)
	dataset.SetProjection(prj)
	for i in range(res_arr.shape[0]):
		dataset.GetRasterBand(i+1).WriteArray(res_arr[i,:,:])
		dataset.GetRasterBand(i+1).SetNoDataValue(ndata)
	dataset.FlushCache()

in_file = sys.argv[1]

in_file_ext = os.path.splitext(in_file)[-1]

if in_file_ext == '.pfb':
	in_arr = pfio.pfread(in_file)
	nz,ny,nx = in_arr.shape
elif in_file_ext == '.txt':
	in_arr = np.loadtxt(in_file,skiprows=1)
	with open(in_file,'r') as fi:
		first_line = fi.readline().strip()
	nx,ny,nz = [int(x) for x in first_line.split()]
	in_arr = np.reshape(in_arr,(nz,ny,nx))
	in_arr = in_arr[:,::-1,:]
else:
	print('does not support this file format..exit')
	sys.exit()

#grid_coor_file = 'Grid_Centers_Short_Deg.format.txt'
grid_coor_file = gzip.GzipFile('grid_coor.npy.gz', "r")
grid_coor = np.load(grid_coor_file)

mask_file = '../CONUS1_inputs/conus_1km_PFmask2.tif'

#open mask file
ds_mask = gdal.Open(mask_file)
mask_arr = ds_mask.ReadAsArray()

##get geo transform
geom = ds_mask.GetGeoTransform()

#in projection
inProj = Proj(init='epsg:4326')


##get mask file projection as outProj
outSRS_converter = osr.SpatialReference()  # makes an empty spatial ref object
outSRS_converter.ImportFromWkt(ds_mask.GetProjection())  # populates the spatial ref object with our WKT SRS
outSRS_forPyProj = outSRS_converter.ExportToProj4()  # Exports an SRS ref as a Proj4 string usable by PyProj
outProj = Proj(outSRS_forPyProj)

##create a zeros array
res_arr = np.ones((nz,mask_arr.shape[0],mask_arr.shape[1]))*-99.

#loop through grid coordinates
out_x, out_y = project_array(grid_coor, inProj, outProj)
out_x = out_x.reshape(ny,nx)
out_y = out_y.reshape(ny,nx)
out_y = out_y[::-1,:]
out_x = out_x[::-1,:]

arr_x = ((out_x - geom[0])/geom[1]).astype(np.int)
arr_y = ((geom[3]-out_y)/geom[1]).astype(np.int)

res_arr[:,arr_y,arr_x] = in_arr

#write tif file
fname = ''.join(os.path.splitext(in_file)[:-1])+'.tif'
if os.path.isfile(fname):
	os.remove(fname)

createTif(fname, res_arr, geom, ds_mask.GetProjection())






