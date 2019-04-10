#!/usr/local/bin/python3

import gdal, ogr
import shapely
import numpy as np
from glob import glob
import os
import sys
import subprocess
import matplotlib.pyplot as plt

def latlon2pixzone(xll, yll, dx, dy, lat0,lat1,lon0,lon1):
	rl = abs((yll - lat0)/dy)
	ru = abs((yll - lat1)/dy)
	cu = abs((lon0 - xll)/dx)
	cl = abs((lon1 - xll)/dx)
	return int(round(rl)), int(round(ru)), int(round(cl)), int(round(cu))

def rasterize(out_raster,in_shape,ds_ref,dtype=gdal.GDT_Int16,ndata=-99):
	#target raster file
	geom_ref = ds_ref.GetGeoTransform()
	target_ds = gdal.GetDriverByName('GTiff').Create(out_raster,
													ds_ref.RasterXSize,
													ds_ref.RasterYSize,
													1, dtype)
	target_ds.SetProjection(ds_ref.GetProjection())
	target_ds.SetGeoTransform(geom_ref)
	target_ds.GetRasterBand(1).SetNoDataValue(ndata)
	#shapefile
	shp_source = ogr.Open(in_shape)
	shp_layer = shp_source.GetLayer()
	#Rasterize layer
	if gdal.RasterizeLayer(target_ds, [1],
							shp_layer,
							options=["ATTRIBUTE=OBJECTID"])  != 0:
		raise Exception("error rasterizing layer: %s" % shp_layer)
	else:
		target_ds.FlushCache()
		return 0

def reshape2d(f,ds_ref):
	arr = np.loadtxt(f,skiprows=1)
	arr = np.reshape(arr,(ds_ref.RasterYSize,ds_ref.RasterXSize),'C')
	return np.flipud(arr)

def reshape3d(f,ds_ref):
	arr = np.loadtxt(f,skiprows=1)
	arr = np.reshape(arr,(5,ds_ref.RasterYSize,ds_ref.RasterXSize),'C')
	return arr[:,::-1,:]

def subset(arr,shp_raster_arr,value,ds_ref, ndata=0):
	arr1 = arr.copy()
	#create new geom
	old_geom = ds_ref.GetGeoTransform()
	#find new up left index
	yy,xx = np.where(shp_raster_arr==value)
	new_x = old_geom[0] + old_geom[1]*(min(xx)+1)
	new_y = old_geom[3] + old_geom[5]*(min(yy)+1)
	new_geom = (new_x,old_geom[1],old_geom[2],new_y,old_geom[4],old_geom[5])
	#start subsetting
	if len(arr.shape) == 2:
		arr1[shp_raster_arr!=value] = ndata
		new_arr = arr1[min(yy):max(yy)+1,min(xx):max(xx)+1]
	else:
		arr1[:,shp_raster_arr!=value] = ndata
		new_arr = arr1[:,min(yy):max(yy)+1,min(xx):max(xx)+1]
	return new_arr, new_geom

def createTif_Single(fname, res_arr, geom, ds_ref, dtype=gdal.GDT_Int16, ndata=-99):
	driver = gdal.GetDriverByName('GTiff')
	dataset = driver.Create(
			fname,
			res_arr.shape[1],
			res_arr.shape[0],
			1,
			dtype,
			options = ['COMPRESS=LZW'])
	dataset.SetGeoTransform(geom)
	dataset.SetProjection(ds_ref.GetProjection())
	outband = dataset.GetRasterBand(1)
	outband.WriteArray(res_arr)
	outband.FlushCache()
	outband.SetNoDataValue(ndata)
	dataset.FlushCache()

def createTif_Multi(fname, res_arr, geom, ds_ref, dtype=gdal.GDT_Int16, ndata=-99):
	driver = gdal.GetDriverByName('GTiff')
	dataset = driver.Create(
			fname,
			res_arr.shape[2],
			res_arr.shape[1],
			res_arr.shape[0],
			dtype,
			options = ['COMPRESS=LZW'])
	dataset.SetGeoTransform(geom)
	dataset.SetProjection(ds_ref.GetProjection())
	for i in range(res_arr.shape[0]):
		outband = dataset.GetRasterBand(i+1)
		outband.WriteArray(res_arr[i,:,:])
		outband.FlushCache()
		outband.SetNoDataValue(ndata)
	dataset.FlushCache()

#needed files
conus_pf_1k_dem = 'HSProc_Stream5LakesSinks_EP0.1_dem.tif'

avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/Topography/HSProc_Stream5LakesSinks_Ep0.1/'

region_shp = 'Regions.shp'

avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Subdomain_Extraction/Shape_Files/Regions_shp/'

region_shps = [region_shp.split('.')[0]+x for x in ['.shp','.dbf','.prj','.shx','.sbx','.sbn']]

region_raster = 'Regions.tif'

#check if file exits
if not os.path.isfile(conus_pf_1k_dem):	
	print(conus_pf_1k_dem+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	os.system('iget -K '+avra_path_tif+conus_pf_1k_dem+' .')

if not os.path.isfile(region_shp):
	print(region_shp+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	for shp_component_file in region_shps:
		os.system('iget -K '+avra_path_shp+shp_component_file+' .')

#basin id
basin_id = 14

##read dem raster
ds_ref = gdal.Open(conus_pf_1k_dem)
arr_ref = ds_ref.ReadAsArray()


##rasterize region shapefile
if os.path.isfile(region_raster):
	os.remove(region_raster)

rasterize(region_raster,region_shp,ds_ref)

shp_raster_arr = gdal.Open(region_raster).ReadAsArray()

#crop to get a tighter mask
new_dem, new_geom = subset(arr_ref,shp_raster_arr,basin_id,ds_ref)
new_mask, _ = subset(shp_raster_arr,shp_raster_arr,basin_id,ds_ref)
new_mask[new_mask==basin_id] = 1

#create dem and mask text file
np.savetxt('UC_dem.txt',new_dem,fmt='%.2f',delimiter=' ')
np.savetxt('UC_mask.txt',new_mask,fmt='%d',delimiter=' ')






