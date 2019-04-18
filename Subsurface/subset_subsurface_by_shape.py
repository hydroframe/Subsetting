#!/usr/local/bin/python3

#required libraries
import gdal, ogr, osr
import numpy as np
import os
import sys
#import matplotlib.pyplot as plt

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
	#start subsetting
	if len(arr.shape) == 2:
		arr1[shp_raster_arr!=value] = ndata
		new_arr = arr1[min(yy):max(yy)+1,min(xx):max(xx)+1]
		###add n extra grid cells to every direction
		n = int(max(new_arr.shape)*0.02) #proportional to new array dimensions
		#print (n)
		return_arr = np.zeros((new_arr.shape[0]+2*n,new_arr.shape[1]+2*n))
		return_arr[n:-n,n:-n] = new_arr
	else:
		arr1[:,shp_raster_arr!=value] = ndata
		new_arr = arr1[:,min(yy):max(yy)+1,min(xx):max(xx)+1]
		###add n extra grid cells to every direction
		n = int(max(new_arr.shape)*0.02) #proportional to new array dimensions
		#print (n)
		return_arr = np.zeros((new_arr.shape[0],new_arr.shape[1]+2*n,new_arr.shape[2]+2*n))
		return_arr[:,n:-n,n:-n] = new_arr
	return return_arr

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

###select feature method: either select feature by id or by point coordinate
select_type = sys.argv[1]
#print (select_type)
if select_type != '-id' and select_type != '-p':
	print('select type should be either: (1) basin_id: -id or (2)point coordinate: -p')
	sys.exit()
elif select_type == '-id':
	basin_id = sys.argv[2]
	#print(basin_id)
	basin_id = int(basin_id)
elif select_type == '-p':
	p_lat = sys.argv[2]
	p_lon = sys.argv[3]
	p_lat = float(p_lat)
	p_lon = float(p_lon)
else:
	print('select type should be either: (1) basin_id: -id or (2)point coordinate: -p')
	sys.exit()

###required raster files
conus1_domain = 'CONUS1_domain.tif'
avra_path_domain1 = '/iplant/home/hoangtv1899/'
conus2_domain = 'conus_1km_PFmask2.tif'
avra_path_domain2 = '/iplant/home/shared/avra/CONUS2.0/Inputs/domain/'

###required subsurface file
sub_file = '3d-grid.v3.txt'
avra_path_sub = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Input_Development/Subsurface/'

###required shapefile
region_shp = 'Regions.shp'

avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Subdomain_Extraction/Shape_Files/Regions_shp/'

region_shps = [region_shp.split('.')[0]+x for x in ['.shp','.dbf','.prj','.shx','.sbx','.sbn']]

region_raster = 'Regions.tif'

###check if file exits, if not we need to login to avra and download. This part requires icommand authorization
if not os.path.isfile(conus1_domain):	
	print(conus1_domain+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	os.system('iget -K '+avra_path_domain1+conus1_domain+' .')

if not os.path.isfile(conus2_domain):	
	print(conus2_domain+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	os.system('iget -K '+avra_path_domain2+conus2_domain+' .')

if not os.path.isfile(sub_file):
	print(sub_file+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	os.system('iget -K '+avra_path_sub+sub_file+' .')

if not os.path.isfile(region_shp):
	print(region_shp+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	for shp_component_file in region_shps:
		os.system('iget -K '+avra_path_shp+shp_component_file+' .')

###read domain raster
ds_ref1 = gdal.Open(conus1_domain)
arr_ref1 = ds_ref1.ReadAsArray()
ds_ref2 = gdal.Open(conus2_domain)
arr_ref2 = ds_ref2.ReadAsArray()
geom2 = ds_ref2.GetGeoTransform()

##rasterize region shapefile
if os.path.isfile(region_raster):
	os.remove(region_raster)

rasterize(region_raster,region_shp,ds_ref2)

shp_raster_arr = gdal.Open(region_raster).ReadAsArray()

###open subsurface data
subsurface_arr = reshape3d(sub_file, ds_ref1)
##create tif file with old projection
sub_file1 = '3d-grid-v1.tif'
if os.path.isfile(sub_file1):
	os.remove(sub_file1)

createTif_Multi(sub_file1, subsurface_arr, 
				ds_ref1.GetGeoTransform(), ds_ref1, 
				dtype=gdal.GDT_Int16, ndata=-99)

###reproject to conus2.0 domain
sub_file2 = '3d-grid-v2.tif'
if os.path.isfile(sub_file2):
	os.remove(sub_file2)

gdal.Warp(sub_file2,sub_file1,
			srcSRS=ds_ref1.GetProjection(),
			dstSRS=ds_ref2.GetProjection(),
			xRes=geom2[1],
			yRes=geom2[1],
			outputBounds=(geom2[0],geom2[3]+geom2[5]*arr_ref2.shape[0],
							geom2[0]+geom2[1]*arr_ref2.shape[1],geom2[3]))

###Open the newly reprojected file
new_subsurface_arr = gdal.Open(sub_file2).ReadAsArray()

###subset
clip_subsurface_arr = subset(new_subsurface_arr,shp_raster_arr,basin_id,ds_ref2)

try:
	if select_type == '-id':
		out_name = sys.argv[3]
	elif select_type == '-p':
		out_name = sys.argv[4]

except:
	out_name = str(basin_id)

###write to .sa file
clip_subsurface_sa = out_name+'.sa'

if os.path.isfile(clip_subsurface_sa):
	os.remove(clip_subsurface_sa)

with open(clip_subsurface_sa,'w') as fo:
	fo.write(str(clip_subsurface_arr.shape[2])+' '+str(clip_subsurface_arr.shape[1])+\
			' '+str(clip_subsurface_arr.shape[0])+'\n')
	np.savetxt(fo,clip_subsurface_arr[:,::-1,:].reshape(-1,1))

###convert to .pfb and .silo files
with open('File_Conversion.tcl','r') as fi:
	tcl_content = fi.read()

tcl_content = tcl_content.split('\n')
new_content = []
for line in tcl_content:
	if '$out_name' in line:
		line = line.replace('$out_name',out_name)
	if '$nrow' in line:
		line = line.replace('$nrow',str(clip_subsurface_arr.shape[1]))
	if '$ncol' in line:
		line = line.replace('$ncol',str(clip_subsurface_arr.shape[2]))
	if '$dim' in line:
		line = line.replace('$dim',str(clip_subsurface_arr.shape[0]))
	new_content.append(line+'\n')

with open('File_Conversion_'+out_name+'.tcl','w') as fo:
	for line in new_content:
		fo.write(line)

os.system('tclsh File_Conversion_'+out_name+'.tcl')

print('done!')
