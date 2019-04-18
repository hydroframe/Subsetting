#!/usr/local/bin/python3

#required libraries
import gdal, ogr, osr
#import shapely
import numpy as np
import os
import sys

def latlon2pixzone(xul, yul, dx, dy, lat0,lon0):
	rl = abs((yul - lat0)/dy)
	cu = abs((lon0 - xul)/dx)
	return int(round(rl)), int(round(cu))

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
	return return_arr, new_geom

###select feature method: either select feature by id or by point coordinate
select_type = sys.argv[1]

#print (select_type)
if (select_type != '-id' and select_type != '-p') :
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
conus_pf_1k_dem = 'HSProc_Stream5LakesSinks_EP0.1_dem.tif'
river_mask_sa = 'HSProc_Stream5LakesSinks_EP0.1_RivMask.sa'

avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/Topography/HSProc_Stream5LakesSinks_Ep0.1/'

###required shapefile
region_shp = 'Regions.shp'

avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Subdomain_Extraction/Shape_Files/Regions_shp/'

region_shps = [region_shp.split('.')[0]+x for x in ['.shp','.dbf','.prj','.shx','.sbx','.sbn']]

region_raster = 'Regions.tif'

###check if file exits, if not we need to login to avra and download. This part requires icommand authorization
if not os.path.isfile(conus_pf_1k_dem):	
	print(conus_pf_1k_dem+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	os.system('iget -K '+avra_path_tif+conus_pf_1k_dem+' .')
	os.system('iget -K '+avra_path_tif+river_mask_sa+' .')

if not os.path.isfile(region_shp):
	print(region_shp+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	for shp_component_file in region_shps:
		os.system('iget -K '+avra_path_shp+shp_component_file+' .')


###read dem raster
ds_ref = gdal.Open(conus_pf_1k_dem)
arr_ref = ds_ref.ReadAsArray()

##read river mask
river_raw = np.loadtxt(river_mask_sa,skiprows=1)
river_arr = river_raw.reshape(arr_ref.shape)[::-1,:]

###rasterize region shapefile
if os.path.isfile(region_raster):
	os.remove(region_raster)

rasterize(region_raster,region_shp,ds_ref)

shp_raster_arr = gdal.Open(region_raster).ReadAsArray()
if select_type == '-p':
	##get current file projection
	outSrs = osr.SpatialReference(wkt=ds_ref.GetProjection())
	ds_geom = ds_ref.GetGeoTransform()
	
	##Transform between EPSG:4326 and current file projection:
	inSrs = osr.SpatialReference()
	inSrs.ImportFromEPSG(4326)
	coordTransform = osr.CoordinateTransformation(inSrs, outSrs)
	#create a geometry from coordinates
	inPoint = ogr.Geometry(ogr.wkbPoint)
	inPoint.AddPoint(p_lon,p_lat)
	#transform point
	inPoint.Transform(coordTransform)
	new_p_lat, new_p_lon = inPoint.GetY(), inPoint.GetX()
	
	##get basin id of which contains point:
	p_y, p_x = latlon2pixzone(ds_geom[0], ds_geom[3], 
								ds_geom[1], ds_geom[5], 
								new_p_lat,new_p_lon)
	basin_id = shp_raster_arr[p_y,p_x]

###crop to get a tighter mask
new_dem, new_geom = subset(arr_ref,shp_raster_arr,basin_id,ds_ref)
new_river, _ = subset(river_arr,shp_raster_arr,basin_id,ds_ref)
new_mask, _ = subset(shp_raster_arr,shp_raster_arr,basin_id,ds_ref)
new_mask[new_mask==basin_id] = 1

try:
	if select_type == '-id':
		out_name = sys.argv[3]
	elif select_type == '-p':
		out_name = sys.argv[4]

	out_dem = out_name+'_dem.txt'
	out_river = out_name+'_river.txt'
	out_mask = out_name+'_mask.txt'
except:
	out_dem = str(basin_id)+'_dem.txt'
	out_river = str(basin_id)+'_river.txt'
	out_mask = str(basin_id)+'_mask.txt'

###create dem, river and mask text file
np.savetxt(out_dem,new_dem,fmt='%.2f',delimiter=' ')
np.savetxt(out_river,new_river,fmt='%d',delimiter=' ')
np.savetxt(out_mask,new_mask,fmt='%d',delimiter=' ')






