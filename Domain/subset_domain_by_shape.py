#!/usr/local/bin/python3

#required libraries
import gdal, ogr, osr
import numpy as np
#from glob import glob
import os
import sys
#import subprocess
#import matplotlib.pyplot as plt

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
conus_pf_1k_mask = 'conus_1km_PFmask2.tif'
conus_pf_1k_sinks = 'conus_1km_PFmask_manualsinks.tif' #1 for cells inside domain, 0 for cells outside domain, 2 for sinks
conus_pf_1k_lakes = 'conus_1km_PFmask_selectLakesmask.tif' #1 for lakes, 0 for everything else
conus_pf_1k_lakes_border = 'conus_1km_PFmask_selectLakesborder.tif'
conus_pf_1k_border_type = '1km_PF_BorderCellType.tif' # A mask marking with 1 for for cells with an ocean border and 2 for cells with a land border

conus_pf_1k_tifs = [conus_pf_1k_mask,conus_pf_1k_sinks,conus_pf_1k_lakes,
					conus_pf_1k_lakes_border,conus_pf_1k_lakes_border,
					conus_pf_1k_border_type]
avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/domain/'

###required shapefile
region_shp = 'Regions.shp'

avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Subdomain_Extraction/Shape_Files/Regions_shp/'

region_shps = [region_shp.split('.')[0]+x for x in ['.shp','.dbf','.prj','.shx','.sbx','.sbn']]

region_raster = 'Regions.tif'

###check if file exits, if not we need to login to avra and download. This part requires icommand authorization
if any([not os.path.isfile(x) for x in conus_pf_1k_tifs]):	
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



###read domain raster
ds_ref = gdal.Open(conus_pf_1k_mask)
arr_ref = ds_ref.ReadAsArray()

arr_border_type = gdal.Open(conus_pf_1k_border_type).ReadAsArray()
lakes_mat = gdal.Open(conus_pf_1k_lakes).ReadAsArray()
sinks_mat = gdal.Open(conus_pf_1k_sinks).ReadAsArray()

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
mask_mat, new_geom = subset(arr_ref,shp_raster_arr,basin_id,ds_ref)
bordt_mat, _ = subset(arr_border_type,shp_raster_arr,basin_id,ds_ref)
lakes_mat_crop, _ = subset(lakes_mat,shp_raster_arr,basin_id,ds_ref)
sinks_mat_crop, _ = subset(sinks_mat,shp_raster_arr,basin_id,ds_ref)

###create back borders
##Back borders occure where mask[y+1]-mask[y] is negative (i.e. the cell above is a zero and the cell is inside the mask, i.e. a 1)
back_mat=np.zeros(mask_mat.shape)
back_mat[1:,:]=mask_mat[:-1,:] - mask_mat[1:,:]
back_mat[0,:]=-1*mask_mat[0,:] #the upper boundary of the top row
back_mat[back_mat>0]=0
back_mat=-1*back_mat*bordt_mat #changing the valuse to the border types 

###create front borders
##Front borders occure where mask[y-1]-mask[y] is negative (i.e. the cell above is a zero and the cell is inside the mask, i.e. a 1)
front_mat=np.zeros(mask_mat.shape)
front_mat[:-1,:]=mask_mat[1:,:] - mask_mat[:-1,:]
front_mat[-1,:]=-1*mask_mat[-1,:] #the lower boundary of the bottom row
front_mat[front_mat>0]=0
front_mat=-1*front_mat*bordt_mat #changing the valuse to the border types 

###create left borders
##Left borders occure where mask[x-1]-mask[x] is negative 
left_mat=np.zeros(mask_mat.shape)
left_mat[:,1:]=mask_mat[:,:-1] - mask_mat[:,1:]
left_mat[:,0]=-1*mask_mat[:,0] #the left boundary of the first column
left_mat[left_mat>0]=0
left_mat=-1*left_mat*bordt_mat #changing the valuse to the border types 

###create right borders
##Right borders occure where mask[x+1]-mask[x] is negative 
right_mat=np.zeros(mask_mat.shape)
right_mat[:,:-1]=mask_mat[:,1:] - mask_mat[:,:-1]
right_mat[:,-1]=-1*mask_mat[:,-1] #the right boundary of the last column
right_mat[right_mat>0]=0
right_mat=-1*right_mat*bordt_mat #changing the valuse to the border types 

###deal with top and bottom patches
# 3 = regular overland boundary
# 4 = Lake
# 5 = Sink
# 6 = bottom

top_mat=sinks_mat_crop.copy()
top_mat[top_mat==1]=3
top_mat[lakes_mat_crop>0]=4
top_mat[top_mat==2]=5

bottom_mat=mask_mat*6

###create header and write ascii files
header = 'ncols '+str(mask_mat.shape[1])+'\n'
header += 'nrows '+str(mask_mat.shape[0])+'\n'
header += 'xllcorner 0.0\n'
header += 'yllcorner 0.0\n'
header += 'cellsize 1000.0\n'
header += 'NODATA_value 0.0\n'

patches = {'Back_Border.asc':back_mat,
			'Front_Border.asc':front_mat,
			'Right_Border.asc':right_mat,
			'Left_Border.asc':left_mat,
			'Bottom_Border.asc':bottom_mat,
			'Top_Border.asc':top_mat}

list_patches = list(patches.keys())

for patch in patches:
	with open(patch,'w') as fo:
		fo.write(header)
		np.savetxt(fo,patches[patch].reshape([-1,1]),'%.3f',' ')

###Create solid domain file
#This part assumes that you have install pf-mask-utilities
try:
	if select_type == '-id':
		out_name = sys.argv[3]
	elif select_type == '-p':
		out_name = sys.argv[4]

	out_vtk = out_name+'.vtk'
	out_pfsol = out_name+'.pfsol'
except:
	out_vtk = str(basin_id)+'.vtk'
	out_pfsol = str(basin_id)+'.pfsol'

if os.path.isfile(out_vtk):
	os.remove(out_vtk)

if os.path.isfile(out_pfsol):
	os.remove(out_pfsol)

cmd_create_sol = 'pf-mask-utilities/mask-to-pfsol --mask-back '+list_patches[0]+\
					' --mask-front '+list_patches[1]+\
					' --mask-right '+list_patches[2]+\
					' --mask-left '+list_patches[3]+\
					' --mask-bottom '+list_patches[4]+\
					' --mask-top '+list_patches[5]+\
					' --vtk '+out_vtk+\
					' --pfsol '+out_pfsol

os.system(cmd_create_sol)



