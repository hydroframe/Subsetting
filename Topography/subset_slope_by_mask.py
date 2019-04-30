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


slope_file = sys.argv[1]

if not slope_file or not os.path.isfile(slope_file) or not  os.path.isfile(slope_file.replace('slopex','slopey')):
	print('the slope files must be exist')
	sys.exit()

mask_file = sys.argv[2]
if not mask_file or not os.path.isfile(mask_file) :
	print('the mask file must be exist')
	sys.exit()

##read mask file
mask_arr = gdal.Open(mask_file).ReadAsArray()

yy,xx = np.where(mask_arr==1)

new_arr = mask_arr[min(yy):max(yy)+1,min(xx):max(xx)+1]

##find the extends
###add n extra grid cells to every direction
n = int(max(new_arr.shape)*0.1)

il, jl, iu, ju = max(min(xx)-n,0), max(min(yy)-n,0), max(xx)+n+1,  max(yy)+n+1

try:
	outnamex = sys.argv[3]
except:
	outnamex = 'sub_'+os.path.splitext(os.path.basename(slope_file))[0]

innamex = os.path.splitext(os.path.basename(slope_file))[0]
innamey = innamex.replace('slopex','slopey')
outnamey = outnamex.replace('slopex','slopey')


##modify the subset tcl file
with open('subset_slopes.tcl','r') as fi:
	tcl_content = fi.read()

tcl_content = tcl_content.split('\n')
new_content = []
for line in tcl_content:
	if '$innamex' in line:
		line = line.replace('$innamex',innamex)
	if '$innamey' in line:
		line = line.replace('$innamey',innamey)
	if '$outnamex' in line:
		line = line.replace('$outnamex',outnamex)
	if '$outnamey' in line:
		line = line.replace('$outnamey',outnamey)
	if '$il' in line:
		line = line.replace('$il',str(il))
	if '$jl' in line:
		line = line.replace('$jl',str(jl))
	if '$iu' in line:
		line = line.replace('$iu',str(iu))
	if '$ju' in line:
		line = line.replace('$ju',str(ju))
	new_content.append(line+'\n')

with open('subset_slopes'+outnamex+'.tcl','w') as fo:
	for line in new_content:
		fo.write(line)

os.system('tclsh subset_slopes'+outnamex+'.tcl')





