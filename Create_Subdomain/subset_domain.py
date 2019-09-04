#!/usr/local/bin/python3

#required libraries
import gdal, ogr, osr
import numpy as np
import pandas as pd
import argparse
#from glob import glob
import os
import sys
import pfio
from Define_Watershed import DelinWatershed
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

def read_from_file(infile):
	#get extension
	ext = os.path.splitext(os.path.basename(infile))[1]
	if ext in ['.tif','.tiff']:
		res_arr = gdal.Open(infile).ReadAsArray()
	elif ext == '.sa': #parflow ascii file
		with open(infile, 'r') as fi:
			header = fi.readline()
		nx, ny, nz = [int(x) for x in header.strip().split(' ')]
		arr = pd.read_csv(infile, skiprows=1, header=None).values
		res_arr = np.reshape(arr,(nz,ny,nx))[:,::-1,:]
	elif ext == '.pfb': #parflow binary file
		res_arr = pfio.pfread(infile)
	else:
		print('can not read file type '+ext)
		sys.exit()
	return res_arr

def subset(arr,mask_arr,ds_ref, ndata=0):
	arr1 = arr.copy()
	#create new geom
	old_geom = ds_ref.GetGeoTransform()
	#find new up left index
	yy,xx = np.where(mask_arr==1)
	new_x = old_geom[0] + old_geom[1]*(min(xx)+1)
	new_y = old_geom[3] + old_geom[5]*(min(yy)+1)
	new_geom = (new_x,old_geom[1],old_geom[2],new_y,old_geom[4],old_geom[5])
	#start subsetting
	if len(arr.shape) == 2:
		arr1[mask_arr!=1] = ndata
		new_arr = arr1[min(yy):max(yy)+1,min(xx):max(xx)+1]
		###add grid cell to make dimensions as multiple of 32 (nicer PxQxR)
		len_y, len_x = new_arr.shape
		new_len_y = ((len_y//32)+1)*32
		n1 = (new_len_y-len_y)//2
		n2 = new_len_y-len_y-n1
		new_len_x = ((len_x//32)+1)*32
		n3 = (new_len_x-len_x)//2
		n4 = new_len_x-len_x-n3
		return_arr = np.zeros((new_len_y,new_len_x))
		return_arr[n1:-n2,n3:-n4] = new_arr
	else:
		arr1[:,mask_arr!=1] = ndata
		new_arr = arr1[:,min(yy):max(yy)+1,min(xx):max(xx)+1]
		###add grid cell to make dimensions as multiple of 32 (nicer PxQxR)
		_, len_y, len_x = new_arr.shape
		new_len_y = ((len_y//32)+1)*32
		n1 = (new_len_y-len_y)//2
		n2 = new_len_y-len_y-n1
		new_len_x = ((len_x//32)+1)*32
		n3 = (new_len_x-len_x)//2
		n4 = new_len_x-len_x-n3
		return_arr = np.zeros((new_arr.shape[0],new_len_y,new_len_x))
		return_arr[:,n1:-n2,n3:-n4] = new_arr
	return return_arr, new_geom


parser = argparse.ArgumentParser(description='Create a solid file of masked domain for ParFlow')
subparsers = parser.add_subparsers(dest='type',help='subset using three options:')

#group 1: using shapefile
parser_a = subparsers.add_parser('shapefile', help='subset using shapefile and the selected id of watershed')
parser_a.add_argument('-shp_file',type=str, help = 'input shapefile')
parser_a.add_argument('-id',type=int, help = 'id of the selected watershed')
parser_a.add_argument('-out_name',type=str, help = 'name of output solidfile (optional)')
parser_a.add_argument('-dx',type=int, help = 'spatial resolution of solidfile (optional). Default is 1000')
parser_a.add_argument('-dz',type=int, help = 'lateral resolution of solidfile (optional). Default is 1000')
parser_a.add_argument('-printmask',type=int, help = 'print mask (optional). Default is 0')
#parser_a.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
#parser_a.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

#group 2: using mask file
parser_b = subparsers.add_parser('mask', help='subset using a mask file')
parser_b.add_argument('-mask_file',type=str, help = 'input mask file')
parser_b.add_argument('-out_name',type=str, help = 'name of output solidfile (optional)')
parser_b.add_argument('-dx',type=int, help = 'spatial resolution of solidfile (optional). Default is 1000')
parser_b.add_argument('-dz',type=int, help = 'lateral resolution of solidfile (optional). Default is 1000')
parser_b.add_argument('-printmask',type=int, help = 'print mask (optional). Default is 0')
#parser_b.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
#parser_b.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

#group 3: using custom watershed
parser_c = subparsers.add_parser('define_watershed', help='subset using a newly created watershed')
parser_c.add_argument('-dir_file',type=str, help = 'input direction file',)
parser_c.add_argument('-outlet_file',type=str, help = 'file contains coordinates of outlet points')
parser_c.add_argument('-out_name',type=str, help = 'name of output solidfile (required)')
parser_c.add_argument('-dx',type=int, help = 'spatial resolution of solidfile (optional). Default is 1000')
parser_c.add_argument('-dz',type=int, help = 'lateral resolution of solidfile (optional). Default is 1000')
parser_c.add_argument('-printmask',type=int, help = 'print mask (optional). Default is 0')
#parser_c.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
#parser_c.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

###Download and install pf-mask-utilities
if not os.path.isdir('pf-mask-utilities'):
	os.system('git clone https://github.com/smithsg84/pf-mask-utilities.git')
	os.chdir('pf-mask-utilities')
	os.system('make')
	os.chdir('..')

###required raster files
conus_pf_1k_mask = 'conus_1km_PFmask2.tif'
conus_pf_1k_sinks = 'conus_1km_PFmask_manualsinks.tif' #1 for cells inside domain, 0 for cells outside domain, 2 for sinks
conus_pf_1k_lakes = 'conus_1km_PFmask_selectLakesmask.tif' #1 for lakes, 0 for everything else
conus_pf_1k_lakes_border = 'conus_1km_PFmask_selectLakesborder.tif'
conus_pf_1k_border_type = '1km_PF_BorderCellType.tif' # A mask marking with 1 for for cells with an ocean border and 2 for cells with a land border

conus_pf_1k_tifs = [conus_pf_1k_mask,conus_pf_1k_sinks,conus_pf_1k_lakes,
					conus_pf_1k_lakes_border,conus_pf_1k_border_type]
avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/domain/'

###check if file exits, if not we need to login to avra and download. This part requires icommand authorization
if any([not os.path.isfile(x) for x in conus_pf_1k_tifs]):	
	print(conus_pf_1k_mask+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	for tif_file in conus_pf_1k_tifs:
		os.system('iget -K '+avra_path_tif+tif_file+' .')

###read domain raster
ds_ref = gdal.Open(conus_pf_1k_mask)
arr_ref = ds_ref.ReadAsArray()

arr_border_type = gdal.Open(conus_pf_1k_border_type).ReadAsArray()
lakes_mat = gdal.Open(conus_pf_1k_lakes).ReadAsArray()
sinks_mat = gdal.Open(conus_pf_1k_sinks).ReadAsArray()


#parsing arguments
args = parser.parse_args()

#deal with optional arguments
if not args.dx:
	dx = 1000
else:
	dx = args.dx

if not args.dz:
	dz = 1000
else:
	dz = args.dz

if not args.printmask:
	printmask = 0
else:
	printmask = 1

if args.type == 'shapefile':
	basin_id = args.id
	region_shp = args.shp_file
	#check if shapefile exits locally
	avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Subdomain_Extraction/Shape_Files/Regions_shp/'
	
	region_shps = [region_shp.split('.')[0]+x for x in ['.shp','.dbf','.prj','.shx','.sbx','.sbn']]
	
	region_raster = 'Regions.tif'
	
	if not os.path.isfile(region_shp):
		print(region_shp+' does not exits...downloading from avra')
		auth = os.system('iinit')
		if auth != 0:
			print('Authentication failed...exit')
			sys.exit()
	
		for shp_component_file in region_shps:
			os.system('iget -K '+avra_path_shp+shp_component_file+' .')
	
	
	###rasterize region shapefile
	if os.path.isfile(region_raster):
		os.remove(region_raster)

	rasterize(region_raster,region_shp,ds_ref)
	
	shp_raster_arr = gdal.Open(region_raster).ReadAsArray()
	
	###mask array
	mask_arr = (shp_raster_arr == basin_id).astype(np.int)

elif args.type == 'mask':
	mask_file = args.mask_file
	if not os.path.isfile(mask_file):
		print (mask_file+' does not exits...please create one')
		sys.exit()
	
	file_ext = os.path.splitext(os.path.basename(mask_file))[1]
	if file_ext == '.tif':
		ds_mask = gdal.Open(mask_file)
	
		#check if mask file has the same projection and extent with the domain mask file
		if any([ds_ref.GetProjection() != ds_mask.GetProjection(),
				sorted(ds_ref.GetGeoTransform()) != sorted(ds_mask.GetGeoTransform())]):
			print ('mask and domain do not match...exit')
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
		
		avra_path_direction = '/iplant/home/shared/avra/CONUS2.0/Inputs/Topography/Str5Ep0/'
		os.system('iget -K '+avra_path_direction+dir_file+' .')
	
	file_ext = os.path.splitext(os.path.basename(dir_file))[1]
	if file_ext == '.tif':
		ds_dir = gdal.Open(dir_file)
	
		#check if direction file has the same projection and extent with the domain mask file
		if any([ds_ref.GetProjection() != ds_dir.GetProjection(),
				sorted(ds_ref.GetGeoTransform()) != sorted(ds_dir.GetGeoTransform())]):
			print ('direction and domain do not match...exit')
			sys.exit()
	
	outlet_file = args.outlet_file
	if not os.path.isfile(outlet_file):
		print (outlet_file+' does not exits...please create one')
		sys.exit()
	
	dir_arr = read_from_file(dir_file)
	queue = np.loadtxt(outlet_file)
	queue = queue.reshape(-1,2)
	
	#get the mask array from DelinWatershed function
	mask_arr = DelinWatershed(queue, dir_arr,printflag=True)

###crop to get a tighter mask
mask_mat, new_geom = subset(arr_ref,mask_arr,ds_ref)
bordt_mat, _ = subset(arr_border_type,mask_arr,ds_ref)
lakes_mat_crop, _ = subset(lakes_mat,mask_arr,ds_ref)
sinks_mat_crop, _ = subset(sinks_mat,mask_arr,ds_ref)

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
header += 'cellsize '+str(dx)+'\n'
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
if args.out_name:
	out_vtk = args.out_name+'.vtk'
	out_pfsol = args.out_name+'.pfsol'
else:
	if args.type == 'shapefile':
		out_name = str(basin_id)
	elif args.type == 'mask':
		out_name = os.path.splitext(os.path.basename(args.mask_file))[0]
	elif args.type == 'define_watershed':
		print ('need to specified name of the solid file')
		sys.exit()
	out_vtk = out_name+'.vtk'
	out_pfsol = out_name+'.pfsol'

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
					' --pfsol '+out_pfsol+\
					' --depth '+str(dz)

os.system(cmd_create_sol)

if printmask:
	import matplotlib.pyplot as plt
	from mpl_toolkits.basemap import Basemap
	from pyproj import Proj, transform
	
	ds_geom = ds_ref.GetGeoTransform()
	
	yy,xx = np.where(mask_arr==1)
	len_y, len_x = max(yy)+1-min(yy),max(xx)+1-min(xx)
	new_len_y = len_y*5
	new_len_x = len_x*5
	n1 = (new_len_y-len_y)//2
	n2 = new_len_y-len_y-n1
	n3 = (new_len_x-len_x)//2
	n4 = new_len_x-len_x-n3
	new_mask = mask_arr[max(min(yy)-n1,0):min(mask_arr.shape[0],max(yy)+n2+1),
						max(0,min(xx)-n3):min(mask_arr.shape[1],max(xx)+1+n4)]
	
	yis, xis = np.nonzero(mask_arr)
	x_centroid, y_centroid = xis.mean(), yis.mean()
	
	lon_centroid = ds_geom[0]+x_centroid*ds_geom[1]
	lat_centroid = ds_geom[3]+y_centroid*ds_geom[5]
	
	lon_0, lat_0 = transform(ds_ref.GetProjection(),Proj(init='epsg:4326'),
							lon_centroid,lat_centroid)
	print(new_mask.shape,lat_0,lon_0)
	m = Basemap(projection='lcc',
			resolution='l',
			rsphere=(6378137.00,6356752.3142),
			width=new_mask.shape[1]*10000,
			height=new_mask.shape[0]*10000,
			lat_1=30,lat_2=60,
			lat_0=lat_0, lon_0=lon_0
			)
	
	m.drawcoastlines()
	m.drawcountries()
	m.drawstates()
	
	x = np.linspace(0, m.urcrnrx, new_mask.shape[1])
	y = np.linspace(0, m.urcrnry, new_mask.shape[0])
	
	xx, yy = np.meshgrid(x, y)
	
	m.pcolormesh(xx, yy, np.ma.masked_where(new_mask<-100,new_mask))
	plt.savefig('mask.png')