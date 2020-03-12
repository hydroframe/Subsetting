#!/usr/local/bin/python3

#required libraries
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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

def isinteger(x):
	return np.equal(np.mod(x, 1), 0)

def latlon2pixzone(xul, yul, dx, dy, lat0,lon0):
	rl = abs((yul - lat0)/dy)
	cu = abs((lon0 - xul)/dx)
	return int(round(rl)), int(round(cu))

def rasterize(out_raster,in_shape,ds_ref,dtype=gdal.GDT_Int32,ndata=-99):
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

def subset(arr,mask_arr,ds_ref, crop_to_domain, ndata=0):
	arr1 = arr.copy()
	#create new geom
	old_geom = ds_ref.GetGeoTransform()
	#find new up left index
	yy,xx = np.where(mask_arr==1)
	len_y = max(yy)-min(yy)+1
	len_x = max(xx)-min(xx)+1
	###add grid cell to make dimensions as multiple of 32 (nicer PxQxR)
	new_len_y = ((len_y//32)+1)*32
	n1 = (new_len_y-len_y)//2
	n2 = new_len_y-len_y-n1
	new_len_x = ((len_x//32)+1)*32
	n3 = (new_len_x-len_x)//2
	n4 = new_len_x-len_x-n3
	###create new geom
	new_x = old_geom[0] + old_geom[1]*(min(xx)+1)
	new_y = old_geom[3] + old_geom[5]*(min(yy)+1)
	new_geom = (new_x,old_geom[1],old_geom[2],new_y,old_geom[4],old_geom[5])
	new_mask = mask_arr[min(yy)-n1:max(yy)+n2+1,min(xx)-n3:max(xx)+n4+1]
	#start subsetting
	if len(arr.shape) == 2:
		if crop_to_domain:
			arr1[mask_arr!=1] = ndata
			return_arr = np.zeros((new_len_y,new_len_x))
			return_arr[n1:-n2,n3:-n4] = arr1[min(yy):max(yy)+1,min(xx):max(xx)+1]
		else:
			return_arr = arr1[min(yy)-n1:max(yy)+n2+1,min(xx)-n3:max(xx)+n4+1]
		return_arr = return_arr[np.newaxis,...]
	else:
		if crop_to_domain:
			arr1[:,mask_arr!=1] = ndata
			return_arr = np.zeros((arr1.shape[0],new_len_y,new_len_x))
			return_arr[:,n1:-n2,n3:-n4] = arr1[:,min(yy):max(yy)+1,min(xx):max(xx)+1]
		else:
			return_arr = arr1[:,min(yy)-n1:max(yy)+n2+1,min(xx)-n3:max(xx)+n4+1]
	bbox = (min(yy)-n1,max(yy)+n2+1,min(xx)-n3,max(xx)+n4+1)
	return return_arr, new_geom, new_mask, bbox


parser = argparse.ArgumentParser(description='Create a clipped input ParFlow binary files')
parser.add_argument('-i','--input_file',type=str, help='input file for subsetting')
parser.add_argument('--crop_to_domain',type=int, help='crop to domain (i.e. value outside of the domain will be assign as nodata -- optional).Default is 1')
parser.add_argument('-x0',type=int, help='x0 (optional).Default is 0')
parser.add_argument('-y0',type=int, help='y0 (optional).Default is 0')
parser.add_argument('-z0',type=int, help='z0 (optional).Default is 0')
subparsers = parser.add_subparsers(dest='type',help='subset using three options:')

#group 1: using shapefile
parser_a = subparsers.add_parser('shapefile', help='subset using shapefile and the selected id of watershed')
parser_a.add_argument('-shp_file',type=str, help = 'input shapefile')
parser_a.add_argument('-id',type=int, help = 'id of the selected watershed')
parser_a.add_argument('-out_name',type=str, help = 'name of output solidfile (optional)')
parser_a.add_argument('-dx',type=int, help = 'spatial resolution of solidfile (optional). Default is 1000')
parser_a.add_argument('-dz',type=int, help = 'lateral resolution of solidfile (optional). Default is 1000')
parser_a.add_argument('-printmask',type=int, help = 'print mask (optional). Default is 0')
parser_a.add_argument('-printbbox',type=int, help = 'print bounding box (optional). Default is 0')
#parser_a.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
#parser_a.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

#group 2: using mask file
parser_b = subparsers.add_parser('mask', help='subset using a mask file')
parser_b.add_argument('-mask_file',type=str, help = 'input mask file')
parser_b.add_argument('-out_name',type=str, help = 'name of output solidfile (optional)')
parser_b.add_argument('-dx',type=int, help = 'spatial resolution of solidfile (optional). Default is 1000')
parser_b.add_argument('-dz',type=int, help = 'lateral resolution of solidfile (optional). Default is 1000')
parser_b.add_argument('-printmask',type=int, help = 'print mask (optional). Default is 0')
parser_b.add_argument('-printbbox',type=int, help = 'print bounding box (optional). Default is 0')
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
parser_c.add_argument('-printbbox',type=int, help = 'print bounding box (optional). Default is 0')
#parser_c.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
#parser_c.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

###required raster files
conus_pf_1k_mask = '../CONUS1_inputs/conus_1km_PFmask2.tif'

avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/domain/'

###check if file exits, if not we need to login to avra and download. This part requires icommand authorization
if not os.path.isfile(conus_pf_1k_mask):	
	print(conus_pf_1k_mask+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	os.system('iget -K '+avra_path_tif+conus_pf_1k_mask+' ../CONUS1_inputs/')

###read domain raster
ds_ref = gdal.Open(conus_pf_1k_mask)
arr_ref = ds_ref.ReadAsArray()

#parsing arguments
args = parser.parse_args()

###read input file
try:
	infile = args.i
except:
	infile = args.input_file

if not os.path.isfile(infile):
	print(infile+' does not exits...exitting')
	sys.exit()

file_ext = os.path.splitext(os.path.basename(infile))[1]
if file_ext == '.tif':
	ds_in = gdal.Open(infile)

	#check if direction file has the same projection and extent with the domain mask file
	if any([ds_ref.GetProjection() != ds_in.GetProjection(),
			sorted(ds_ref.GetGeoTransform()) != sorted(ds_in.GetGeoTransform())]):
		print ('input and domain do not match...exit')
		sys.exit()

arr_in = read_from_file(infile)

#deal with optional arguments
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

if not args.printmask:
	printmask = 0
else:
	printmask = 1

if not args.printbbox:
	printbbox = 0
else:
	printbbox = 1

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

#main arguments
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
				np.allclose(np.array(ds_ref.GetGeoTransform()),
							np.array(ds_mask.GetGeoTransform()),atol=1e-5)==False]):
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
clip_arr, new_geom, new_mask_x, bbox = subset(arr_in,mask_arr,ds_ref,crop_to_domain)

###create clipped outputs
if args.out_name:
	out_name = args.out_name
	out_name = os.path.basename(out_name).split('.')
	out_name[-1] = 'pfb'
	out_pfb = '.'.join(out_name)
else:
	if args.type == 'shapefile':
		out_name = str(basin_id)
	elif args.type == 'mask':
		out_name = os.path.splitext(os.path.basename(args.mask_file))[0]
	elif args.type == 'define_watershed':
		print ('need to specified name of the solid file')
		sys.exit()
	out_pfb = out_name+'.pfb'

if os.path.isfile(out_pfb):
	os.remove(out_pfb)


pfio.pfwrite(clip_arr,out_pfb,float(x0),float(y0),float(z0),
								float(dx),float(dx),float(dz))

if printbbox:
	out_bbox = out_pfb.replace('.pfb','.txt')
	with open(out_bbox,'w') as fo:
		fo.write('y1\ty2\tx1\tx2\n')
		fo.write('\t'.join('%d' % x for x in bbox))


if printmask:
	from pyproj import Proj, transform
		
	ds_geom = ds_ref.GetGeoTransform()
	
	soil_idx = {1:'s1',2:'s2',3:'s3',4:'s4',5:'s5',6:'s6',7:'s7',8:'s8',
				9:'s9',10:'s10',11:'s11',12:'s12',13:'s13',19:'b1',20:'b2',
				21:'g1',22:'g2',23:'g3',24:'g4',25:'g5',26:'g6',27:'g7',28:'g8'}
	
	yy,xx = np.where(mask_arr==1)
	len_y, len_x = max(yy)+1-min(yy),max(xx)+1-min(xx)
	if len_y < 500 or len_x < 500:
		new_len_y = 800
		new_len_x = 1000
	else:
		new_len_y = len_y*3
		new_len_x = len_x*3
	n1 = (new_len_y-len_y)//2
	n2 = new_len_y-len_y-n1
	n3 = (new_len_x-len_x)//2
	n4 = new_len_x-len_x-n3
	new_mask = mask_arr[max(min(yy)-n1,0):min(mask_arr.shape[0],max(yy)+n2+1),
						max(0,min(xx)-n3):min(mask_arr.shape[1],max(xx)+1+n4)]
	yy1,xx1 = np.where(new_mask==1)
	if 'slope' in out_pfb:
		new_mask[new_mask==1] = clip_arr[0,new_mask_x==1]
		cmap = plt.get_cmap('RdBu')
		norm = mpl.colors.Normalize(vmin = -0.25, vmax = 0.25)
		show_arr = clip_arr[0,:,:]
		pre_name = '_slope'
		ticks = None
		yticklabels = None

	else:
		new_mask[new_mask==1] = clip_arr[-1,new_mask_x==1]
		show_arr = clip_arr[-1,:,:]
		if isinteger(show_arr).all():
			show_arr = np.ma.masked_where(show_arr==0,show_arr)
			vmin, vmax = np.min(show_arr), np.max(show_arr)
			cmap = plt.get_cmap('tab10',vmax-vmin+1)
			norm = mpl.colors.Normalize(vmin = vmin-0.5, vmax = vmax+0.5)
			ticks = np.arange(vmin,vmax+1,1)
			#print(ticks)
			pre_name = '_indi'
			yticklabels=[soil_idx[x] if int(x) in soil_idx else ' ' \
							for x in ticks]
		else:
			show_arr = np.ma.masked_where(show_arr==0,show_arr)
			vmin, vmax = np.min(show_arr), np.max(show_arr)
			cmap = plt.get_cmap('Blues')
			norm = mpl.colors.SymLogNorm(vmin = vmin, vmax = vmax,linthresh=0.03)
			pre_name = '_press'
			ticks = None
			yticklabels = None
		
	#yis, xis = np.nonzero(mask_arr)
	x_centroid, y_centroid = max(0,min(xx)-n3)+(new_mask.shape[1]//2), \
						max(min(yy)-n1,0)+new_mask.shape[0]//2
	
	lon_centroid = ds_geom[0]+x_centroid*ds_geom[1]
	lat_centroid = ds_geom[3]+y_centroid*ds_geom[5]
	
	lon_0, lat_0 = transform(ds_ref.GetProjection(),Proj(init='epsg:4326'),
							lon_centroid,lat_centroid)
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	ax1 = inset_axes(ax,
                    width="30%", # width = 30% of parent_bbox
                    height=1., # height : 1 inch
                    loc=2)
	im1 = ax.imshow(show_arr[:,:],cmap=cmap,norm=norm)
	cb = plt.colorbar(im1,fraction=0.046, pad=0.04, cmap=cmap,
							norm=norm,ticks=ticks)
	#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
	#						norm=mpl.colors.Normalize(vmin = -0.5, vmax = 6.5),
	#						 ticks=np.arange(-0,7,1))
	if ticks is not None:
		cb.ax.set_yticklabels(yticklabels)
	#print(new_mask.shape,lat_0,lon_0)
	m = Basemap(projection='lcc',
			resolution='l',
			rsphere=(6378137.00,6356752.3142),
			width=new_mask.shape[1]*1000,
			height=new_mask.shape[0]*1000,
			lat_1=30,lat_2=60,
			lat_0=lat_0, lon_0=lon_0,ax=ax1
			)
	
	m.drawcoastlines()
	m.drawcountries()
	m.drawstates()
	
	x = np.linspace(0, m.urcrnrx, new_mask.shape[1])
	y = np.linspace(0, m.urcrnry, new_mask.shape[0])
	
	xx, yy = np.meshgrid(x, y)
	
	m.pcolormesh(xx, yy, np.ma.masked_where(new_mask[::-1,:]==0,new_mask[::-1,:]),
				cmap=cmap,norm=norm)
	#plt.show()
	
	plt.savefig(out_pfb.replace('.pfb',pre_name+'.png'),dpi=300)
	

