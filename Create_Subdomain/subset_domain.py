#!/usr/bin/env python3
#required libraries
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gdal
import ogr
from os.path import join, abspath, dirname
import numpy as np
import argparse
import pandas as pd
#from glob import glob
import os
import sys
import pfio
from Define_Watershed import DelinWatershed
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

def latlon2pixzone(xul, yul, dx, dy, lat0, lon0):
    rl = abs((yul - lat0)/dy)
    cu = abs((lon0 - xul)/dx)
    return int(round(rl)), int(round(cu))


def rasterize(out_raster, in_shape, ds_ref,
              dtype=gdal.GDT_Int32, ndata=-99,
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
                                    "ALL_TOUCHED=TRUE"],
                           ) != 0:
        raise Exception("error rasterizing layer: %s" % shp_layer)
    else:
        target_ds.FlushCache()
        return 0


def read_from_file(infile):
    # get extension
    ext = os.path.splitext(os.path.basename(infile))[1]
    if ext in ['.tif', '.tiff']:
        res_arr = gdal.Open(infile).ReadAsArray()

    # parflow ascii file
    elif ext == '.sa':
        with open(infile, 'r') as fi:
            header = fi.readline()
        nx, ny, nz = [int(x) for x in header.strip().split(' ')]
        arr = pd.read_csv(infile, skiprows=1, header=None).values
        res_arr = np.reshape(arr, (nz, ny, nx))[:, ::-1, :]
    # parflow binary file
    elif ext == '.pfb':
        res_arr = pfio.pfread(infile)
    else:
        print('can not read file type '+ext)
        sys.exit()
    return res_arr


def subset(arr, mask_arr, ds_ref, ndata=0):
    arr1 = arr.copy()

    # create new geom
    old_geom = ds_ref.GetGeoTransform()

    # find new up left index
    yy, xx = np.where(mask_arr == 1)
    new_x = old_geom[0] + old_geom[1]*(min(xx)+1)
    new_y = old_geom[3] + old_geom[5]*(min(yy)+1)
    new_geom = (new_x, old_geom[1], old_geom[2],
                new_y, old_geom[4], old_geom[5])

    # start subsetting
    if len(arr.shape) == 2:
        arr1[mask_arr != 1] = ndata
        new_arr = arr1[min(yy):max(yy)+1, min(xx):max(xx)+1]

        # add grid cell to make dimensions as multiple of 32 (nicer PxQxR)
        len_y, len_x = new_arr.shape
        new_len_y = ((len_y//32)+1)*32
        n1 = (new_len_y-len_y)//2
        n2 = new_len_y-len_y-n1
        new_len_x = ((len_x//32)+1)*32
        n3 = (new_len_x-len_x)//2
        n4 = new_len_x-len_x-n3
        return_arr = np.zeros((new_len_y, new_len_x))
        return_arr[n1:-n2, n3:-n4] = new_arr
    else:
        arr1[:, mask_arr != 1] = ndata
        new_arr = arr1[:, min(yy):max(yy)+1, min(xx):max(xx)+1]

        # add grid cell to make dimensions as multiple of 32 (nicer PxQxR)
        _, len_y, len_x = new_arr.shape
        new_len_y = ((len_y//32)+1)*32
        n1 = (new_len_y-len_y)//2
        n2 = new_len_y-len_y-n1
        new_len_x = ((len_x//32)+1)*32
        n3 = (new_len_x-len_x)//2
        n4 = new_len_x-len_x-n3
        return_arr = np.zeros((new_arr.shape[0], new_len_y, new_len_x))
        return_arr[:, n1:-n2, n3:-n4] = new_arr
    return return_arr, new_geom


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a solid file of'
                                                 'masked domain for ParFlow')
    parser.add_argument('-out_name', type=str,
                        help='name of output solidfile (optional)')
    parser.add_argument('-dx', type=int, default=1000,
                        help='spatial resolution of solidfile (optional). '
                             'Default is 1000')
    parser.add_argument('-dz', type=int, default=1000,
                        help='lateral resolution of solidfile (optional). '
                             'Default is 1000')
    parser.add_argument('-printmask', type=int, default=0,
                        help='print mask (optional). Default is 0')
    parser.add_argument('-pfmask', '--conus_pf_1k_mask', type=str,
                        default='conus_1km_PFmask2.tif',
                        help='ParFlow CONUS 1km mask (tif)')
    parser.add_argument('-pflakesmask', '--conus_pf_1k_lakes', type=str,
                        default='conus_1km_PFmask_selectLakesmask.tif',
                        help='ParFlow CONUS 1km lakes mask (tif)')
    parser.add_argument('-pflakesborder', '--conus_pf_1k_lakes_border',
                        type=str,
                        default='conus_1km_PFmask_selectLakesborder.tif',
                        help='ParFlow CONUS 1km lake borders (tif)')
    parser.add_argument('-pfbordertype', '--conus_pf_1k_border_type',
                        type=str, default='1km_PF_BorderCellType.tif',
                        help='ParFlow CONUS 1km border cell type (tif)')
    parser.add_argument('-pfsinks', '--conus_pf_1k_sinks',
                        type=str, default='conus_1km_PFmask_manualsinks.tif',
                        help='ParFlow CONUS 1km sinks (tif)')

    subparsers = parser.add_subparsers(dest='type',
                                       help='subset using three options:')

    # group 1: using shapefile
    parser_a = subparsers.add_parser('shapefile',
                                     help='subset using shapefile and the '
                                          'selected id of watershed')
    parser_a.add_argument('-shp_file', type=str, required=True,
                          help='input shapefile')
    parser_a.add_argument('-id', type=int, nargs='+', required=True,
                          help='id of the selected watershed')
    parser_a.add_argument('-att', type=str, required=False, default='OBJECTID',
                          help='Column name of the shape attribute to use '
                               'during rasterize')

    # group 2: using mask file
    parser_b = subparsers.add_parser('mask', help='subset using a mask file')
    parser_b.add_argument('-mask_file', type=str, required=True,
                          help='input mask file')

    # group 3: using custom watershed
    parser_c = subparsers.add_parser('define_watershed',
                                     help='subset using a newly created'
                                          'watershed')
    parser_c.add_argument('-dir_file', type=str, required=True,
                          help='input direction file',)
    parser_c.add_argument('-outlet_file', type=str, required=True,
                          help='file contains coordinates of outlet points')

    # parsing arguments
    args = parser.parse_args()

    outdir = abspath(dirname(args.out_name))

    # make sure that type is not None. Exit early if it is.
    if args.type is None:
        parser.print_usage()
        print(f'{os.path.basename(__file__)}: error: No method specified. '
              'Please choose either: shapefile, mask, or define_watershed')
        sys.exit(1)

    # Download and install pf-mask-utilities
    if not os.path.isdir('pf-mask-utilities'):
        os.system('git clone https://github.com/smithsg84/pf-mask-utilities.git')
        os.chdir('pf-mask-utilities')
        os.system('make')
        os.chdir('..')

    # required raster files
    conus_pf_1k_mask = args.conus_pf_1k_mask

    # 1 for cells inside domain, 0 for cells outside domain, 2 for sinks
    conus_pf_1k_sinks = args.conus_pf_1k_sinks

    # 1 for lakes, 0 for everything else
    conus_pf_1k_lakes = args.conus_pf_1k_lakes
    conus_pf_1k_lakes_border = args.conus_pf_1k_lakes_border

    # A mask marking with 1 for for cells with an ocean border and
    # 2 for cells with a land border
    conus_pf_1k_border_type = args.conus_pf_1k_border_type

    conus_pf_1k_tifs = [conus_pf_1k_mask, conus_pf_1k_sinks, conus_pf_1k_lakes,
                        conus_pf_1k_lakes_border, conus_pf_1k_border_type]

    # check if file exist, if not we need to login to avra and download.
    # This part requires icommand authorization
    if any([not os.path.isfile(x) for x in conus_pf_1k_tifs]):
        avra_path_tif = '/iplant/home/shared/avra/CONUS2.0/Inputs/domain/'
        print(conus_pf_1k_mask+' does not exist...downloading from avra')
        auth = os.system('iinit')
        if auth != 0:
            print('Authentication failed...exit')
            sys.exit()

        for tif_file in conus_pf_1k_tifs:
            os.system('iget -K '+avra_path_tif+tif_file+' .')

    # read domain raster
    ds_ref = gdal.Open(conus_pf_1k_mask)
    arr_ref = ds_ref.ReadAsArray()

    arr_border_type = gdal.Open(conus_pf_1k_border_type).ReadAsArray()
    lakes_mat = gdal.Open(conus_pf_1k_lakes).ReadAsArray()
    sinks_mat = gdal.Open(conus_pf_1k_sinks).ReadAsArray()

    # deal with optional arguments
    dx = args.dx
    dz = args.dz
    printmask = args.printmask

    if args.type == 'shapefile':
        basin_ids = args.id
        region_shp = args.shp_file

        # check if shapefile exist locally
        avra_path_shp = '/iplant/home/shared/avra/CONUS_1.0/' + \
                        'SteadyState_Final/Subdomain_Extraction/' + \
                        'Shape_Files/Regions_shp/'

        region_shps = [region_shp.split('.')[0]+x for x in
                       ['.shp', '.dbf', '.prj', '.shx', '.sbx', '.sbn']]

        region_raster = f'{outdir}/Regions.tif'

        if not os.path.isfile(region_shp):
            print(region_shp+' does not exist...downloading from avra')
            auth = os.system('iinit')
            if auth != 0:
                print('Authentication failed...exit')
                sys.exit()

            for shp_component_file in region_shps:
                os.system('iget -K '+avra_path_shp+shp_component_file+' .')

        # rasterize region shapefile
        if os.path.isfile(region_raster):
            os.remove(region_raster)

        rasterize(region_raster, region_shp, ds_ref, attribute_name=args.att)

        shp_raster_arr = gdal.Open(region_raster).ReadAsArray()

        # mask array
        mask_arr = np.isin(shp_raster_arr, basin_ids).astype(np.int)

    elif args.type == 'mask':
        mask_file = args.mask_file
        if not os.path.isfile(mask_file):
            print(mask_file+' does not exist...please create one')
            sys.exit()

        file_ext = os.path.splitext(os.path.basename(mask_file))[1]
        if file_ext == '.tif':
            ds_mask = gdal.Open(mask_file)

            # check if mask file has the same projection and extent
            # with the domain mask file
            if any([ds_ref.GetProjection() != ds_mask.GetProjection(),
                    sorted(ds_ref.GetGeoTransform()) !=
                    sorted(ds_mask.GetGeoTransform())]):
                print('mask and domain do not match...exit')
                sys.exit()

        mask_arr = read_from_file(mask_file)

    elif args.type == 'define_watershed':
        dir_file = args.dir_file
        if not os.path.isfile(dir_file):
            print(dir_file+' does not exist...downloading from avra')
            auth = os.system('iinit')
            if auth != 0:
                print('Authentication failed...exit')
                sys.exit()

            avra_path_direction = '/iplant/home/shared/avra/CONUS2.0/' + \
                                  'Inputs/Topography/Str5Ep0/'
            os.system('iget -K '+avra_path_direction+dir_file+' .')

        file_ext = os.path.splitext(os.path.basename(dir_file))[1]
        if file_ext == '.tif':
            ds_dir = gdal.Open(dir_file)

            # check if direction file has the same projection and
            # extent with the domain mask file
            if any([ds_ref.GetProjection() != ds_dir.GetProjection(),
                    sorted(ds_ref.GetGeoTransform()) !=
                    sorted(ds_dir.GetGeoTransform())]):
                print('direction and domain do not match...exit')
                sys.exit()

        outlet_file = args.outlet_file
        if not os.path.isfile(outlet_file):
            print(outlet_file+' does not exist...please create one')
            sys.exit()

        dir_arr = read_from_file(dir_file)
        queue = np.loadtxt(outlet_file)
        queue = queue.reshape(-1, 2)

        # get the mask array from DelinWatershed function
        mask_arr = DelinWatershed(queue, dir_arr, printflag=True)

    if printmask:
        import matplotlib.pyplot as plt
        plt.imshow(arr_ref+mask_arr*2)
        plt.savefig('mask.png')

    # crop to get a tighter mask
    mask_mat, new_geom = subset(arr_ref, mask_arr, ds_ref)
    bordt_mat, _ = subset(arr_border_type, mask_arr, ds_ref)
    lakes_mat_crop, _ = subset(lakes_mat, mask_arr, ds_ref)
    sinks_mat_crop, _ = subset(sinks_mat, mask_arr, ds_ref)

    # create back borders
    # Back borders occur where mask[y+1]-mask[y] is negative
    # (i.e. the cell above is a zero and the cell is inside the mask, i.e. a 1)
    back_mat = np.zeros(mask_mat.shape)
    back_mat[1:, :] = mask_mat[:-1, :] - mask_mat[1:, :]

    # the upper boundary of the top row
    back_mat[0, :] = -1*mask_mat[0, :]
    back_mat[back_mat > 0] = 0

    # changing the valuse to the border types
    back_mat = -1*back_mat*bordt_mat

    # create front borders
    # Front borders occure where mask[y-1]-mask[y] is negative
    # (i.e. the cell above is a zero and the cell is inside the mask, i.e. a 1)
    front_mat = np.zeros(mask_mat.shape)
    front_mat[:-1, :] = mask_mat[1:, :] - mask_mat[:-1, :]

    # the lower boundary of the bottom row
    front_mat[-1, :] = -1*mask_mat[-1, :]
    front_mat[front_mat > 0] = 0

    # changing the valuse to the border types
    front_mat = -1*front_mat*bordt_mat

    # create left borders
    # Left borders occure where mask[x-1]-mask[x] is negative
    left_mat = np.zeros(mask_mat.shape)
    left_mat[:, 1:] = mask_mat[:, :-1] - mask_mat[:, 1:]

    # the left boundary of the first column
    left_mat[:, 0] = -1*mask_mat[:, 0]
    left_mat[left_mat > 0] = 0

    # changing the values to the border types
    left_mat = -1*left_mat*bordt_mat

    # create right borders
    # Right borders occure where mask[x+1]-mask[x] is negative
    right_mat = np.zeros(mask_mat.shape)
    right_mat[:, :-1] = mask_mat[:, 1:] - mask_mat[:, :-1]

    # the right boundary of the last column
    right_mat[:, -1] = -1*mask_mat[:, -1]
    right_mat[right_mat > 0] = 0

    # changing the valuse to the border types
    right_mat = -1*right_mat*bordt_mat

    # deal with top and bottom patches
    # 3 = regular overland boundary
    # 4 = Lake
    # 5 = Sink
    # 6 = bottom

    top_mat = sinks_mat_crop.copy()
    top_mat[top_mat == 1] = 3
    top_mat[lakes_mat_crop > 0] = 4
    top_mat[top_mat == 2] = 5

    bottom_mat = mask_mat * 6

    # create header and write ascii files
    header = 'ncols '+str(mask_mat.shape[1]) + '\n'
    header += 'nrows '+str(mask_mat.shape[0]) + '\n'
    header += 'xllcorner 0.0\n'
    header += 'yllcorner 0.0\n'
    header += 'cellsize '+str(dx) + '\n'
    header += 'NODATA_value 0.0\n'

    patches = {f'{outdir}/Back_Border.asc': back_mat,
               f'{outdir}/Front_Border.asc': front_mat,
               f'{outdir}/Right_Border.asc': right_mat,
               f'{outdir}/Left_Border.asc': left_mat,
               f'{outdir}/Bottom_Border.asc': bottom_mat,
               f'{outdir}/Top_Border.asc': top_mat}

    list_patches = list(patches.keys())

    for patch in patches:
        with open(patch, 'w') as fo:
            fo.write(header)
            np.savetxt(fo, patches[patch].reshape([-1, 1]), '%.3f', ' ')

    # Create solid domain file
    # This part assumes that you have install pf-mask-utilities
    if args.out_name:
        out_vtk = args.out_name+'.vtk'
        out_pfsol = args.out_name+'.pfsol'
    else:
        if args.type == 'shapefile':
            out_name = '_'.join([str(s) for s in basin_ids])
        elif args.type == 'mask':
            out_name = os.path.splitext(os.path.basename(args.mask_file))[0]
        elif args.type == 'define_watershed':
            print('need to specified name of the solid file')
            sys.exit()
        out_vtk = out_name+'.vtk'
        out_pfsol = out_name+'.pfsol'

    if os.path.isfile(out_vtk):
        os.remove(out_vtk)

    if os.path.isfile(out_pfsol):
        os.remove(out_pfsol)

    cmd_create_sol = 'pf-mask-utilities/mask-to-pfsol --mask-back ' + \
                     list_patches[0] + \
                     ' --mask-front '+list_patches[1] + \
                     ' --mask-right '+list_patches[2] + \
                     ' --mask-left '+list_patches[3] + \
                     ' --mask-bottom '+list_patches[4] + \
                     ' --mask-top '+list_patches[5] + \
                     ' --vtk '+out_vtk + \
                     ' --pfsol '+out_pfsol + \
                     ' --depth '+str(dz)

    os.system(cmd_create_sol)

if printmask:
	from pyproj import Proj, transform
		
	ds_geom = ds_ref.GetGeoTransform()
	
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
	shape_1 = min(mask_arr.shape[1],max(xx)+1+n4)-max(0,min(xx)-n3)
	shape_0 = min(mask_arr.shape[0],max(yy)+n2+1)-max(min(yy)-n1,0)
	new_mask = np.zeros((shape_0,shape_1))
	#new_mask[new_mask==1] = top_mat[top_mat!=0]
	#yis, xis = np.nonzero(mask_arr)
	x_centroid, y_centroid = max(0,min(xx)-n3)+shape_1//2, \
						max(min(yy)-n1,0)+shape_0//2
	
	new_mask[n1:n1+top_mat.shape[0],n3:n3+top_mat.shape[1]] = top_mat
	
	lon_centroid = ds_geom[0]+x_centroid*ds_geom[1]
	lat_centroid = ds_geom[3]+y_centroid*ds_geom[5]
	
	lon_0, lat_0 = transform(ds_ref.GetProjection(),Proj(init='epsg:4326'),
							lon_centroid,lat_centroid)
	cmap = plt.get_cmap('tab10', 7)
	
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	ax1 = inset_axes(ax,
                    width="30%", # width = 30% of parent_bbox
                    height=1., # height : 1 inch
                    loc=2)
	im1 = ax.imshow(np.ma.masked_where(top_mat==0,top_mat),cmap=cmap,
					vmin = 0-.5, vmax = 6+.5)
	cb = plt.colorbar(im1,fraction=0.046, pad=0.04, cmap=cmap,
							norm=mpl.colors.Normalize(vmin = -0.5, vmax = 6.5),
							 ticks=np.arange(-0,7,1))
	#cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
	#						norm=mpl.colors.Normalize(vmin = -0.5, vmax = 6.5),
	#						 ticks=np.arange(-0,7,1))
	cb.ax.set_yticklabels(['land','ocean',' ','top','lake','sink','bottom'])
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
				cmap=cmap,vmin = 0-.5, vmax = 6+.5)
	#plt.show()
	
	plt.savefig(out_pfsol.replace('.pfsol','_mask.png'),dpi=300)

