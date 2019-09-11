#!/usr/local/bin/python3

#required libraries
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gdal, ogr, osr
import numpy as np
import argparse
import pandas as pd
#from glob import glob
import os
import sys
import pfio

#Parsing arguments

parser = argparse.ArgumentParser(description='Create subset input files for ParFlow simulation')
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

###required slope files


###required subsurface file
subsurface_tif = '3d-grid.v3.tif'

#check if subsurface_tif is exists
if not os.path.isfile(subsurface_tif):
	print(subsurface_tif+' does not exits...download and process from avra')
	grid_3d_file = '3d-grid.v3.txt'
	avra_path_subsurface = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Input_Development/Subsurface/'+\
		grid_3d_file
	os.system('iget -K '+avra_path_subsurface+' .')
	os.chdir('utils')
	os.system('python3 map_conus_1_to_2.py ..'+grid_3d_file)
	os.chdir('..')
	

###required PME file
pme_tif = 'PME.tif'

if not os.path.isfile(pme_tif):
	print(pme_tif+' does not exits...download and process from avra')
	pme_file = 'PME.txt'
	avra_path_pme = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Input_Development/PME/'+\
		pme_file
	os.system('iget -K '+avra_path_pme+' .')
	os.chdir('utils')
	os.system('python3 map_conus_1_to_2.py ..'+pme_file)
	os.chdir('..')


















