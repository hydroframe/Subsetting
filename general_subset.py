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
import shutil
import pfio

#Parsing arguments

parser = argparse.ArgumentParser(description='Create subset input files for ParFlow simulation')
subparsers = parser.add_subparsers(dest='type',help='subset using three options:')

#group 1: using shapefile
parser_a = subparsers.add_parser('shapefile', help='subset using shapefile and the selected id of watershed')
parser_a.add_argument('-shp_file',type=str, help = 'input shapefile')
parser_a.add_argument('-id',type=int, help = 'id of the selected watershed')
parser_a.add_argument('-out_name',type=str, help = 'name of output solidfile (required)')
parser_a.add_argument('-dx',type=int, help = 'spatial resolution of solidfile (optional). Default is 1000')
parser_a.add_argument('-dz',type=int, help = 'lateral resolution of solidfile (optional). Default is 1000')
parser_a.add_argument('-printmask',type=int, help = 'print mask (optional). Default is 0')
#parser_a.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
#parser_a.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

#group 2: using mask file
parser_b = subparsers.add_parser('mask', help='subset using a mask file')
parser_b.add_argument('-mask_file',type=str, help = 'input mask file')
parser_b.add_argument('-out_name',type=str, help = 'name of output solidfile (required)')
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
conus_pf_1k_mask = 'CONUS1_inputs/conus_1km_PFmask2.tif'
conus_pf_1k_sinks = 'CONUS1_inputs/conus_1km_PFmask_manualsinks.tif' #1 for cells inside domain, 0 for cells outside domain, 2 for sinks
conus_pf_1k_lakes = 'CONUS1_inputs/conus_1km_PFmask_selectLakesmask.tif' #1 for lakes, 0 for everything else
conus_pf_1k_lakes_border = 'CONUS1_inputs/conus_1km_PFmask_selectLakesborder.tif'
conus_pf_1k_border_type = 'CONUS1_inputs/1km_PF_BorderCellType.tif' # A mask marking with 1 for for cells with an ocean border and 2 for cells with a land border

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
slopex_tif = 'CONUS1_inputs/Str3Ep0_smth.rvth_1500.mx0.5.mn5.sec0.up_slopex.tif'
avra_path_slope = '/iplant/home/shared/avra/CONUS2.0/Inputs/Topography/Str5Ep0/'
if not os.path.isfile(slopex_tif):
	os.system('iget -K '+avra_path_slope+os.path.basename(slopex_tif)+' CONUS1_inputs/')

slopey_tif = 'CONUS1_inputs/Str3Ep0_smth.rvth_1500.mx0.5.mn5.sec0.up_slopey.tif'
if not os.path.isfile(slopey_tif):
	os.system('iget -K '+avra_path_slope+os.path.basename(slopey_tif)+' CONUS1_inputs/')

###required subsurface file
subsurface_tif = 'CONUS1_inputs/3d-grid.v3.tif'

#check if subsurface_tif is exists
if not os.path.isfile(subsurface_tif):
	print(subsurface_tif+' does not exits...download and process from avra')
	grid_3d_file = '3d-grid.v3.txt'
	avra_path_subsurface = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Input_Development/Subsurface/'+\
		grid_3d_file
	os.system('iget -K '+avra_path_subsurface+' CONUS1_inputs/')
	os.chdir('utils')
	os.system('python3 map_conus_1_to_2.py ../CONUS1_inputs/'+grid_3d_file)
	os.chdir('..')

###required PME file
pme_tif = 'CONUS1_inputs/PME.tif'

if not os.path.isfile(pme_tif):
	print(pme_tif+' does not exits...download and process from avra')
	pme_file = 'PME.txt'
	avra_path_pme = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Input_Development/PME/'+\
		pme_file
	os.system('iget -K '+avra_path_pme+' CONUS1_inputs/')
	os.chdir('utils')
	os.system('python3 map_conus_1_to_2.py ../CONUS1_inputs/'+pme_file)
	os.chdir('..')

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

if not args.out_name:
	print ('need to specified out_name')
	sys.exit()
else:
	out_name = args.out_name

list_conus_inputs = [subsurface_tif,pme_tif,slopex_tif,slopey_tif]

if args.type == 'shapefile':
	basin_id = args.id
	region_shp = args.shp_file
	#create domain
	os.chdir('Create_Subdomain')
	os.system('python3 subset_domain.py shapefile -shp_file '+region_shp+\
					' -id '+str(basin_id)+' -out_name '+out_name+' -printmask '+printmask)
	os.chdir('..')
	#subset input
	os.chdir('Clip_Inputs')
	for input in list_conus_inputs:
		os.system('python3 clip_inputs.py -i ../'+\
					input+' shapefile -shp_file '+region_shp+\
					' -id '+str(basin_id)+' -out_name '+out_name+'_'+\
					os.path.basename(input))
	os.chdir('..')

elif args.type == 'mask':
	mask_file = args.mask_file
	if not os.path.isfile(mask_file):
		print (mask_file+' does not exits...please create one')
		sys.exit()
	#create domain
	os.chdir('Create_Subdomain')
	os.system('python3 subset_domain.py mask -mask_file '+mask_file+\
					' -out_name '+out_name+' -printmask '+printmask)
	os.chdir('..')
	#subset input
	os.chdir('Clip_Inputs')
	for input in list_conus_inputs:
		os.system('python3 clip_inputs.py -i ../'+\
					input+' mask -mask_file '+mask_file+\
					' -out_name '+out_name+'_'+\
					os.path.basename(input))
	os.chdir('..')

elif args.type == 'define_watershed':
	dir_file = args.dir_file
	outlet_file = args.outlet_file
	if not os.path.isfile(outlet_file):
		print (outlet_file+' does not exits...please create one')
		sys.exit()
	
	#create domain
	os.chdir('Create_Subdomain')
	os.system('python3 subset_domain.py define_watershed -dir_file '+dir_file+\
					' -outlet_file '+outlet_file+\
					' -out_name '+out_name+' -printmask '+printmask)
	os.chdir('..')
	#subset input
	os.chdir('Clip_Inputs')
	for input in list_conus_inputs:
		os.system('python3 clip_inputs.py -i ../'+\
					input+' define_watershed -dir_file '+dir_file+\
					' -outlet_file '+outlet_file+\
					' -out_name '+out_name+'_'+\
					os.path.basename(input))
	os.chdir('..')
#move newly created files to input_files folder
if os.path.isdir('input_files/'):
	shutil.rmtree('input_files/')

os.mkdir('input_files/')
os.system('cp Create_Subdomain/'+out_name+'.pfsol input_files/')
os.system('cp Clip_Inputs/'+out_name+'_*.pfb input_files/')

#generate tcl script and run
input_files = sorted(glob('input_files/*'))
os.chdir('Make_Tcl')

os.system('python3 generate_tcl.py -o '+out_name+'.tcl '+\
			'-i parking_lot_template.tcl --runname '+out_name+\
			' -sl ../'+input_files[-1]+\
			' -so ../'+input_files[0]+' -evap 1 '+
			'--evap_file ../'+input_files[2]+' -e 10 --batches 0 3 6'

os.mkdir('run_output')
os.system('cp '+out_name+'.tcl run_output/')
os.chdir('run_output')
os.system('tclsh '+out_name+'.tcl')


