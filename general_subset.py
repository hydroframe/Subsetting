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
from glob import glob
import os
import subprocess
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
parser_a.add_argument('-printbbox',type=int, help = 'print bounding box (optional). Default is 0')
#parser_a.add_argument('-z_bottom',type=int, help = 'bottom of domain (optional). Default is 0')
#parser_a.add_argument('-z_top',type=int, help = 'top of domain (optional). Default is 1000')

###required raster files

if not os.path.isdir('CONUS1_inputs/'):
	os.mkdir('CONUS1_inputs/')

conus_pf_1k_mask = 'CONUS1_inputs/Domain_Blank_Mask.tif'

avra_path_tif = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Other_Domain_Files/'

###check if file exits, if not we need to login to avra and download. This part requires icommand authorization
auth = None
if not os.path.isfile(conus_pf_1k_mask):	
	print(conus_pf_1k_mask+' does not exits...downloading from avra')
	auth = os.system('iinit')
	if auth != 0:
		print('Authentication failed...exit')
		sys.exit()
	
	os.system('iget -vK '+avra_path_tif+os.path.basename(conus_pf_1k_mask)+' CONUS1_inputs/')

###required slope files
slopex_tif = 'CONUS1_inputs/slopex.pfb'
avra_path_slope = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/RunScript_and_Inputs/'
if not os.path.isfile(slopex_tif):
	print(slopex_tif+' does not exits...downloading from avra')
	if auth is None:
		auth = os.system('iinit')
		if auth != 0:
			print('Authentication failed...exit')
			sys.exit()
	os.system('iget -vK '+avra_path_slope+os.path.basename(slopex_tif)+' CONUS1_inputs/')

slopey_tif = 'CONUS1_inputs/slopey.pfb'
if not os.path.isfile(slopey_tif):
	print(slopey_tif+' does not exits...downloading from avra')
	if auth is None:
		auth = os.system('iinit')
		if auth != 0:
			print('Authentication failed...exit')
			sys.exit()
	os.system('iget -vK '+avra_path_slope+os.path.basename(slopey_tif)+' CONUS1_inputs/')

###required subsurface file
subsurface_tif = 'CONUS1_inputs/grid3d.v3.pfb'

#check if subsurface_tif is exists
if not os.path.isfile(subsurface_tif):
	print(subsurface_tif+' does not exits...download and process from avra')
	if auth is None:
		auth = os.system('iinit')
		if auth != 0:
			print('Authentication failed...exit')
			sys.exit()
	os.system('iget -vK '+avra_path_slope+os.path.basename(subsurface_tif)+' CONUS1_inputs/')

	"""
	grid_3d_file = '3d-grid.v3.txt'
	avra_path_subsurface = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Input_Development/Subsurface/'+\
		grid_3d_file
	os.system('iget -K '+avra_path_subsurface+' CONUS1_inputs/')
	os.chdir('utils')
	os.system('python3 map_conus_1_to_2.py ../CONUS1_inputs/'+grid_3d_file)
	os.chdir('..')
	"""

###required PME file
pme_tif = 'CONUS1_inputs/PmE.flux.pfb'

if not os.path.isfile(pme_tif):
	print(subsurface_tif+' does not exits...download and process from avra')
	if auth is None:
		auth = os.system('iinit')
		if auth != 0:
			print('Authentication failed...exit')
			sys.exit()
	print(pme_tif+' does not exits...download and process from avra')
	
	os.system('iget -vK '+avra_path_slope+os.path.basename(pme_tif)+' CONUS1_inputs/')
	"""
	pme_file = 'PME.txt'
	avra_path_pme = '/iplant/home/shared/avra/CONUS_1.0/SteadyState_Final/Input_Development/PME/'+\
		pme_file
	os.system('iget -K '+avra_path_pme+' CONUS1_inputs/')
	os.chdir('utils')
	os.system('python3 map_conus_1_to_2.py ../CONUS1_inputs/'+pme_file)
	os.chdir('..')
	"""

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

if not args.printbbox:
	printbbox = 0
else:
	printbbox = 1

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
	create_sub = subprocess.run(['python3', 'subset_domain.py',
							'shapefile','-shp_file',region_shp,
							'-id',str(basin_id),
							'-out_name',out_name,
							'-printmask',str(printmask),
							'-printbbox', str(printbbox)], stdout=subprocess.PIPE)
	temp_list = create_sub.stdout.decode('utf-8').split('\n')
	batches = ''
	for line in temp_list:
		if 'Number of triangles in patch' in line:
			line = line.strip()
			batches += line.split()[-3]+' '
	#os.system('python3 subset_domain.py shapefile -shp_file '+region_shp+\
	#				' -id '+str(basin_id)+' -out_name '+out_name+' -printmask '+str(printmask))
	os.chdir('..')
	#subset input
	os.chdir('Clip_Inputs')
	for input in list_conus_inputs:
		os.system('python3 clip_inputs.py -i ../'+\
					input+' shapefile -shp_file '+region_shp+\
					' -id '+str(basin_id)+' -out_name '+out_name+'_'+\
					os.path.basename(input)+' -printmask '+str(printmask)+\
					' -printbbox '+str(printbbox))
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
			'--evap_file ../'+input_files[2]+' -e 10 --batches '+batches)

os.chdir('..')

if os.path.isdir('run_output/'):
	shutil.rmtree('run_output/')

os.mkdir('run_output')
os.system('cp Make_Tcl/'+out_name+'.tcl run_output/')
os.chdir('run_output')
os.system('tclsh '+out_name+'.tcl')


