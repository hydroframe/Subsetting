#!/usr/local/bin/python3

#required libraries
import gdal, ogr, osr
import numpy as np
#from glob import glob
import os
import sys
#import subprocess
#import matplotlib.pyplot as plt

###mask name
mask_name = sys.argv[1]
###slope name
slope_name = sys.argv[2]
###template file
temp_file = sys.argv[3]
###run name
try:
	run_name = sys.argv[4]
except:
	run_name = mask_name+'_'+slope_name

###Check if pfsol, pfb file exits
domain_path = '../Domain/'
topo_path = '../Topography/'

pfsol_file = domain_path+mask_name+'.pfsol'
asc_file = domain_path+'Top_Border.asc'
pfb_slopex = topo_path+mask_name+'_'+slope_name+'.pfb'
pfb_slopey = topo_path+mask_name+'_'+slope_name.replace('slopex','slopey')+'.pfb'

for req_file in [pfsol_file,pfb_slopex,pfb_slopey,asc_file]:
	if not os.path.isfile(req_file):
		print('missing '+req_file+'...please check in Domain or Topography generation')
		sys.exit()

###Check if template file exits
if not os.path.isfile(temp_file):
	print ('template file is missing')
	sys.exit()

###If all the file exits, copy them to current folder
for req_file in [pfsol_file,pfb_slopex,pfb_slopey,asc_file]:
	os.system('cp '+req_file+' .')


###Get dimensions of the files
with open(asc_file,'r') as fi:
	head = [next(fi) for x in range(6)]

for line in head:
	line = line.strip()
	if 'ncols' in line:
		nx = int(line.split()[1])
	if 'nrows' in line:
		ny = int(line.split()[1])


###read current template file
with open(temp_file,'r') as fi:
	tcl_content = fi.read()

tcl_content = tcl_content.split('\n')

new_content = []
for line in tcl_content:
	if 'set runname' in line: #change runname
		line = line.replace('CONUS2_Parkinglot',run_name)
	if 'pfset ComputationalGrid.NX' in line: #set computational grid nx
		line = line.replace('4442',str(nx))
	if 'pfset ComputationalGrid.NY' in line: #set computational grid ny
		line = line.replace('3256',str(ny))
	if 'pfset GeomInput.domaininput.FileName' in line: #set name of the solid file
		line = line.replace('CONUS2.0_test1.pfsol',os.path.basename(pfsol_file))
	if 'pfset TopoSlopesY.Type' in line: #set slope .pfb file as input
		line = 'pfset TopoSlopesY.Type PFBFile'
	if 'pfset TopoSlopesX.Type' in line: #set slope .pfb file as input
		line = 'pfset TopoSlopesX.Type PFBFile'
	if 'pfset TopoSlopesX.FileName' in line: #set name of the pfb file
		line = 'pfset TopoSlopesX.FileName '+os.path.basename(pfb_slopex)
	if 'pfset TopoSlopesY.FileName' in line: #set name of the pfb file
		line = 'pfset TopoSlopesY.FileName '+os.path.basename(pfb_slopey)
	if 'pfdist CONUS.slopex.pfb' in line:
		line = line.replace('CONUS.slopex.pfb',os.path.basename(pfb_slopex))
	if 'pfdist CONUS.slopey.pfb' in line:
		line = line.replace('CONUS.slopey.pfb',os.path.basename(pfb_slopey))
	new_content.append(line+'\n')


new_content.append('pfundist '+os.path.basename(pfb_slopex)+'\n')
new_content.append('pfundist '+os.path.basename(pfb_slopey)+'\n')

with open(run_name+'.tcl','w') as fo:
	for line in new_content:
		fo.write(line)

#os.system('tclsh '+run_name+'.tcl')


