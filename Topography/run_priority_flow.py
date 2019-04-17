#!/usr/local/bin/python3

#This code is used to create a new version of Workflow_Example.R in order to run the subset elevation and mask

import os
import sys
import numpy as np
from glob import glob

#try to get name of the dem and mask
try:
	out_name = sys.argv[1]
	out_dem = out_name+'_dem.txt'
	out_river = out_name+'_river.txt'
	out_mask = out_name+'_mask.txt'
except:
	#if out_name is not specified, this code will grab the most recent created text files
	list_text_dem = glob('*_dem.txt')
	list_text_mask = glob('*_mask.txt')
	list_text_river = glob('*_river.txt')
	if not list_text_dem:
		print('error...please check if text files were correctly generated')
		sys.exit()
	out_dem = list_text_dem[-1]
	out_river = list_text_river[-1]
	out_mask = list_text_mask[-1]
	out_name = os.path.basename(out_dem).split('_')[0]

dem_arr = np.loadtxt(out_dem)
ncol = dem_arr.shape[1] #get number of column of the text file
nrow = dem_arr.shape[0] #get number of row of the text file

with open('Workflow_Example4.R','r') as fi:
	content = fi.read()

content = content.split('\n')

#replace R content
new_content = []
for line in content:
	#if 'image.plot' in line:
	#	line = line.replace('image.plot','#image.plot') #stop making plot
	if 'runname=' in line:
		line = line.replace('Test',out_name) #define new runname
	if 'dem_test.txt' in line:
		line = line.replace('dem_test.txt',out_dem) #open the newly subsetted dem file
		line = line.replace('215',str(ncol)) #place the correct number of column
	if 'river_mask_test.txt' in line:
		line = line.replace('river_mask_test.txt',out_river) #open the newly subsetted dem file
		line = line.replace('215',str(ncol)) #place the correct number of column
	if 'mask_test.txt' in line:
		line = line.replace('mask_test.txt',out_mask) #open the newly subsetted dem file
		line = line.replace('215',str(ncol)) #place the correct number of column
	new_content.append(line+'\n')

#write to new R code
with open('PriorityFlow_Irregular_Domain.R','w') as fo:
	for line in new_content:
		fo.write(line)

#run the newly create .R code
os.system('chmod 755 PriorityFlow_Irregular_Domain.R')
os.system('Rscript PriorityFlow_Irregular_Domain.R')

#run .tcl conversion file to convert slope .sa files to .pfb files
with open('File_Conversion.tcl','r') as fi:
	tcl_content = fi.read()

tcl_content = tcl_content.split('\n')
new_content = []
for line in tcl_content:
	if '$out_name' in line:
		line = line.replace('$out_name',out_name)
	if '$nrow' in line:
		line = line.replace('$nrow',str(nrow))
	if '$ncol' in line:
		line = line.replace('$ncol',str(ncol))
	new_content.append(line+'\n')

with open('File_Conversion_'+out_name+'.tcl','w') as fo:
	for line in new_content:
		fo.write(line)

os.system('tclsh File_Conversion_'+out_name+'.tcl')

print('done!')