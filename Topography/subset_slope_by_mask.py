#!/usr/local/bin/python3

#required libraries
import gdal, ogr, osr
#import shapely
import numpy as np
import os
import sys

#Read slope file 
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
n = int(max(new_arr.shape)*0.02)

il, jl, iu, ju = max(min(xx)-n,0), max(min(yy)-n,0), max(xx)+n+1,  max(yy)+n+1

try:
	outnamex = sys.argv[3]
except:
	outnamex = os.path.splitext(os.path.basename(mask_file))[0]+\
				'_'+os.path.splitext(os.path.basename(slope_file))[0]

innamex = os.path.splitext(os.path.basename(slope_file))[0]
innamey = innamex.replace('slopex','slopey')
outnamey = outnamex.replace('slopex','slopey')

if os.path.isfile(outnamex):
	os.remove(outnamex)

if os.path.isfile(outnamey):
	os.remove(outnamey)

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





