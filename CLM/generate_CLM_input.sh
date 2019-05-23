#!/bin/sh
#Steps to generate inputs for CLM
# 1. Extract lat long (in degree - epsg:4326) of the domain
# 2. Preprocess NLDAS input
# 3. Prepare drv_clmin (Index file)
# 4. Copy drv_vegp.dat file
# 5. Create clm vegm file (designates the land cover fractions for every cell)
# 6. Distribute the forcing input

###Parse inputs




###Step 1: Extract lat long of the domain
#### python file: domain_extract_latlon.py (need to fix the file to meet Laura's requirement!!!)

python3 domain_extract_latlon.py 14

###Step 2: Preprocess NLDAS input (process in powell server): for more information please refer to the file preproc/get_nldas_workflow.txt
####Step 2.1. Edit the batch.get_nldas file to get specific information
####Step 2.2. Copy the newly created latlon file to preproc folder
####Step 2.3. Run the executable

###Step3: Edit the index file (drv_clmin.dat)

###Step4: Copy the drv_vegp.dat file

###Step5: Create the clm_vegm file
####Python file: create_vegm_latlon.py (still need to fix for standard coding requirement)
####Required inputs:
#	a) the domain latlon file
#	b) igbp land cover file for North America
python3 create_vegm_latlon.py

###Step6: Distribute the forcing input
tclsh Dist_Forcings.tcl