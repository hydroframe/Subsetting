#!/bin/bash
shopt -s expand_aliases

export PYTHONWARNINGS="ignore:Unverified HTTPS request"
export GDAL_DATA="/Users/castro/miniconda2/envs/parflow/share/gdal"
export PROJ_LIB="/Users/castro/miniconda2/envs/parflow/share/proj"
alias py="/Users/castro/miniconda2/envs/parflow/bin/python"


## clip slopex
#echo clip slopex
#py Clip_Inputs/clip_inputs.py \
#    -i test/slopex.tif \
#    -out_name output_slopex \
#    -pfmask test/conus_1km_PFmask2.tif \
#    shapefile \
#    -shp_file test/shpfile/watershed.shp \
#    -id 1 2 3 \
#    -att ID
#
## clip slopey
##echo clip slopey
#py Clip_Inputs/clip_inputs.py \
#    -i test/slopey.tif \
#    -out_name output_slopey \
#    -pfmask test/conus_1km_PFmask2.tif \
#    shapefile \
#    -shp_file test/shpfile/watershed.shp \
#    -id 1 2 3 \
#    -att ID 
#
## clip pressure
#echo clip pressure
#py Clip_Inputs/clip_inputs.py \
#    -i test/press.in.tif \
#    -out_name output_press \
#    -pfmask test/conus_1km_PFmask2.tif \
#    shapefile \
#    -shp_file test/shpfile/watershed.shp \
#    -id 1 2 3 \
#    -att ID 
#
## clip grid3d
#echo clip grid3d
#py Clip_Inputs/clip_inputs.py \
#    -i test/grid3d.v3.tif \
#    -out_name output_grid \
#    -pfmask test/conus_1km_PFmask2.tif \
#    shapefile \
#    -shp_file test/shpfile/watershed.shp \
#    -id 1 2 3 \
#    -att ID 
#
#echo subset domain
#py Create_Subdomain/subset_domain.py \
#    -out_name subset.pfsol \
#    -pfmask test/conus_1km_PFmask2.tif \
#    -pflakesmask test/conus_1km_PFmask_selectLakesmask.tif \
#    -pflakesborder test/conus_1km_PFmask_selectLakesborder.tif \
#    -pfbordertype test/1km_PF_BorderCellType.tif \
#    -pfsinks test/conus_1km_PFmask_manualsinks.tif shapefile \
#    -shp_file test/shpfile/watershed.shp \
#    -id 1 2 3 \
#    -att ID
#
echo domain extract
py CLM/domain_extract_latlon.py \
    -shp_file shpfile/watershed.shp \
    -id 1 2 3 \
    -att ID \ 
    -pfmask test/conus_1km_PFmask2.tif \
    -out_name latlon.txt 








