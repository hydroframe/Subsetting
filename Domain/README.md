# Project Title

Create a solid file of masked domain for ParFlow

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run this project on your local machine, you will need:
1. Since the code will download all the pre-processed files for the CONUS2.0 domain from Cyverse, you will need to:
	* Sign up for Cyverse account
	* Setting up icommands, more detail can be found [here](https://wiki.cyverse.org/wiki/display/DS/Setting+Up+iCommands) 
2. Python3 installed with these specificed packages: numpy, gdal, ogr, osr, and shapely
3. [pf-mask-utilities](https://github.com/smithsg84/pf-mask-utilities) (C++ complier required). The bash file (subset_domain.sh) already has a part to download and install the ParFlow mask utilities.


## Usage

You will run the subset_domain.sh file in order to install pf-mask-utilities, set permission for the python file (subset_domain_by_shape.py), and run the python file.
subset_domain.sh takes three inputs:

subset_domain.sh [-sel_type sel_type] [-s specific] [-out_name out_name]

### Description

**-sel_type sel_type:** (Required) You can choose subsetting a domain by either:
 * (1) stating the basin id (referenced to HydroSHED watershed id): -id
 * or (2) choosing one point which lies inside the basin: -p

**-s specific** (Required) If you choose *-id* for *-sel_type*, you need to specify the basin id. If you choose *-p* for *-sel_type*, please specify the point coordinates (i.e. lat long) in EPSG:4326 projection

**-out_name out_name** (Optional) Name of the output solid file. If not specified, subsetting basin id will be used 

```
Give an example
```
Subsetting by basin id

./subset_domain.sh -sel_type -id -s 14 -out_name Upper_Colorado

```
Give an example
```
Subsetting by point

./subset_domain.sh -sel_type -p -s 32.8 -108.3 -out_name Upper_Colorado


## Python workflow

A note about workflow for the python file (subset_domain_by_shape.py):
1. Download all the required rasters and shapefile from Cyverse
2. Rasterize the shapefile with the same extent and projection of the *conus_1km_PFmask2.tif*
3. Crop all the rasters to extents which only contain the target basin
4. Create borders (back, front, left, right, top, bottom) ascii files. **Notes:** this part adopted from Prof. Condon R codes.
5. Create a solid domain file by calling mask-to-pfsol function

## Authors

See the list of [contributors](https://github.com/orgs/hydroframe/people) who participated in this project.

## Acknowledgments

* ParFlow raster domain files created by Prof. Condon (UA)
* R codes for subsetting and creating patch borders by Prof. Condon (UA)
* ParFlow Mask utilities codes by Mr. Steven Smith (LLNL)

