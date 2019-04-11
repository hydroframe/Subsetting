# Project Title

Create ParFlow binary files of topography from elevation

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run this project on your local machine, you will need:
1. Since the code will download all the pre-processed files for the CONUS2.0 domain from Cyverse, you will need to:
	* Sign up for Cyverse account
	* Setting up icommands, more detail can be found [here](https://wiki.cyverse.org/wiki/display/DS/Setting+Up+iCommands) 
2. Python3 installed with these specificed packages: numpy, gdal, ogr, osr, and shapely
3. R installed with a required package: fields
4. Tcl installed

## Usage

You will run the subset_elevation.sh file in order to download PriorityFlow codes from Prof. Condon github, set permission, and run the python files (subset_elevation_by_shape.py and run_priority_flow.py).

subset_elevation.sh takes three inputs:

```
subset_elevation.sh [-sel_type sel_type] [-s specific] [-out_name out_name]
```

### Description

**-sel_type sel_type:** (Required) You can choose subsetting a domain by either:
 * (1) stating the basin id (referenced to HydroSHED watershed id): -id
 * or (2) choosing one point which lies inside the basin: -p

**-s specific** (Required) If you choose *-id* for *-sel_type*, you need to specify the basin id. If you choose *-p* for *-sel_type*, please specify the point coordinates (i.e. lat long) in EPSG:4326 projection

**-out_name out_name** (Optional) Name of the output solid file. If not specified, subsetting basin id will be used 

### Examples

Subsetting by basin id

```
./subset_elevation.sh -sel_type -id -s 14 -out_name Upper_Colorado
```

Subsetting by point

```
./subset_elevation.sh -sel_type -p -s 32.8 -108.3 -out_name Upper_Colorado
```


## Workflows

### PriorityFlow:
*PriorityFlow* is a toolkit for topographic processing for hydrologic models, more detail can be found in [Condon and Maxwell 2019](https://doi.org/10.1016/j.cageo.2019.01.020) and in Prof. Condon PriorityFlow Github [repository](https://github.com/lecondon/PriorityFlow).

### Python workflow #1:
A note about workflow for the python file (subset_elevation_by_shape.py):
1. Download all the preprocessed elevation raster *HSProc_Stream5LakesSinks_EP0.1_dem.tif* and shapefile from Cyverse
2. Rasterize the shapefile with the same extent and projection of the *HSProc_Stream5LakesSinks_EP0.1_dem.tif*
3. Crop all the rasters to extents which only contain the target basin
4. Create dem and mask text files

### Python workflow #2:
A note about workflow for the python file (run_priority_flow.py):
1. Make change to the Workflow_Examole2.R (Irregular domain with no river network) in order to take newly subsetted dem and mask text files as inputs
2. Run the new R code (*PriorityFlow_Irregular_Domain.R*)
3. Make change to the File_Conversion.tcl to take newly created .sa slope files and output them in .pfb format

## Authors

See the list of [contributors](https://github.com/orgs/hydroframe/people) who participated in this project.

## Acknowledgments

* ParFlow elevation raster created by Prof. Condon (UA)
* R codes for DEM processing by Prof. Condon (UA)

