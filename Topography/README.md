# Project Title

Subset ParFlow binary files of slope from CONUS2.0 slope files

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run this project on your local machine, you will need:
1. Since the code will download all the pre-processed files for the CONUS2.0 domain from Cyverse, you will need to:
	* Sign up for Cyverse account
	* Setting up icommands, more detail can be found [here](https://wiki.cyverse.org/wiki/display/DS/Setting+Up+iCommands) 
2. Python3 installed with these specificed packages: numpy, gdal, ogr, osr, and shapely
3. Tcl installed
4. Mask file and slope files already downloaded

## Usage

You will run the subset_slope_by_mask.py file in order to subset the CONUS2.0 slope files to the maskfile extent.

subset_slope_by_mask.py takes three inputs:

```
subset_slope_by_mask.py slope_file mask_file output_name
```

### Description

**slope_file:** (Required) Preprocessed slope file from: [avra](/iplant/home/shared/avra/CONUS2.0/Inputs/Topography/Str5Ep0)

**mask_file_** (Required) Test mask file from [avra](/iplant/home/shared/avra/CONUS2.0/Test_Domains/Small_BC_Tests)

**output_name** (Optional) Name of the output slope files. If not specified, name of the mask file and slope file will be concatenated. 

### Examples

```
./subset_slope_by_mask.py Str5Ep0_unsmth.mx0.5.mn5.sec0.up_slopex.pfb Coast1_mask.tif
```


## Workflows

A note about workflow for the python file (subset_slope_by_mask.py):
1. Read the mask file and extract extents.
2. Clip CONUS slope files using **pfgetsubbox** command

## Authors

See the list of [contributors](https://github.com/orgs/hydroframe/people) who participated in this project.

## Acknowledgments

* TCL subset code provided by Prof. Maxwell

### Side note:
The old codes: *File_Conversion.tcl*, *run_priority_flow.py*, *subset_elevation_by_shape.py*, *subset_elevation.sh*, and *Workflow_Example4.R* are still being kept for future use.