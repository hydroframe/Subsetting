# Project Title

Create tcl script from a parkinglot template

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run this project on your local machine, you will need:
1. ParFlow installed
2. Python3 installed with these specificed packages: numpy, gdal, ogr, osr, and shapely

## Usage

You will run the generate_tcl_script.py file in order to copy all the necessary solid and pfb files, generate a new tcl script, and run the simulation.

generate_tcl_script.py takes three inputs:

```
generate_tcl_script.py name_of_the_mask name_of_the_slope name_of_the_template_file name_of_simulation
```

### Description

**name_of_the_mask:** (Required) Name of the mask that you used for subsetting domain and slopes

**name_of_the_slope** (Required) Name of the slope that you used in subsetting slopes. As mentioned in Domain and Topography sub-directories, unless you specify the output name, outputs of the subset files will be named based on the input files.

**name_of_the_template_file** (Required) Name of the input tcl file.

**name_of_simulation** (Optional) Name of the simulation. If not specified, it will be concatenated between **name_of_the_mask:** and **name_of_the_slope**.

### Examples

```
./generate_tcl_script.py Coast1_mask Str5Ep0_unsmth.mx0.5.mn5.sec0.up_slopex CONUS2.0_Parkinglot.tcl Coast1

```


## Workflows

A note about workflow for the python file (generate_tcl_script.py):
1. Check if necessary files exits and copy them into the current working dir.
2. Get dimensions of the solidfile
3. Find how many patch needs to run the simulation
4. Read the tcl template file and generate a new one by making changes in
	a) 'set runname'
	b) 'pfset ComputationalGrid.NX'
	c) 'pfset ComputationalGrid.NY'
	d) 'pfset GeomInput.domaininput.FileName'
	e) 'pfset Geom.domain.Patches'
	f) 'pfset BCPressure.PatchNames'
	g) 'pfset TopoSlopesY.Type'
	h) 'pfset TopoSlopesX.Type'
	i) 'pfset TopoSlopesX.FileName'
	j) 'pfset TopoSlopesY.FileName'

## Authors

See the list of [contributors](https://github.com/orgs/hydroframe/people) who participated in this project.

## Acknowledgments

* Parkinglot simulation code provided by Prof. Maxwell (Mines) and Prof. Condon (UA)

