# Subsetting tools

Subsetting inputs of smaller domains within the Contiguous US (CONUS) for ParFlow simulations. 

## Directories

1. *CLM*: Subset NLDAS forcing data and generate *.dat files for transient simulations.
2. *Clip_Inputs*: Subset required/optional inputs for ParFlow model. Inputs to be subsetted can be topography, fluxes, initial pressure or subsurface indicator.
3. *Create_Subdomain*: Create a solid file (**.pfsol**) of the subsetted domain.
4. *Make_Tcl*: Generate a tcl script (**.tcl**) from:
	a) A template tcl script.
	b) Newly subsetted inputs and domain.
5. *utils*: Helper codes to reproject CONUS1.0 domain into CONUS2.0 domain.

### Prerequisites

To run this project on your local machine, you will need:
1. Since the code will download all the pre-processed files for the CONUS2.0 domain from Cyverse, you will need to:
	* Sign up for Cyverse account
	* Setting up icommands, more detail can be found [here](https://wiki.cyverse.org/wiki/display/DS/Setting+Up+iCommands) 
2. Python3 installed with these specificed packages: numpy, gdal, ogr, osr, mpl_toolkits.basemap and shapely
3. *pfio* tool for reading and writing .pfb files. To download and install this tool please use this [link](https://github.com/hydroframe/tools/tree/master/pfio)

## Usage of the general_subset.py

### Synopsis

```
python3 general_subset.py {shapefile|mask|define_watershed} [-shp_file] [-id] [-out_name] [-dx] [-dz] [-printmask] [-mask_file] [-dir_file] [-outlet_file]
```

### Description

**{shapefile|mask|define_watershed}** (Required) You need to choose one of the three subsetting methods:
	1. Using a shapefile with an ID of the selected feature.
	2. Using a mask file of the domain with value 1 inside the domain and value 0 elsewhere.
	3. Using an output mask file from *Define_Watershed.py*. *Define_Watershed.py* defines define the watershed for a point or set of outlet points based on the flow direction file.

**-shp_file** (Used conjunctionally with *shapefile*) Name of the shapefile in *.shp* format. Please note that the function also requires other files in *.dbf* and *.prj* along with *.shp* file.   

**-id** (Used conjunctionally with *shapefile*) ID of the selected feature within the shapefile.

**-out_name** (Required) Name of the simulation.

**-dx** (Optional) Spatial resolution of the subset file. Default value is 1000. 

**-dz** (Optional) Vertical resolution of the subset file. Default value is 1000. 

**-printmask** (Optional) Output an image of domain position. Default value is 0.

**-mask_file** (Used conjunctionally with *mask*) Name of the mask file. Please note that the mask file must have same extent and projection with the *input_file*.

**-dir_file** (Used conjunctionally with *define_watershed*) Name of the direction file. Please note that the mask file must have same extent and projection with the *input_file*.

**-outlet_file** (Used conjunctionally with *define_watershed*) Name of the outlet file. The file contains coordinates of each point (y x) by each line.

### Examples

Subsetting by mask

```
python3 general_subset.py mask -mask_file ../../mask/InternalLake1_mask.tif -out_name InternalLake1
```

Subsetting by shapefile

```
python3 general_subset.py shapefile -shp_file ../../shp/Regions.shp -id 14 -out_name UpperCO
```

Subsetting by defined watershed

```
python3 general_subset.py define_watershed -dir_file ../../mask/Str5Ep0_direction.tif -outlet_file outlet.txt -out_name Lake2 -printmask 1
```

### Python workflow

A note about workflow for the python file (general_subset.py):
1. Download all the required domain, slope, subsurface, flux files from Cyverse
2. Call subset_domain.py in Create_Subdomain directory to create solid file of the domain
3. Call clip_inputs_.py in Clip_Inputs directory to create to subset: (1) slope, (2) subsurface and (3) fluxes into new domain extent
3. Call generate_tcl.py in Make_Tcl directory to create new tcl scrip from the parking_lot_template.tcl
4. Move all the necessary inputs files into a folder name 'input_files'
5. Run simulation and outputs into a folder name 'run_output'

### Development Team

See the list of [contributors](https://github.com/orgs/hydroframe/people) who participated in this project.