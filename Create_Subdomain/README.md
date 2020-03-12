# Create Subdomain

Create a solid file of masked domain for ParFlow

### Development Team

See the list of [contributors](https://github.com/orgs/hydroframe/people) who participated in this project.

## Usage

### Prerequisites

To run this project on your local machine, you will need:
1. Since the code will download all the pre-processed files for the CONUS2.0 domain from Cyverse, you will need to:
	* Sign up for Cyverse account
	* Setting up icommands, more detail can be found [here](https://wiki.cyverse.org/wiki/display/DS/Setting+Up+iCommands) 
2. Python3 installed with these specificed packages: numpy, gdal, ogr, osr, and shapely
3. *pfio* tool for reading and writing .pfb files. To download and install this tool please use this [link](https://github.com/hydroframe/tools/tree/master/pfio)

### Synopsis

```
python3 subset_domain.py {shapefile|mask|define_watershed} [-shp_file] [-id] [-out_name] [-dx] [-dz] [-printmask] [-printbbox] [-mask_file] [-dir_file] [-outlet_file]
```

### Description

**{shapefile|mask|define_watershed}** (Required) You need to choose one of the three subsetting methods:
	1. Using a shapefile with an ID of the selected feature.
	2. Using a mask file of the domain with value 1 inside the domain and value 0 elsewhere.
	3. Using an output mask file from *Define_Watershed.py*. *Define_Watershed.py* defines define the watershed for a point or set of outlet points based on the flow direction file.

**-shp_file** (Used conjunctionally with *shapefile*) Name of the shapefile in *.shp* format. Please note that the function also requires other files in *.dbf* and *.prj* along with *.shp* file.   

**-id** (Used conjunctionally with *shapefile*) ID of the selected feature within the shapefile.

**-out_name** (Optional for *shapefile* and *mask*, required for *define_watershed*) Name of the output from the function. If not defined when using *shapefile* or *mask*, *out_name* will be given from selected *id* (for *shapefile*) or from *mask name* (for *mask*). When using *define_watershed*, *out_name* needed to be defined.

**-dx** (Optional) Spatial resolution of the subset file. Default value is 1000. 

**-dz** (Optional) Vertical resolution of the subset file. Default value is 1000. 

**-printmask** (Optional) Output an image of domain position. Default value is 0.

**-printbbox** (Optional) Output a file which contains the bounding box of the domain. Default value is 0.

**-mask_file** (Used conjunctionally with *mask*) Name of the mask file. Please note that the mask file must have same extent and projection with the *input_file*.

**-dir_file** (Used conjunctionally with *define_watershed*) Name of the direction file. Please note that the mask file must have same extent and projection with the *input_file*.

**-outlet_file** (Used conjunctionally with *define_watershed*) Name of the outlet file. The file contains coordinates of each point (y x) by each line.

### Examples
Subsetting by mask

```
python3 subset_domain.py mask -mask_file ../../mask/Coast1_mask.tif
```


Subsetting by shapefile
```
python3 subset_domain.py shapefile -shp_file ../../shp/Regions.shp -id 14
```

Subsetting by a delineated watershed
```
python3 subset_domain.py define_watershed -dir_file ../../mask/Str5Ep0_direction.tif -outlet_file outlet.txt -out_name Lake2 -printmask 1
```


### Python workflow

A note about workflow for the python file (subset_domain.py):
1. Download all the required domain rasters from Cyverse
2. Read the mask file
3. Crop all the rasters to extents which only contain the target basin
4. Create borders (back, front, left, right, top, bottom) ascii files. **Notes:** this part adopted from Prof. Condon R codes.
5. Create a solid domain file by calling mask-to-pfsol function

