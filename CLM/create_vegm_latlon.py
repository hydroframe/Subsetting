#!/usr/local/bin/python3

import numpy as np
import gdal, ogr, osr
from pyproj import Proj, transform
#import matplotlib.pyplot as plt

igbp_file = 'naigbpl20.tif'
latlon_file = 'UC.latlon.txt'
out_file = 'drv_vegm.UC.dat'
nx = 637
ny = 903
sand = 0.16
clay = 0.26
color = 2

#open land cover file
lc_ds = gdal.Open(igbp_file)
##get land cover file projection as outProj
outSRS_converter = osr.SpatialReference()  # makes an empty spatial ref object
outSRS_converter.ImportFromWkt(lc_ds.GetProjection())  # populates the spatial ref object with our WKT SRS
outSRS_forPyProj = outSRS_converter.ExportToProj4()  # Exports an SRS ref as a Proj4 string usable by PyProj
outProj = Proj(outSRS_forPyProj)
##get land cover array
lc_arr = lc_ds.ReadAsArray()
##get land cover geo transform
xul, xres, _, yul, _, yres = lc_ds.GetGeoTransform()

#open lat long file
latlon_arr = np.loadtxt(latlon_file)
inProj = Proj(init='epsg:4326')

#get value of land cover for each coordinate
npoints = latlon_arr.shape[0]
## make output matrix
output = np.zeros((npoints,25))
output[:,4] = sand
output[:,5] = clay
output[:,6] = color

#lc_vals = []
for ii in range(npoints):
	in_lat = latlon_arr[ii,0]
	in_lon = latlon_arr[ii,1]
	out_lon, out_lat = transform(inProj,outProj,in_lon,in_lat)
	lc_val = lc_arr[int(float(out_lat - yul)/yres), int(float(out_lon - xul)/xres)]
	#lc_vals.append(lc_val)
	if (ii+1) % nx == 0:
		x = nx
	else:
		x = (ii+1) % nx
	y = ii // nx +1
	## start writing to output matrix
	output[ii,0:4] = x,y,in_lat,in_lon
	output[ii,lc_val+6] = 1


heading="x y lat lon sand clay color fractional coverage of grid, by vegetation class (Must/Should Add to 1.0)"
col_names=['',  '',  '(Deg)', '(Deg)', '(%/100)', '', 'index', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18']

#write to out file
with open(out_file,'w') as fo:
	fo.write(heading+'\n')
	fo.write(' '.join(col_names)+'\n')
	np.savetxt(fo,output,fmt=['%d']*2+['%.5f']*2+['%.2f']*2+['%d']*19)


#Check result <-- seem right!
'''
min_lat_in = np.min(latlon_arr[:,0])
max_lat_in = np.max(latlon_arr[:,0])
min_lon_in = np.min(latlon_arr[:,1])
max_lon_in = np.max(latlon_arr[:,1])

min_lon_out, min_lat_out = transform(inProj,outProj,min_lon_in,min_lat_in)
max_lon_out, max_lat_out = transform(inProj,outProj,max_lon_in,max_lat_in)

fig, (ax1, ax2) = plt.subplots(2)

ax1.imshow(lc_arr[int(float(max_lat_out-yul)/yres):int(float(min_lat_out-yul)/yres),
					int(float(min_lon_out-xul)/xres):int(float(max_lon_out-xul)/xres)])
#ax1.colorbar()

lats = np.reshape(latlon_arr[:,0],(ny,nx))
lons = np.reshape(latlon_arr[:,1],(ny,nx))
vals = np.reshape(np.array(lc_vals),(ny,nx))

ax2.contourf(lons,lats,vals)
#ax2.colorbar()
plt.show()
'''














