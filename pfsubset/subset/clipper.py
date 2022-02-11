"""Classes for clipping gridded inputs"""
import os
import logging
from abc import ABC, abstractmethod
import numpy as np
from numpy.core.fromnumeric import clip
from numpy.core.numeric import full
import numpy.ma as ma
from parflow.tools.io import read_clm
from pfsubset.subset import TIF_NO_DATA_VALUE_OUT as NO_DATA
from pfsubset.subset.utils import io as file_io_tools
from pfsubset.subset.mask import SubsetMask
from parflowio.pyParflowio import PFData

class Clipper(ABC):

    """Abstract Clipper Class"""

    @abstractmethod
    def subset(self, data):
        """Clip the data

        Parameters
        ----------
        data : pfb files or ndarray
            3d array of data to be clipped

        Returns
        -------
        ndarray
            clipped portion of pfb files or data_array
        """
        pass


class BoxClipper(Clipper):
    """Clip a rectangular data region specified by a bounding box
    """
    def __repr__(self):
        if self.box:
            return f"{self.__class__.__name__}(x_0:{self.x_0}, x_end:{self.x_end}, y_0:{self.y_0}, y_end:{self.y_end}, nx:{self.nx}, ny:{self.ny}, nz:{self.nz!r}, " \
                f"ref_array:{self.ref_array!r}, padding:{self.padding!r}, no_data:{self.no_data!r}"
        else:
            return f"no_data:{self.no_data!r}"

        
    def __init__(self,x=1, y=1, z=1, nx=None, ny=None, nz=None, padding=(0, 0, 0, 0), no_data=NO_DATA):
        """

        Parameters
        ----------
        x : int, optional
            the starting x value (1 based index)
        y : int, optional
            the starting y value (1 based index)
        z : int, optional
            the starting z value (1 based index)
        nx : int, optional
            the number of x cells to clip (default extents of 'ref_array')
        ny : int, optional
            the number of y cells to clip (default extents of 'ref_array')
        nz : int, optional
            the number of z cells to clip (default extents of 'ref_array')
        padding : tuple, optional
            no_data padding to add around data. specified clockwise from top. (top,right,bot,left)
        no_data: int
            the no data value to use (default package NO_DATA_VALUE)

        Returns
        -------
        BoxClipper

        Raises
        -------
        Exception
            Invalid dimensions will raise an exception
        """
  
        self.padding = padding
        self.no_data = no_data

        self.box = False

        self.x = x
        self.y = y
        self.z = z
        self.nx = nx
        self.ny = ny
        self.nz = nz
    
        self.padding = padding
        

        if  self.nx is not None:
            if nx < 1 or ny < 1  or x < 1 or y < 1:
                raise Exception("Error: invalid dimension, x,y,z nx, ny, nz must be >=1")
                
            self.update_bbox(x, y, z, nx, ny, nz, padding)

    def update_bbox(self, x=None, y=None, z=None, nx=None, ny=None, nz=None, padding=(0, 0, 0, 0)):
        """update the x,y,z, nx, ny, nz and padding values

        Parameters
        ----------
        x : int, optional
             Starting X value (Default value = None)
        y : int, optional
             Staring Y value (Default value = None)
        z : int, optional
             Starting Z value (Default value = None)
        nx : int, optional
             number of X cells to clip (Default value = None)
        ny : int, optional
             number of Y cells to clip (Default value = None)
        nz : int, optional
             number of Z cells to clip (Default value = None)
        padding : tuple, optional
             number of no_data cells to add around box, clockwise starting from top (CSS style)
             (default value = (0,0,0,0)

        Returns
        -------
        None
        """
        if nx is not None:
            self.nx = nx
        if ny is not None:
            self.ny = ny
        if nz is not None:
            self.nz = nz
        if x is not None:
            self.x_0 = x - 1
            self.x_end = self.x_0 + nx
        if z is not None and nz is not None:
            self.z_0 = z - 1
            self.z_end = self.z_0 + nz
        if y is not None:
            self.y_0 = y - 1
            self.y_end = self.y_0 + ny
        self.padding = padding
        self.box = True 

    

    def get_latest_data_array(self):
        return(self.return_array)
    
    def write_to_file(self,out_dir,file_name, ref_proj, no_data):
        if file_name.endswith('.pfb'):
            file_io_tools.write_pfb(self.return_array, os.path.join(out_dir, f'{file_name}'))
        if file_name.endswith('.tif') and self.clipped_geom is not None and ref_proj is not None:
            file_io_tools.write_array_to_geotiff(os.path.join(out_dir, f'{file_name}'),
                                                self.return_array, self.clipped_geom, ref_proj, no_data=no_data)

    def subset(self, data=None):
        """ Clip the data_array to the region specified by the bounding box

        Parameters
        ----------
        data : pfb file, tif file, ndarray
            pfb file or tif file or a 3d array of data to clip

        Returns
        -------
        ret_array : ndarray
            the clipped `data_array`

        """
        data_file = None
        data_array = None
        if (not isinstance(data, np.ndarray)):
            data_file = str(data)

            if data_file.endswith('.tif'):
                data_array = file_io_tools.read_file(data_file)
    
            elif data_file.endswith('.pfb'):  # pfsubset binary file
                pfdata = PFData((data_file))
                pfdata.loadHeader()
                if self.box:
                    pfdata.loadClipOfData(clip_x=self.x_0,clip_y=self.y_0,extent_x=self.nx,extent_y=self.ny)
                else:
                    pfdata.loadData()
                data_array = pfdata.moveDataArray()
                pfdata.close()
        else:
            data_array = data
        

        if(not self.box):
            return data_array, None, None, None
            #self.update_bbox(self.x, self.y, self.z, self.nx, self.ny, self.nz, self.padding)
    
        if any(self.padding):
            #print("This case (padding) is not supported at the moment\n")
            # create a full dimensioned array of no_data_values
            ret_array = np.full(shape=(data_array.shape[0],
                            self.ny + self.padding[0] + self.padding[2],
                            self.nx + self.padding[1] + self.padding[3]), fill_value=self.no_data, dtype=data_array.dtype)
            # assign values from the data_array into the return array, mind the padding
            ret_array[:, self.padding[2]:self.ny + self.padding[2], self.padding[3]:self.nx + self.padding[3]] = data_array[:, self.y_0:self.y_end,self.x_0:self.x_end]
        else:
            ret_array = data_array[:,self.y_0:self.y_end, self.x_0:self.x_end]
            return ret_array , None, None, None
        self.return_array = ret_array
        return ret_array, None, None, None


class MaskClipper(Clipper):
    """Clip an irregular data region specified by a mask"""

    def __repr__(self):
        return f"{self.__class__.__name__}(subset_mask:{self.subset_mask!r}, bbox:{self.bbox!r}, " \
               f"clipped_geom:{self.clipped_geom!r}, clipped_mask:{self.clipped_mask!r}"

    def __init__(self, subset_mask, no_data_threshold=NO_DATA):
        """Assumes input mask_array has 1's written to valid data, 0's for bounding box,
             and <=no_data_threshold for no_data val, no_data_threshold must be < 0 or bounding box will
             not be identifiable

        Parameters
        ----------
        subset_mask : SubsetMask or Mask file
            instantiated mask object
        no_data_threshold : int
            upper bound to which all values are no_data

        Returns
        -------
        MaskClipper
        """
    
        if type(subset_mask) is str:
            self.subset_mask = SubsetMask(subset_mask)
        else:
            self.subset_mask = subset_mask


        min_y, max_y, min_x, max_x = self.subset_mask.bbox_edges
        self.bbox = [min_y, max_y + 1, min_x, max_x + 1]
        self.clipped_geom = self.subset_mask.calculate_new_geom(min_x, min_y,
                                                                self.subset_mask.mask_tif.GetGeoTransform())
        self.clipped_mask = (self.subset_mask.bbox_mask.filled(fill_value=no_data_threshold) == 1)[:,
                            self.bbox[0]:self.bbox[1],
                            self.bbox[2]:self.bbox[3]]

    def get_latest_data_array(self):
        return(self.return_arr)
    
        
    def write_to_file(self,out_dir,file_name, ref_proj, no_data):
        if file_name.endswith('.pfb'):
            file_io_tools.write_pfb(self.return_array, os.path.join(out_dir, f'{file_name}'))
        if file_name.endswith('.tif') and self.clipped_geom is not None and ref_proj is not None:
            file_io_tools.write_array_to_geotiff(os.path.join(out_dir, f'{file_name}'),
                                                self.return_array, self.clipped_geom, ref_proj, no_data=no_data)

    def subset(self, data, no_data=NO_DATA, crop_inner=1):
        """subset the data from data_array in the shape and extents of the clipper's clipped subset_mask

        Parameters
        ----------
        data: pfb or tif file or numpy.ndarray
            pfb file or 3d array of data to subset
        no_data : int
            no data value for outputs (Default = NO_DATA)
        crop_inner : int
            crop the data to the bbox(0) or the mask(1) (default = 1)

        Returns
        -------
        return_arr : numpy.ndarray
            the subset data as a 3d array
        clipped_geom : list
            gdal geometry for the clip
        clipped_mask : numpy.ndarray
            mask of 1/0 showing clipped area
        bounding box : tuple
            x, y, nx, ny values indicating region clipped

        """
        data_file = None
        data_array = None
        if (type(data) is str):
            data_file = data
            if data_file.endswith('.tif'):
                data_array = file_io_tools.read_file(data_file)
    
            elif data_file.endswith('.pfb'):  # pfsubset binary file
                pfdata = PFData((data_file))
                pfdata.loadHeader()

                x = self.bbox[2]
                nx = self.bbox[3] - self.bbox[2]

                y = self.bbox[0]
                ny = self.bbox[1] - self.bbox[0]

                pfdata.loadClipOfData(int(x), int(y), int(nx), int(ny))
                data_array = pfdata.viewDataArray()

                #return data_array, self.clipped_geom, self.clipped_mask, self.subset_mask.get_human_bbox()
                pfdata.close()
                del pfdata
                full_mask = self.subset_mask.bbox_mask.mask
                full_mask = full_mask[:,self.bbox[0]: self.bbox[1],self.bbox[2]: self.bbox[3]]  #data_array.shape[1], data_array.shape[2]))

                clip_mask = ~self.clipped_mask
                # Handle multi-layered files, such as subsurface or forcings
                if data_array.shape[0] > 1:
                    full_mask = np.broadcast_to(full_mask, data_array.shape)
                    clip_mask = np.broadcast_to(clip_mask, (data_array.shape[0], clip_mask.shape[1], clip_mask.shape[2]))
                    logging.info(f'clipper: broadcast subset_mask to match input data z layers: {data_array.shape[0]}')
                # full_dim_mask the input data using numpy masked array module (True=InvalidData, False=ValidData)
                masked_data = ma.masked_array(data=data_array, mask=full_mask)

                return_arr = ma.masked_array(masked_data.filled(fill_value=no_data),
                                         mask=clip_mask).filled(fill_value=no_data)

                return return_arr, self.clipped_geom, self.clipped_mask, self.subset_mask.get_human_bbox()

            else:
                raise Exception("Error: only pfb or tif files can be used directly with clippers")
        else:
            data_array = data

        full_mask = self.subset_mask.bbox_mask.mask
        clip_mask = ~self.clipped_mask
        # Handle multi-layered files, such as subsurface or forcings
        if data_array.shape[0] > 1:
            full_mask = np.broadcast_to(full_mask, data_array.shape)
            clip_mask = np.broadcast_to(clip_mask, (data_array.shape[0], clip_mask.shape[1], clip_mask.shape[2]))
            logging.info(f'clipper: broadcast subset_mask to match input data z layers: {data_array.shape[0]}')
        # full_dim_mask the input data using numpy masked array module (True=InvalidData, False=ValidData)
        masked_data = ma.masked_array(data=data_array, mask=full_mask)

        if crop_inner:
            # return an array that includes all of the z data, and x and y no_data outside of the full_dim_mask area
            return_arr = ma.masked_array(masked_data[:,
                                         self.bbox[0]: self.bbox[1],
                                         self.bbox[2]: self.bbox[3]].filled(fill_value=no_data),
                                         mask=clip_mask).filled(fill_value=no_data)
            # logging.info(f'clipped data with (z,y,x) shape {data_array.shape} to {return_arr.shape} '
            #              f'using bbox (top, bot, left, right) {self.printable_bbox}')
        else:
            # return an array that includes all of the z data, and x and y inside the bounding box
            return_arr = ma.masked_array(masked_data[:,
                                         self.bbox[0]: self.bbox[1],
                                         self.bbox[2]: self.bbox[3]],
                                         mask=np.zeros(clip_mask.shape)).filled()
            # logging.info(f'clipped data with (z,y,x) shape {data_array.shape} to {return_arr.shape} '
            #              f'using bbox (top, bot, left, right) {self.printable_bbox}')

        self.return_arr = return_arr
        return self.return_arr, self.clipped_geom, self.clipped_mask, self.subset_mask.get_human_bbox()


class ClmClipper:
    """Specialized clipper for CLM input files"""

    def __init__(self, bbox):
        """Clips CLM datafiles lat/lon and land cover

        Parameters
        ----------
        bbox : BBox
            BBox object describing the mask bounds
        """
        self.bbox = bbox.get_human_bbox()
        padded_extents = bbox.get_padded_extents()
        self.clipper = BoxClipper(x=padded_extents[2]+1, y=padded_extents[0]+1, nx=padded_extents[3]-padded_extents[2], ny=padded_extents[1]-padded_extents[0],
                                  nz=1)

    def clip_latlon(self, lat_lon_file):
        """Clip the domain lat/lon data to the bounding box of the mask

        Parameters
        ----------
        lat_lon_file : str
            lat/lon data for the domain

        Returns
        -------
        sa_formatted : numpy.ndarray
            the clipped data formatted in a 1d output arrray for writing to a .sa file
        clipped_data : numpy.ndarray
            the raw clipped data (3d) read from 'lat_lon_file'

        """
        data = file_io_tools.read_file(lat_lon_file)
        clipped_data, _, clipped_mask, bbox = self.clipper.subset(data=data)
        #sa_formatted = np.flip(clipped_data, axis=1).flatten()
        sa_formatted = clipped_data.flatten()
        return sa_formatted, clipped_data

    def clip_land_cover(self, lat_lon_array, land_cover_file):
        """Clip the domain land cover data to the bounding box of the mask

        Parameters
        ----------
        lat_lon_array : numpy.ndarray
            clipped lat/lon tuple data for the masked area
        land_cover_file : str
            land cover file for the domain

        Returns
        -------
        sa_formatted : numpy.ndarray
            the clipped data formatted in a 1d output arrray for writing to a .sa file
        output : numpy.ndarray
            vegm formatted representation of the data (2d)
        """
        lat_lon_proper = np.char.split(lat_lon_array.astype(str), ' ')
        data = file_io_tools.read_file(land_cover_file)
        clipped_data, _, clipped_mask, bbox = self.clipper.subset(data=data)
        #sa_formatted = np.flip(clipped_data, axis=1).flatten()
        sa_formatted = clipped_data.flatten()
        sand = 0.16
        clay = 0.26
        color = 2
        # get value of land cover for each coordinate
        npoints = sa_formatted.shape[0]
        # make output matrix
        output = np.zeros((npoints, 25))
        output[:, 4] = sand
        output[:, 5] = clay
        output[:, 6] = color
        # assign x values, looping from 1 to x extent, holding y constant
        output[:, 0] = list(range(1, clipped_data.shape[2] + 1)) * clipped_data.shape[1]
        # assign y values, repeating each y value from 1 to y extent for every x
        output[:, 1] = np.repeat(range(1, clipped_data.shape[1] + 1), clipped_data.shape[2])
        # assign lat values
        output[:, 2] = [latlon[0] for latlon in lat_lon_proper]
        # assign lon values
        output[:, 3] = [latlon[1] for latlon in lat_lon_proper]
        cols = [int(i) + 6 for i in sa_formatted]
        rows = list(range(npoints))
        output[rows, cols] = 1
        return sa_formatted, output

    def clip_vegm(self, vegm_file):
        """Clip a vegm .dat file to the bounding box of the mask

        Parameters
        ----------
        vegm_file : str
            Path to the vegm .dat file to clip

        Returns
        -------
        output : numpy.ndarray
            vegm formatted representation of the data (2d)
        """
        vegm_data = read_clm(vegm_file, type='vegm')

        # pfsubset returns an ndarray from read_clm which is (x, y) in dimensions
        # transpose so that we're dealing with (z, y, x) dimensions
        vegm_data = np.transpose(vegm_data)

        # We will be subsetting exactly 23 z values corresponding to all vegm columns
        # Save old values from self.clipper so we can reapply them later
        _z_0 = self.clipper.z_0
        _nz = self.clipper.nz

        self.clipper.update_bbox(z=1, nz=23)
        clipped_data, _, _, _ = self.clipper.subset(data=vegm_data)

        # Reapply old values - note that the z-parameter in update_bbox is 1-indexed, so _z_0 + 1
        self.clipper.update_bbox(z=_z_0+1, nz=_nz)

        # The following lines of code have to do with generating the first 2 (x, y) columns needed
        # for the vegm .dat file
        ny, nx = clipped_data.shape[1:]
        # y-indices, then x-indices
        indices = np.indices((ny, nx)) + 1
        # x-indices, then y-indices
        indices = indices[::-1, :, :]
        # prepend x-indices, y-indices to transform shape (23, y, x) => (25, y, x)
        clipped_data = np.vstack([indices, clipped_data])

        # Convert to (x, y, 25) in shape and merge the first 2 axes to create a (x*y, 25) ndarray.
        # The result corresponds to the parameters in other methods of this class:
        #    the 2D output 'output' of clip_land_cover
        #    the 2D input 'land_cover_data' of write_land_cover
        clipped_data = clipped_data.transpose().reshape(-1, 25)

        return clipped_data

    def write_land_cover(self, land_cover_data, out_file):
        """Write the land cover file in vegm format

        Parameters
        ----------
        land_cover_data : ndarray
            formatted vegm data (2d array)
        out_file : str
            path and name to write output

        Returns
        -------
        None
        """
        heading = "x y lat lon sand clay color fractional coverage of grid, by vegetation class (Must/Should Add to " \
                  "1.0) "
        col_names = ['', '', '(Deg)', '(Deg)', '(%/100)', '', 'index', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                     '10', '11',
                     '12', '13', '14', '15', '16', '17', '18']
        header = '\n'.join([heading, ' '.join(col_names)])

        file_io_tools.write_array_to_text_file(out_file=out_file, data=land_cover_data, header=header,
                                               fmt=['%d'] * 2 + ['%.6f'] * 2 + ['%.2f'] * 2 + ['%d'] * 19)

    def write_lat_lon(self, lat_lon_data, out_file, x=0, y=0, z=0):
        """Write the lat/lon data to a ParFlow simple ascii formatted file

        Parameters
        ----------
        lat_lon_data : ndarray
            lat/lon data in formatted 1d array
        out_file : str 
            path and name for output file
        x : int
            size of x dimension (Default value = 0)
        y : int
            size of y dimension (Default value = 0)
        z : int
            size of z dimension (Default value = 0)

        Returns
        -------
        None

        """
        file_io_tools.write_array_to_text_file(out_file=out_file, data=lat_lon_data, fmt='%s', header=f'{x} {y} {z}')
