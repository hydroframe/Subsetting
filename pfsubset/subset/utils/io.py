"""Common disk I/O operations"""

import os
import sys
import logging
from pathlib import Path
import pandas as pd
import numpy as np
try:
    from osgeo import gdal
except ImportError:
    import gdal
from pfsubset.subset import TIF_NO_DATA_VALUE_OUT as NO_DATA
from parflowio.pyParflowio import PFData
from pfsubset.subset.bbox import BBox


def read_file(infile,min_x=None,min_y=None,nx=None,ny=None):
    """read an input file and return a 3d numpy array

    Parameters
    ----------
    infile : str
        file to open (.pfb, .sa, .tif, .tiff)

    Returns
    -------
    res_arr : ndarray
        a 3d numpy array with data from file in (z,y,x) format with y axis 0 at bottom

    """
    infile_path = Path(infile)
    # get extension
    ext = infile_path.suffix
    file_string_path = os.fspath(infile_path)
    if ext in ['.tif', '.tiff']:
        res_arr = gdal.Open(file_string_path).ReadAsArray()
        if len(res_arr.shape) == 2:
            res_arr = res_arr[np.newaxis, ...]
        # flip y axis so tiff aligns with PFB native alignment
        res_arr = np.flip(res_arr, axis=1)
    elif ext == '.sa':  # parflow ascii file
        with open(file_string_path, 'r') as fi:
            header = fi.readline()
        nx, ny, nz = [int(x) for x in header.strip().split(' ')]
        arr = pd.read_csv(file_string_path, skiprows=1, header=None).values
        res_arr = np.reshape(arr, (nz, ny, nx))[:, :, :]
    elif ext == '.pfb':  # parflow binary file
        pfdata = PFData(file_string_path)
        pfdata.loadHeader()
        if not (min_x is None):
            print("Attemping to load clip\n");
            pfdata.loadClipOfData(min_x,min_y,nx,ny)
        else:
            pfdata.loadData()
        res_arr = pfdata.moveDataArray()
        pfdata.close()
        del pfdata
    else:
        raise ValueError('can not read file type ' + ext)

    return res_arr


def read_geotiff(infile):
    """wrapper for reading geotifs with gdal

    Parameters
    ----------
    infile : str
        the geotif to open

    Returns
    -------
    dataset
        gdal dataset object

    """
    file_path = Path(infile)

    return gdal.Open(os.fspath(file_path))


def write_pfb(data, outfile, x0=0, y0=0, z0=0, dx=1000, dz=1000):
    """Write a 3d numpy array to a PFB output file

    Parameters
    ----------
    data : ndarray
        3d numpy data array to write to pfb !(x,y,z)!
    outfile : str
        filename and path to write output
    x0 : int, optional
        initial x location (Default value = 0)
    y0 : int, optional
        initial y location (Default value = 0)
    z0 : int, optional
        initial z location (Default value = 0)
    dx : int, optional
        horizontal resolution (Default value = 1000)
    dz : int, optional
        vertical resolution (Default value = 1000)

    Returns
    -------
    None

    """
    logging.info(f'wrote pfb file {outfile}, (z,y,x)={data.shape}')
    pf_data = PFData()
    pf_data.setDataArray(data)
    pf_data.setDX(dx)
    pf_data.setDY(dx)
    pf_data.setDZ(dz)
    pf_data.setX(x0)
    pf_data.setY(y0)
    pf_data.setZ(z0)
    pf_data.writeFile(outfile)
    del pf_data


def write_bbox(bbox, outfile):
    """Write bounding box values to tab separated text file

    Parameters
    ----------
    bbox : list of ints
        array of bounding box values [x1, y1, nx, ny]
    outfile :
        where to write the file

    Returns
    -------
    None

    """
    logging.info(f'wrote bbox file {outfile}, {bbox}')
    with open(outfile, 'w') as fp:
        fp.write('x1\ty1\tnx\tny\n')
        fp.write('\t'.join('%d' % x for x in bbox))


def read_bbox(bbox_file):
    """Parse a tab separated bounding box text file and return the array of values as integers

    Parameters
    ----------
    bbox_file : str
        the file to read

    Returns
    -------
    list of ints
        an array of integers representing the bounding box [x1, y1, nx, ny]

    """
    with open(bbox_file, 'r') as bbox:
        lines = bbox.readlines()
        bbox_vals = [int(s) for s in lines[1].split('\t')]
        bbox = BBox(x_1=bbox_vals[0], y_1=bbox_vals[1], nx=bbox_vals[2], ny=bbox_vals[3])
        return bbox.get_human_bbox()


def write_array_to_geotiff(out_raster_path, data, geo_transform, projection, dtype=gdal.GDT_Float64, no_data=NO_DATA):
    """write a numpy array to a geotiff

    Parameters
    ----------
    out_raster_path : str
        where to write the output file
    data : ndarray
        3d array of data to write
    geo_transform : list of ints
        gdal formatted geoTransform to use for the geoTif
    projection : str
        srs wkt formatted Projection to use for the geoTif
    dtype : gdal.datatype, optional
        gdal datatype to use for the geoTif (Default value = gdal.GDT_Float64)
    no_data : int, optional
        no data value to encode in the geoTif (Default value = NO_DATA)

    Returns
    -------
    None

    """
    # flip the tif y axis back to tif standard (Tif 0's start at top left, PFB 0's at bottom left)
    data = np.flip(data, axis=1)
    data = np.ascontiguousarray(data, dtype=np.float64)
    driver = gdal.GetDriverByName('GTiff')
    no_bands, rows, cols = data.shape
    data_set = driver.Create(out_raster_path, xsize=cols, ysize=rows, bands=no_bands, eType=dtype,
                             options=['COMPRESS=LZW', 'NUM_THREADS=ALL_CPUS'])
    data_set.SetGeoTransform(geo_transform)
    data_set.SetProjection(projection)
    for i, image in enumerate(data, 1):
        data_set.GetRasterBand(i).WriteArray(image)
        data_set.GetRasterBand(i).SetNoDataValue(no_data)
    logging.info(f'wrote geotif {out_raster_path}, (bands,rows,cols)=({no_bands}, {rows}, {cols})')
    # noinspection PyUnusedLocal
    data_set = None


def write_array_to_simple_ascii(data, out_file):
    """write an array to ParFlow simple ascii (.sa) format

    Parameters
    ----------
    data : ndarray
        3d numpy array of data
    out_file : str
        filename to write data

    Returns
    -------
    None

    """
    write_array_to_text_file(out_file=out_file, data=data.flatten(), fmt='%s',
                             header=f'{data.shape[2]} {data.shape[1]} {data.shape[0]}')


def write_array_to_text_file(data, out_file, header, fmt, delimiter=' ', comments=''):
    """write a flattened array to text output file

    Parameters
    ----------
    data : ndarray
        the 1d numpy array of data to write
    out_file : str
        where to write the data
    header : str
        the header to use in the file
    fmt : str
        the python string format to use for printing each element of the array
    delimiter : str, optional
        the delimiter character to use when writing the elements (optional) (Default value = ' ')
    comments : str, optional
        the comment character to write for header (optional) (Default value = '')

    Returns
    -------
    None

    """
    np.savetxt(fname=out_file, X=data, delimiter=delimiter, comments=comments, header=header, fmt=fmt)
