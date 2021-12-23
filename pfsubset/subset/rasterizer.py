"""Classes for converting inputs to gridded masks"""
try:
    from osgeo import gdal
except ImportError:
    import gdal

try:
    from osgeo import ogr
except ImportError:
    import ogr

import os
import logging
from pfsubset.subset import TIF_NO_DATA_VALUE_OUT as NO_DATA
from pfsubset.subset.mask import SubsetMask


class ShapefileRasterizer:
    """Class for converting shapefile to raster for use as mask"""

    def __init__(self, input_path, shapefile_name, reference_dataset, no_data=NO_DATA, output_path='.'):
        """
        Parameters
        ----------
        input_path : str
            path to input files (shapefile set)
        shapefile_name : str
            name of shapefile dataset
        reference_dataset : gdal.dataset
            gdal dataset defining the overall domain
        no_data : int, optional
            value to write for no_data cells (default -999)
        output_path : str, optional
            where to write the outputs (default '.')

        Returns
        -------
        ShapefileRasterizer
        """
        if no_data in [0, 1]:
            raise Exception(f'ShapfileRasterizer: '
                            f'Do not used reserved values 1 or 0 for no_data value: got no_data={no_data}')
        self.shapefile_path = input_path
        self.output_path = output_path
        self.shapefile_name = shapefile_name
        self.ds_ref = reference_dataset
        self.no_data = no_data
        # TODO Handle shape extension using Pathlib.path
        self.full_shapefile_path = os.path.join(self.shapefile_path, '.'.join((self.shapefile_name, 'shp')))
        self.check_shapefile_parts()
        self.subset_mask = None

    def __repr__(self):
        return f"{self.__class__.__name__}(shapefile_path:{self.shapefile_path!r}, " \
               f"shapefile_name:{self.shapefile_name!r}, output_path:{self.output_path!r}, ds_ref:{self.ds_ref!r}, " \
               f"no_data:{self.no_data!r}, full_shapefile_path:{self.full_shapefile_path!r}, " \
               f"subset_mask:{self.subset_mask!r}"

    def check_shapefile_parts(self):
        """verify the required parts of a shapefile are present in the same folder
        logs a warning
        Returns
        -------
        None
        """
        shape_parts = [".".join((self.shapefile_name, ext)) for ext in ['shp', 'dbf', 'prj', 'shx']]
        for shp_component_file in shape_parts:
            if not os.path.isfile(os.path.join(self.shapefile_path, shp_component_file)):
                logging.warning(f'Shapefile path missing {shp_component_file}')

    def reproject_and_mask(self, dtype=gdal.GDT_Int32, no_data=None, attribute_name='OBJECTID', attribute_ids=None):
        """

        Parameters
        ----------
        attribute_ids : list of ints
            list of attribute ID values to select (Default value = None)
        dtype : gdal.datatype
            the datatype to write (Default value = gdal.GDT_Int32)
        no_data : str
            no_data value to use (Default value = None)
        attribute_name : str
            field in the shapefile to trace (Default value = 'OBJECTID')

        Returns
        -------
        str
            path (virtual mem) to the reprojected full_dim_mask

        """
        if attribute_ids is None:
            attribute_ids = [1]
        if no_data is None:
            no_data = self.no_data
        geom_ref = self.ds_ref.GetGeoTransform()
        tif_path = f'/vsimem/{self.shapefile_name}.tif'
        target_ds = gdal.GetDriverByName('GTiff').Create(tif_path,
                                                         self.ds_ref.RasterXSize,
                                                         self.ds_ref.RasterYSize,
                                                         1, dtype)
        target_ds.SetProjection(self.ds_ref.GetProjection())
        target_ds.SetGeoTransform(geom_ref)
        target_ds.GetRasterBand(1).SetNoDataValue(no_data)
        target_ds.GetRasterBand(1).Fill(no_data)
        # shapefile
        shp_source = ogr.Open(self.full_shapefile_path)
        shp_layer = shp_source.GetLayer()
        # TODO: How to detect if the shape geometries extend beyond our reference bounds?
        # Filter by the shapefile attribute IDs we want
        shp_layer.SetAttributeFilter(f'{attribute_name} in ({",".join([str(i) for i in attribute_ids])})')
        # Rasterize layer
        rtn_code = gdal.RasterizeLayer(target_ds, [1], shp_layer, burn_values=[1])
        if rtn_code == 0:
            target_ds.FlushCache()
            logging.info(f'reprojected shapefile from {str(shp_layer.GetSpatialRef()).replace(chr(10), "")} '
                         f'with extents {shp_layer.GetExtent()} '
                         f'to {self.ds_ref.GetProjectionRef()} with transform {self.ds_ref.GetGeoTransform()}')
        else:
            msg = f'error rasterizing layer: {shp_layer}, gdal returned non-zero value: {rtn_code}'
            logging.exception(msg)
            raise Exception(msg)
        self.subset_mask = SubsetMask(tif_path)
        return tif_path

    def rasterize_shapefile_to_disk(self, out_dir=None, out_name=None, padding=(0, 0, 0, 0), attribute_name='OBJECTID',
                                    attribute_ids=None):
        """rasterize a shapefile to disk in the projection and extents of the reference dataset

        Parameters
        ----------
        out_dir : str
            directory to write outputs (Default value = None)
        out_name : str
            filename for outputs (Default value = None)
        padding : tuple
            optional padding to add 0's around full_dim_mask (Default value = (0,0,0,0))
        attribute_name : str
            optional name of shapefile attribute to select on (Default value = 'OBJECTID')
        attribute_ids : list
            optional list of attribute ids in shapefile to select for full_dim_mask (Default value = None)

        Returns
        -------
        ndarray
            3d array with no_data to extents, 0 in bounding box, 1 in full_dim_mask region

        """
        if attribute_ids is None:
            attribute_ids = [1]
        if out_name is None:
            out_name = f'{self.shapefile_name}.tif'
        if out_dir is None:
            out_dir = self.output_path
        self.reproject_and_mask(attribute_ids=attribute_ids, attribute_name=attribute_name)
        self.subset_mask.add_bbox_to_mask(padding=padding)
        self.subset_mask.write_mask_to_tif(filename=os.path.join(out_dir, out_name))
        self.subset_mask.write_bbox(os.path.join(out_dir, 'bbox.txt'))
        return self.subset_mask.mask_array
