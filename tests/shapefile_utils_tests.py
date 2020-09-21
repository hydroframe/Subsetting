import os
import unittest
from pathlib import Path
import numpy as np
import parflow.subset.utils.io as file_io_tools
import tests.test_files as test_files
from parflow.subset.rasterizer import ShapefileRasterizer


class ShapefileReprojectCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.conus1_mask_datset = file_io_tools.read_geotiff(test_files.conus1_mask.as_posix())
        cls.conus2_mask_dataset = file_io_tools.read_geotiff(test_files.conus2_mask.as_posix())
        cls.shape_path = Path(test_files.huc10190004.get('shapefile')).parent
        cls.shape_name = Path(test_files.huc10190004.get('shapefile')).stem

    def test_reproject_conus1(self):
        rasterizer = ShapefileRasterizer(self.shape_path, self.shape_name, self.conus1_mask_datset)
        rasterizer.reproject_and_mask()
        subset_mask = rasterizer.subset_mask
        subset_mask.add_bbox_to_mask()
        subset_mask.write_mask_to_tif(filename='testout.tif')
        self.assertIsNone(
            np.testing.assert_array_equal(file_io_tools.read_file(test_files.huc10190004.get('conus1_mask').as_posix()),
                                          file_io_tools.read_file('testout.tif')),
            'Should create a mask from CONUS1 with 1/0s')
        os.remove('testout.tif')

    def test_reproject_conus2(self):
        rasterizer = ShapefileRasterizer(self.shape_path, self.shape_name, self.conus2_mask_dataset)
        rasterizer.reproject_and_mask()
        subset_mask = rasterizer.subset_mask
        subset_mask.add_bbox_to_mask()
        subset_mask.write_mask_to_tif(filename='testout.tif')
        self.assertIsNone(
            np.testing.assert_array_equal(file_io_tools.read_file(test_files.huc10190004.get('conus2_mask').as_posix()),
                                          file_io_tools.read_file('testout.tif')),
            'Should create a mask from CONUS2 with 1/0s')
        os.remove('testout.tif')

    def test_rasterize_no_data_values(self):
        rasterizer = ShapefileRasterizer(self.shape_path, shapefile_name=self.shape_name,
                                         reference_dataset=self.conus1_mask_datset, no_data=-9999999)
        rasterizer.reproject_and_mask()
        rasterizer.subset_mask.write_mask_to_tif(filename='testout.tif')
        self.assertIsNone(
            np.testing.assert_array_equal(file_io_tools.read_file(
                test_files.huc10190004.get('conus1_mask_-9999999').as_posix()),
                file_io_tools.read_file('testout.tif')),
            'Should create a mask from CONUS1 with 1/-9999999')
        os.remove('testout.tif')

    def test_rasterizer_bad_no_data_val1(self):
        with self.assertRaises(Exception):
            ShapefileRasterizer(self.shape_path, shapefile_name=self.shape_name,
                                reference_dataset=self.conus1_mask_datset, no_data=1)

    def test_rasterizer_bad_no_data_val0(self):
        with self.assertRaises(Exception):
            ShapefileRasterizer(self.shape_path, shapefile_name=self.shape_name,
                                reference_dataset=self.conus1_mask_datset, no_data=0)

    def test_rasterizer_print_no_except(self):
        rasterizer = ShapefileRasterizer(self.shape_path, self.shape_name, self.conus2_mask_dataset)
        self.assertIsNone(print(rasterizer))

if __name__ == '__main__':
    unittest.main()
