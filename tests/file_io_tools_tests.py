import unittest
import os
import numpy as np
import gdal
from osgeo import osr
import parflow.subset.utils.io as file_io_tools
import tests.test_files as test_files


class FileIOToolBasicTestCase(unittest.TestCase):

    def test_read_sa(self):
        results = file_io_tools.read_file(test_files.forcings_sa)
        self.assertEqual((24, 41, 41), results.shape)
        self.assertEqual(290.4338, results[0, 0, 0])
        self.assertEqual(292.9594, results[-1, 40, 40])

    def test_read_pfb(self):
        results = file_io_tools.read_file(test_files.forcings_pfb)
        self.assertEqual((24, 41, 41), results.shape)
        self.assertEqual(290.4337549299417, results[0, 0, 0])
        self.assertEqual(292.95937295652664, results[-1, 40, 40])

    def test_read_tif(self):
        results = file_io_tools.read_file(test_files.regression_truth_tif)
        self.assertEqual(3, len(results.shape), 'read a 2d tiff always returns a 3d array')
        results3d = file_io_tools.read_file(test_files.forcings_tif)
        self.assertEqual(3, len(results3d.shape), 'read a 3d tiff always returns a 3d array')

    def test_write_read_bbox(self):
        bbox = [10, 15, 20, 25]
        file_io_tools.write_bbox(bbox, 'bbox_test.txt')
        self.assertSequenceEqual(bbox, file_io_tools.read_bbox('bbox_test.txt'),
                                 'writing and reading a bbox does not change values')
        os.remove('bbox_test.txt')

    def test_write_tiff(self):
        ds_ref = gdal.Open(os.fspath(test_files.regression_truth_tif))
        file_io_tools.write_array_to_geotiff('test_write_tif_out.tif',
                                             file_io_tools.read_file(test_files.regression_truth_tif),
                                             ds_ref.GetGeoTransform(), ds_ref.GetProjection())
        data_array = file_io_tools.read_file('test_write_tif_out.tif')
        self.assertIsNone(np.testing.assert_array_equal(data_array,
                                                        file_io_tools.read_file(
                                                            test_files.regression_truth_tif)),
                          'writing and reading a tif gives back the same array values')
        os.remove('test_write_tif_out.tif')

    def test_write_read_pfb(self):
        forcings_data = file_io_tools.read_file(test_files.forcings_pfb)
        file_io_tools.write_pfb(forcings_data, 'test_pfb_out.pfb')
        read_data = file_io_tools.read_file('test_pfb_out.pfb')
        self.assertIsNone(np.testing.assert_array_equal(forcings_data, read_data),
                          'writing and reading a pfb gives back the same array values')
        os.remove('test_pfb_out.pfb')

    def test_read_pfb_sa_tif(self):
        sa_array = file_io_tools.read_file(test_files.forcings_sa)
        pfb_array = file_io_tools.read_file(test_files.forcings_pfb)
        tif_array = file_io_tools.read_file(test_files.forcings_tif)
        self.assertIsNone(np.testing.assert_array_almost_equal(sa_array, pfb_array, decimal=3),
                          'reading a .sa file and a .pfb file result in same array values')
        self.assertIsNone(np.testing.assert_array_almost_equal(sa_array, tif_array, decimal=3),
                          'reading a .sa file and a .pfb file result in same array values')
        self.assertIsNone(np.testing.assert_array_almost_equal(tif_array, pfb_array, decimal=12),
                          'reading a .sa file and a .pfb file result in same array values')

    def test_read_write_sa_sanity_check(self):
        sa_ref_array = file_io_tools.read_file(test_files.forcings_sa)
        file_io_tools.write_array_to_simple_ascii(data=sa_ref_array, out_file='sa_out_test_file.sa')
        sa_read_array = file_io_tools.read_file('sa_out_test_file.sa')
        self.assertIsNone(np.testing.assert_array_almost_equal(sa_ref_array, sa_read_array, decimal=12),
                          'should be able to write to .sa and get the same data back')
        os.remove('sa_out_test_file.sa')

    def test_write_pfb_to_tif_to_sa(self):
        pfb_array = file_io_tools.read_file(test_files.forcings_pfb)
        srs = osr.SpatialReference()
        srs.SetWellKnownGeogCS("WGS84")
        file_io_tools.write_array_to_geotiff(data=pfb_array, out_raster_path='tif_out_test_forcings_file.tif',
                                             geo_transform=[0, 1000, 0, 0, 0, -1000], projection=srs.ExportToWkt())
        sa_array = file_io_tools.read_file(test_files.forcings_sa)
        tif_array = file_io_tools.read_file('tif_out_test_forcings_file.tif')
        self.assertIsNone(np.testing.assert_array_equal(tif_array, pfb_array,
                                                        'Converting from multi-layer pfb to tif gives back same data'))
        self.assertIsNone(np.testing.assert_array_almost_equal(tif_array, sa_array, decimal=4))
        os.remove('tif_out_test_forcings_file.tif')

    def test_read_write_tif_to_pfb(self):
        tif_array = file_io_tools.read_file(test_files.conus1_dem)
        file_io_tools.write_pfb(tif_array, outfile='conus1_dem_tif_to_pfb_test.pfb')
        pfb_array = file_io_tools.read_file('conus1_dem_tif_to_pfb_test.pfb')
        self.assertIsNone(np.testing.assert_array_equal(tif_array, pfb_array,
                                                        'Converting from tif to pfb gives back same data'))
        os.remove('conus1_dem_tif_to_pfb_test.pfb')

    def test_read_compare_tif_to_pfb(self):
        tif_array = file_io_tools.read_file(test_files.conus1_dem)
        pfb_array = file_io_tools.read_file(test_files.conus1_dem_pfb)
        self.assertIsNone(np.testing.assert_array_equal(tif_array, pfb_array))

    def test_read_write_pfb_to_tif(self):
        tif_array = file_io_tools.read_file(test_files.conus1_dem)
        file_io_tools.write_pfb(tif_array, outfile='conus1_dem_tif_to_pfb_test.pfb')
        pfb_array = file_io_tools.read_file('conus1_dem_tif_to_pfb_test.pfb')
        self.assertIsNone(np.testing.assert_array_equal(tif_array, pfb_array,
                                                        'Converting from single-layer tif to pfb gives back same data'))
        os.remove('conus1_dem_tif_to_pfb_test.pfb')

    def test_read_write_multilayer_tif_pfb(self):
        tif_array = file_io_tools.read_file(test_files.conus2_subsurface)
        file_io_tools.write_pfb(tif_array, outfile='conus2_subsurface_tif_to_pfb_test.pfb')
        pfb_array = file_io_tools.read_file('conus2_subsurface_tif_to_pfb_test.pfb')
        self.assertIsNone(np.testing.assert_array_equal(tif_array, pfb_array,
                                                        'Converting from multi-layer tif to pfb gives back same data'))
        os.remove('conus2_subsurface_tif_to_pfb_test.pfb')

    def test_read_write_tif_sanity_check(self):
        tif_array = file_io_tools.read_file(test_files.conus2_subsurface)
        tif_ref = file_io_tools.read_geotiff(test_files.conus2_subsurface)
        file_io_tools.write_array_to_geotiff('conus2_subsurface_tif_read_write_test.tif',
                                             tif_array, tif_ref.GetGeoTransform(), tif_ref.GetProjectionRef())
        tif_read_array = file_io_tools.read_file('conus2_subsurface_tif_read_write_test.tif')
        self.assertIsNone(np.testing.assert_array_equal(tif_array, tif_read_array),
                          'reading and writing a tif should not change the data')
        os.remove('conus2_subsurface_tif_read_write_test.tif')

    def test_read_invalid_file_type_check(self):
        with self.assertRaises(ValueError):
            file_io_tools.read_file(test_files.huc10190004.get('shapefile'))


if __name__ == '__main__':
    unittest.main()
