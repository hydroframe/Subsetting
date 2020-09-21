import contextlib
import unittest
from parflow.subset.clipper import MaskClipper, BoxClipper
import parflow.subset.utils.io as file_io_tools
from parflow.subset.mask import SubsetMask
import numpy as np
import tests.test_files as test_files
import os


class RegressionClipTests(unittest.TestCase):
    """
    Regression tests1 to verify subsetting can correctly clip a data file,
    correctly produces the subset clip,
    and correctly writes the bounding box file
    """

    def test_subset_dem_to_tif_conus1(self):
        data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())
        my_mask = SubsetMask(test_files.huc10190004.get('conus1_mask').as_posix())
        clipper = MaskClipper(subset_mask=my_mask, no_data_threshold=-1)
        return_arr, new_geom, new_mask, bbox = clipper.subset(data_array)
        file_io_tools.write_array_to_geotiff("conus_1_clip_dem_test.tif",
                                             return_arr, new_geom, my_mask.mask_tif.GetProjection())

        self.assertIsNone(
            np.testing.assert_array_equal(file_io_tools.read_file(test_files.huc10190004.get('conus1_dem').as_posix()),
                                          file_io_tools.read_file('conus_1_clip_dem_test.tif')),
            'Clipping DEM matches reference')
        os.remove('conus_1_clip_dem_test.tif')

        file_io_tools.write_bbox(bbox, 'bbox_conus1.txt')

        self.assertSequenceEqual(file_io_tools.read_bbox('bbox_conus1.txt'), test_files.huc10190004.get('conus1_bbox'),
                                 'Subset writes correct bounding box file')
        os.remove('bbox_conus1.txt')

    def test_subset_tif_conus2(self):
        data_array = file_io_tools.read_file(test_files.conus2_dem.as_posix())
        my_mask = SubsetMask(test_files.huc10190004.get('conus2_mask').as_posix())
        clipper = MaskClipper(subset_mask=my_mask, no_data_threshold=-1)
        return_arr, new_geom, new_mask, bbox = clipper.subset(data_array)
        file_io_tools.write_array_to_geotiff("conus_2_clip_dem_test.tif",
                                             return_arr, new_geom, my_mask.mask_tif.GetProjection())
        self.assertIsNone(
            np.testing.assert_array_equal(file_io_tools.read_file(test_files.huc10190004.get('conus2_dem').as_posix()),
                                          file_io_tools.read_file('conus_2_clip_dem_test.tif')),
            'Clipping DEM matches reference')
        os.remove('conus_2_clip_dem_test.tif')

        file_io_tools.write_bbox(bbox, 'bbox_conus2_full.txt')
        self.assertSequenceEqual(file_io_tools.read_bbox('bbox_conus2_full.txt'),
                                 test_files.huc10190004.get('conus2_bbox'),
                                 'Subset writes correct bounding box file')
        os.remove('bbox_conus2_full.txt')

    def test_compare_box_clips(self):
        data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())
        my_mask = SubsetMask(test_files.huc10190004.get('conus1_mask').as_posix())
        clipper = MaskClipper(subset_mask=my_mask, no_data_threshold=-1)
        mask_subset, _, _, bbox = clipper.subset(data_array, crop_inner=0)

        box_clipper = BoxClipper(ref_array=data_array, x=bbox[0], y=bbox[1], nx=bbox[2], ny=bbox[3])
        box_subset, _, _, _ = box_clipper.subset()
        self.assertEqual(mask_subset.shape[0], box_subset.shape[0])
        self.assertEqual(mask_subset.shape[1], box_subset.shape[1])
        self.assertEqual(mask_subset.shape[2], box_subset.shape[2])
        self.assertIsNone(np.testing.assert_array_equal(mask_subset, box_subset))

    def test_box_clip(self):
        data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())
        box_clipper = BoxClipper(ref_array=data_array)
        subset, _, _, _ = box_clipper.subset()
        self.assertEqual(1, subset.shape[0])
        self.assertEqual(3342, subset.shape[2])
        self.assertEqual(1888, subset.shape[1])
        self.assertIsNone(np.testing.assert_array_equal(data_array, subset),
                          'selecting the whole region should return exactly what you would expect')

        box_clipper.update_bbox(x=10, y=10, nx=3332, ny=1878)
        subset2, _, _, _ = box_clipper.subset()
        self.assertEqual(1, subset2.shape[0])
        self.assertEqual(3332, subset2.shape[2])
        self.assertEqual(1878, subset2.shape[1])

        box_clipper.update_bbox(x=10, y=10, nx=201, ny=20)
        subset3, _, _, _ = box_clipper.subset()
        self.assertEqual(1, subset3.shape[0])
        self.assertEqual(201, subset3.shape[2])
        self.assertEqual(20, subset3.shape[1])

        box_clipper.update_bbox(x=1, y=1, nx=500, ny=300)
        subset4, _, _, _ = box_clipper.subset()
        self.assertEqual(1, subset4.shape[0])
        self.assertEqual(500, subset4.shape[2])
        self.assertEqual(300, subset4.shape[1])

        # create a 3d array for testing, z=4, y=3, x=2
        data_array2 = np.array([[[1, 2], [3, 4, ], [5, 6]],
                                [[7, 8], [9, 10], [11, 12]],
                                [[13, 14], [15, 16], [17, 18]],
                                [[19, 20], [21, 22], [23, 24]]])

        box_clipper2 = BoxClipper(ref_array=data_array2)
        subset5, _, _, _ = box_clipper2.subset()
        self.assertIsNone(np.testing.assert_array_equal(data_array2, subset5))
        self.assertEqual(1, subset5[0, 0, 0])
        self.assertEqual(22, subset5[3, 1, 1])

        box_clipper2.update_bbox(x=1, y=1, nx=1, ny=2)
        subset6, _, _, _ = box_clipper2.subset()
        self.assertEqual(1, subset6[0, 0, 0])
        self.assertEqual(13, subset6[2, 0, 0])
        self.assertEqual(15, subset6[2, 1, 0])

    def test_box_clip_with_padding(self):
        data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())
        # css-like padding (top,right,bot,left)
        bbox = test_files.huc10190004.get('conus1_bbox')
        box_clipper = BoxClipper(ref_array=data_array, x=bbox[0], y=bbox[1], nx=bbox[2], ny=bbox[3],
                                 padding=(1, 6, 1, 5))
        subset, _, _, _ = box_clipper.subset()
        self.assertEqual(1, subset.shape[0])
        self.assertEqual(32, subset.shape[1])
        self.assertEqual(96, subset.shape[2])
        file_io_tools.write_pfb(subset, 'WBDHU8_conus1_dem_padded_test.pfb')
        padded_subset_ref = file_io_tools.read_file(test_files.huc10190004.get('conus1_dem_padded_box').as_posix())
        self.assertIsNone(np.testing.assert_array_equal(padded_subset_ref, subset))
        os.remove('WBDHU8_conus1_dem_padded_test.pfb')

    def test_box_clip_invalid_nx_dim(self):
        with self.assertRaises(Exception):
            data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())
            BoxClipper(ref_array=data_array, nx=0)

    def test_box_clip_invalid_x_dim(self):
        with self.assertRaises(Exception):
            data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())
            BoxClipper(ref_array=data_array, x=0)

    def test_box_print_no_exception(self):
        data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())
        clipper = BoxClipper(ref_array=data_array)
        self.assertEqual(-999, clipper.no_data)
        self.assertIsNone(print(clipper))

    def test_mask_print_no_exception(self):
        my_mask = SubsetMask(test_files.huc10190004.get('conus1_mask').as_posix())
        clipper = MaskClipper(subset_mask=my_mask, no_data_threshold=-1)
        self.assertIsNone(print(clipper))


if __name__ == '__main__':
    unittest.main()
