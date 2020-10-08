import unittest
from parflow.subset.mask import SubsetMask
from tests.test_files import huc10190004, test_all_zeros_and_ones_mask, test_padded_icom_mask


class SubsetMaskClassTests(unittest.TestCase):

    def test_normal_startup(self):
        my_mask = SubsetMask(huc10190004.get('conus1_mask'))
        self.assertEqual(-999, my_mask.no_data_value)
        self.assertEqual(0, my_mask.bbox_val)
        self.assertSequenceEqual((0, 0, 0, 0), my_mask.get_bbox().get_padding())

    def test_resize_bbox(self):
        my_mask = SubsetMask(huc10190004.get('conus1_mask'))
        self.assertSequenceEqual((716, 745, 1039, 1123), my_mask.inner_mask_edges)
        self.assertSequenceEqual((716, 745, 1039, 1123), my_mask.bbox_edges)
        self.assertSequenceEqual((30, 85), my_mask.bbox_shape)
        self.assertSequenceEqual((30, 85), my_mask.inner_mask_shape)
        my_mask.add_bbox_to_mask(padding=(9, 9, 9, 9))
        self.assertSequenceEqual((707, 754, 1030, 1132), my_mask.bbox_edges)
        self.assertSequenceEqual((716, 745, 1039, 1123), my_mask.inner_mask_edges)
        self.assertSequenceEqual((30, 85), my_mask.inner_mask_shape)
        self.assertSequenceEqual((48, 103), my_mask.bbox_shape)
        self.assertEqual(-999, my_mask.no_data_value)
        self.assertEqual(0, my_mask.bbox_val)
        self.assertSequenceEqual((9, 9, 9, 9), my_mask.get_bbox().get_padding())
        my_mask.write_mask_to_tif("WBDHU8_conus1_mask_padded.tif")

    def test_all_zeroes_ones_mask(self):
        my_mask = SubsetMask(test_all_zeros_and_ones_mask)
        self.assertSequenceEqual((1752, 2222, 3672, 4078), my_mask.inner_mask_edges)
        self.assertSequenceEqual((1752, 2222, 3672, 4078), my_mask.bbox_edges)
        self.assertSequenceEqual((471, 407), my_mask.bbox_shape)
        self.assertSequenceEqual((471, 407), my_mask.inner_mask_shape)
        my_mask.add_bbox_to_mask(padding=(4, 5, 5, 4))
        self.assertSequenceEqual((1752, 2222, 3672, 4078), my_mask.inner_mask_edges)
        self.assertSequenceEqual((1747, 2226, 3668, 4083), my_mask.bbox_edges)
        self.assertSequenceEqual((480, 416), my_mask.bbox_shape)
        self.assertSequenceEqual((471, 407), my_mask.inner_mask_shape)

    def test_get_padding_from_mask(self):
        my_mask = SubsetMask(test_all_zeros_and_ones_mask)
        self.assertSequenceEqual((1752, 2222, 3672, 4078), my_mask.inner_mask_edges)
        self.assertSequenceEqual((1752, 2222, 3672, 4078), my_mask.bbox_edges)
        self.assertSequenceEqual((471, 407), my_mask.bbox_shape)
        self.assertSequenceEqual((471, 407), my_mask.inner_mask_shape)
        self.assertSequenceEqual((0, 0, 0, 0), my_mask.get_padding())
        my_mask.add_bbox_to_mask(padding=(4, 5, 5, 4))
        self.assertSequenceEqual((1752, 2222, 3672, 4078), my_mask.inner_mask_edges)
        self.assertSequenceEqual((1747, 2226, 3668, 4083), my_mask.bbox_edges)
        self.assertSequenceEqual((480, 416), my_mask.bbox_shape)
        self.assertSequenceEqual((471, 407), my_mask.inner_mask_shape)
        self.assertSequenceEqual((4, 5, 5, 4), my_mask.get_padding())
        self.assertSequenceEqual((4, 5, 5, 4), my_mask.get_bbox().get_padding())

    def test_open_padded_mask(self):
        my_mask = SubsetMask(test_padded_icom_mask)
        self.assertSequenceEqual((4, 5, 5, 4), my_mask.get_padding())
        self.assertSequenceEqual((1752, 2222, 3672, 4078), my_mask.inner_mask_edges)
        self.assertSequenceEqual((1747, 2226, 3668, 4083), my_mask.bbox_edges)
        self.assertSequenceEqual((480, 416), my_mask.bbox_shape)
        self.assertSequenceEqual((471, 407), my_mask.inner_mask_shape)


if __name__ == '__main__':
    unittest.main()
