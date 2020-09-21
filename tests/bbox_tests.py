import unittest
from parflow.subset.bbox import BBox


class BBoxTestCases(unittest.TestCase):
    def test_bbox(self):
        bbox = BBox(1125, 717, 85, 30)
        self.assertSequenceEqual(bbox.get_inner_extents(), (716, 746, 1124, 1209))
        self.assertSequenceEqual(bbox.get_human_bbox(), (1125, 717, 85, 30))
        self.assertSequenceEqual(bbox.get_system_bbox(), (1124, 716, 85, 30))
        self.assertSequenceEqual(bbox.get_padded_extents(), (716, 746, 1124, 1209))

    def test_bbox_init(self):
        bbox = BBox(1, 1)
        self.assertSequenceEqual(bbox.get_inner_extents(), (0, 0, 0, 0))
        self.assertSequenceEqual(bbox.get_human_bbox(), (1, 1, 0, 0))
        self.assertSequenceEqual(bbox.get_system_bbox(), (0, 0, 0, 0))
        self.assertSequenceEqual(bbox.get_padded_extents(), (0, 0, 0, 0))

    def test_bbox_single_cell(self):
        bbox = BBox(1, 1, 1, 1)
        self.assertSequenceEqual(bbox.get_inner_extents(), (0, 1, 0, 1))
        self.assertSequenceEqual(bbox.get_human_bbox(), (1, 1, 1, 1))
        self.assertSequenceEqual(bbox.get_system_bbox(), (0, 0, 1, 1))
        self.assertSequenceEqual(bbox.get_padded_extents(), (0, 1, 0, 1))

    def test_bbox_padding(self):
        # bbox padding sequence matches CSS, (top, right, bottom, left) = (y_upper, x_upper, y_lower, x_lower)
        bbox = BBox(4, 4, 20, 10, (2, 2, 1, 1))
        self.assertSequenceEqual(bbox.get_inner_extents(), (3, 13, 3, 23))
        self.assertSequenceEqual(bbox.get_human_bbox(), (4, 4, 20, 10))
        self.assertSequenceEqual(bbox.get_system_bbox(), (3, 3, 20, 10))
        self.assertSequenceEqual(bbox.get_padded_extents(), (2, 15, 2, 25))


if __name__ == '__main__':
    unittest.main()
