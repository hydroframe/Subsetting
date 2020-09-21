import os
import numpy as np
from unittest import TestCase
from parflow.subset.tools import rasterize_shape
from parflow.subset.utils.io import read_file, read_bbox
import test_files


class RaserizeShapeCLITests(TestCase):

    def setUp(self) -> None:
        self.good_shape_file_path = test_files.huc10190004.get('shapefile').parent
        self.good_shape_file_name = test_files.huc10190004.get('shapefile').stem

    def test_parse_args_empty(self):
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args([])

    def test_parse_path_and_ref_args(self):
        argstring = f'-i {self.good_shape_file_path} -r {test_files.conus1_mask}'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_path_and_shape_args(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name}'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_shape_and_ref_args(self):
        argstring = f'-s {self.good_shape_file_name} -r {test_files.conus1_mask}'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_bad_input_path(self):
        argstring = f'-i path_no_exists -s {self.good_shape_file_name}-r {test_files.conus1_mask}'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_bad_ref_file(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name}-r tif_no_exists.tif'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_bad_input_val_out_dir(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name}-r {test_files.conus1_mask} ' \
                    f'-o path_no_exists'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_bad_padding_char(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name}-r {test_files.conus1_mask} ' \
                    f'-p A 0 1 1'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_bad_padding_neg(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name}-r {test_files.conus1_mask} ' \
                    f'-p -1 0 1 1'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_bad_padding_too_many_vals(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name}-r {test_files.conus1_mask} ' \
                    f'-p 1 0 1 1 5'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_bad_padding_too_few_vals(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name}-r {test_files.conus1_mask} ' \
                    f'-p 1 1 5'
        with self.assertRaises(SystemExit):
            rasterize_shape.parse_args(argstring.split(' '))

    def test_parse_args_min(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name} -r {test_files.conus1_mask}'
        args = rasterize_shape.parse_args(argstring.split(' '))
        self.assertEqual(os.fspath(self.good_shape_file_path), args.input_path)
        self.assertEqual(self.good_shape_file_name, args.shapefile)
        self.assertEqual(os.fspath(test_files.conus1_mask), args.ref_file)
        self.assertEqual('.', args.out_dir)
        self.assertIsNone(args.out_file)
        self.assertSequenceEqual((0, 0, 0, 0), args.padding)
        self.assertSequenceEqual([1], args.attribute_ids)
        self.assertEqual('OBJECTID', args.attribute_name)

    def test_parse_args_all_non_default(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name} -r {test_files.conus1_mask} ' \
                    f'-o .. -n test_subset_name -p 1 2 3 4 -a 5 6 -e ID'
        args = rasterize_shape.parse_args(argstring.split(' '))
        self.assertEqual(os.fspath(self.good_shape_file_path), args.input_path)
        self.assertEqual(self.good_shape_file_name, args.shapefile)
        self.assertEqual(os.fspath(test_files.conus1_mask), args.ref_file)
        self.assertEqual('..', args.out_dir)
        self.assertEqual('test_subset_name', args.out_file)
        self.assertSequenceEqual((1, 2, 3, 4), args.padding)
        self.assertSequenceEqual([5, 6], args.attribute_ids)
        self.assertEqual('ID', args.attribute_name)


class RaserizeShapeRegressionTests(TestCase):

    def setUp(self) -> None:
        self.good_shape_file_path = test_files.huc10190004.get('shapefile').parent
        self.good_shape_file_name = test_files.huc10190004.get('shapefile').stem

    def test_conus1_raster_defaults(self):
        rasterize_shape.rasterize_shape(input_path=self.good_shape_file_path, shapefile=self.good_shape_file_name,
                                        ref_file=test_files.conus1_mask)
        ref_mask = read_file(test_files.huc10190004.get('conus1_mask'))
        ref_bbox = test_files.huc10190004.get('conus1_bbox')
        self.assertSequenceEqual(ref_bbox, read_bbox('bbox.txt'))
        self.assertIsNone(np.testing.assert_array_equal(ref_mask, read_file('WBDHU8.tif')))
        self.assertIsNone(os.remove('bbox.txt'))
        self.assertIsNone(os.remove('WBDHU8.tif'))
