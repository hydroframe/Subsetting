import os
import shutil
import unittest
from pathlib import Path
import numpy as np
from tests import test_files
from parflow.subset.tools import subset_conus
from parflow.subset.utils.io import read_file


class SubsetConusCLITests(unittest.TestCase):

    def setUp(self) -> None:
        self.good_mask_file = test_files.huc10190004.get('conus1_mask')
        self.bad_mask_file = './mask_file_no_exists.tif'
        self.good_input_file_list = [test_files.conus1_dem, test_files.conus1_mask]
        self.bad_input_file_list = './input_file_to_clip_no_exists.pfb'
        self.good_bbox_file = test_files.test_bbox_input
        self.good_shape_file_path = test_files.huc10190004.get('shapefile').parent
        self.good_shape_file_name = test_files.huc10190004.get('shapefile').stem
        self.partial_shapefile_path = test_files.partial_shape_file.parent
        self.partial_shapefile_name = test_files.partial_shape_file.stem

    def test_good_start(self):
        args = subset_conus.parse_args(
            ['-i', os.fspath(self.good_shape_file_path), '-s', self.good_shape_file_name, '-f', '.'])
        self.assertEqual(os.fspath(self.good_shape_file_path), args.input_path)
        self.assertEqual(self.good_shape_file_name, args.shapefile)
        self.assertEqual('.', args.conus_files)
        self.assertEqual(1, args.conus_version)
        self.assertEqual('.', args.out_dir)
        self.assertIsNone(args.out_name)
        self.assertFalse(args.clip_clm)
        self.assertFalse(args.write_tcl)
        self.assertSequenceEqual((0, 0, 0, 0), args.padding)
        self.assertSequenceEqual([1], args.attribute_ids)
        self.assertEqual('OBJECTID', args.attribute_name)
        self.assertFalse(args.write_tifs)

    def test_cli_no_args(self):
        """should error without arguments"""
        with self.assertRaises(SystemExit):
            subset_conus.parse_args([])

    def test_cli_min_args_bad_path(self):
        """Should error with bad shapfile path"""
        with self.assertRaises(SystemExit):
            subset_conus.parse_args(['-i', 'path_no_exists', '-s', self.good_shape_file_name, '-f', '.'])

    def test_cli_min_args_bad_shape_name(self):
        """Should not error with bad shapefile name, not because we don't want to.
        Separating the input path from the file name means we can't check for the file while parsing args
        """
        args = subset_conus.parse_args(
            ['-i', os.fspath(self.good_shape_file_path), '-s', 'shape_name_no_exists', '-f', '.'])
        self.assertEqual('shape_name_no_exists', args.shapefile)
        self.assertEqual(os.fspath(self.good_shape_file_path), args.input_path)

    def test_cli_min_args_missing_shape_parts(self):
        """A shapefile with missing components will not raise an exception"""
        args = subset_conus.parse_args(
            ['-i', os.fspath(self.partial_shapefile_path), '-s', self.partial_shapefile_name, '-f', '.'])
        self.assertEqual(os.fspath(self.partial_shapefile_path), args.input_path)
        self.assertEqual(self.partial_shapefile_name, args.shapefile)

    def test_cli_all_options_non_default(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name} ' \
                    f'-f . -v 2 -o .. -n output_name -c -w -p 1 2 3 4 -a 2 3 -e ID -t'.split(' ')
        args = subset_conus.parse_args(argstring)
        self.assertTrue(args.write_tifs)
        self.assertTrue(args.write_tcl)
        self.assertTrue(args.clip_clm)
        self.assertEqual(self.good_shape_file_name, args.shapefile)
        self.assertEqual(os.fspath(self.good_shape_file_path), args.input_path)
        self.assertEqual(2, args.conus_version)
        self.assertEqual('.', args.conus_files)
        self.assertEqual('..', args.out_dir)
        self.assertEqual('output_name', args.out_name)
        self.assertSequenceEqual((1, 2, 3, 4), args.padding)
        self.assertSequenceEqual([2, 3], args.attribute_ids)
        self.assertEqual('ID', args.attribute_name)

    def test_cli_alt_manifest_file(self):
        argstring = f'-i {self.good_shape_file_path} -s {self.good_shape_file_name} -f . ' \
                    f'-m {test_files.test_domain_manifest}'
        args = subset_conus.parse_args(argstring.split(' '))
        self.assertEqual(os.fspath(test_files.test_domain_manifest), args.manifest_file,
                         'should be able to provide an alternate manifest file')


class SubsetConusRegressionTests(unittest.TestCase):

    def setUp(self) -> None:
        self.good_mask_file = test_files.huc10190004.get('conus1_mask')
        self.bad_mask_file = './mask_file_no_exists.tif'
        self.good_input_file_list = [test_files.conus1_dem, test_files.conus1_mask]
        self.bad_input_file_list = './input_file_to_clip_no_exists.pfb'
        self.good_bbox_file = test_files.test_bbox_input
        self.good_shape_file_path = test_files.huc10190004.get('shapefile').parent
        self.good_shape_file_name = test_files.huc10190004.get('shapefile').stem
        self.partial_shapefile_path = test_files.partial_shape_file.parent
        self.partial_shapefile_name = test_files.partial_shape_file.stem

    def test_conus1_subset_regression(self):
        if os.environ.get('TRAVIS'):
            pass
        else:
            test_dir = Path('test_outputs')
            test_dir.mkdir(exist_ok=True)
            subset_conus.subset_conus(input_path=self.good_shape_file_path,
                                      shapefile=self.good_shape_file_name,
                                      conus_files='/home/arezaii/git/subset_1/CONUS1_inputs',
                                      out_dir=test_dir,
                                      out_name='test_conus1_subset')
            ref_subsurface = read_file(test_files.huc10190004.get('conus1_subsurface'))
            ref_mask = read_file(test_files.huc10190004.get('conus1_mask'))
            ref_slopex = read_file(test_files.huc10190004.get('conus1_slopex'))
            ref_slopey = read_file(test_files.huc10190004.get('conus1_slopey'))

            self.assertIsNone(np.testing.assert_array_equal(ref_subsurface, read_file(test_dir / 'grid3d.v3_clip.pfb')))
            self.assertIsNone(np.testing.assert_array_equal(ref_mask, read_file(
                test_dir / 'test_conus1_subset_raster_from_shapefile.tif')))
            self.assertIsNone(np.testing.assert_array_equal(ref_slopex, read_file(test_dir / 'slopex_clip.pfb')))
            self.assertIsNone(np.testing.assert_array_equal(ref_slopey, read_file(test_dir / 'slopey_clip.pfb')))

            with open(test_dir / 'test_conus1_subset.pfsol', 'r') as test_file:
                with open(test_files.huc10190004.get('conus1_sol'), 'r') as ref_file:
                    self.assertEqual(test_file.read(), ref_file.read(),
                                     'Writing PFSOL file matches reference for conus1')
            """
            make sure the vtk files match as well. skip the first two lines as they contain a version and filename that
            may vary
            """
            with open(test_dir / 'test_conus1_subset.vtk', 'r') as test_file:
                with open(test_files.huc10190004.get('conus1_vtk'), 'r') as ref_file:
                    self.assertEqual(test_file.read().split('\n')[2:], ref_file.read().split('\n')[2:],
                                     'Writing vtk file matches reference for conus1')
            shutil.rmtree(test_dir)


if __name__ == '__main__':
    unittest.main()
