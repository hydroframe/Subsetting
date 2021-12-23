import os
import unittest
import numpy as np
import pfsubset.subset.utils.io as file_io_tools
import tests.test_files as test_files
from pfsubset.subset.clipper import ClmClipper
from pfsubset.subset.bbox import BBox
from pfsubset.subset.mask import SubsetMask


class ClmUtilsClipperRegressionTests(unittest.TestCase):

    # The input file for LAT/LON values is over 100MB, so can't be added to github...

    def test_clm_clip_land_cover(self):
        if os.environ.get('TRAVIS'):
            pass
        elif os.path.isfile(test_files.conus1_latlon):
            bbox_list = test_files.huc10190004.get('conus1_bbox')
            bbox = BBox(bbox_list[0], bbox_list[1], bbox_list[2], bbox_list[3])
            clm_clipper = ClmClipper(bbox)

            latlon_data, _ = clm_clipper.clip_latlon(test_files.conus1_latlon)
            land_cover_data, vegm_data = clm_clipper.clip_land_cover(lat_lon_array=latlon_data,
                                                                     land_cover_file=test_files.conus1_landcover)
            clm_clipper.write_land_cover(vegm_data, 'WBDHU8_vegm_test.dat')
            with open('WBDHU8_vegm_test.dat', 'r') as test_file:
                with open(test_files.huc10190004.get('conus1_vegm').as_posix(), 'r') as ref_file:
                    self.assertEqual(test_file.read().split('\n'), ref_file.read().split('\n'),
                                     'Writing vegm file matches reference for conus1')

            os.remove('WBDHU8_vegm_test.dat')
        else:
            print('WARNING! Unable to run test test_clm_clip_latlon because source file not found. '
                  'copy conus1_Grid_Centers_Short_Deg.format.sa into test_inputs/CONUS1_Inputs to enable test!')
            pass

    def test_clm_clip_latlon(self):
        if os.environ.get('TRAVIS'):
            pass
        elif os.path.isfile(test_files.conus1_latlon):
            bbox_list = test_files.huc10190004.get('conus1_bbox')
            bbox = BBox(bbox_list[0], bbox_list[1], bbox_list[2], bbox_list[3])
            clm_clipper = ClmClipper(bbox)
            latlon_formatted, latlon_data = clm_clipper.clip_latlon(test_files.conus1_latlon)
            self.assertSequenceEqual((1, 30, 85), latlon_data.shape, 'check that clm clip when mask has no padding works')
            clm_clipper.write_lat_lon(latlon_formatted, 'WBDHU8_latlon_test.sa', x=latlon_data.shape[2],
                                      y=latlon_data.shape[1], z=latlon_data.shape[0])
            self.assertIsNone(np.testing.assert_array_equal(file_io_tools.read_file('WBDHU8_latlon_test.sa'),
                                                            file_io_tools.read_file(
                                                                test_files.huc10190004.get('conus1_latlon'))),
                              'clipping conus1 lat/lon values regression test works')
            os.remove('WBDHU8_latlon_test.sa')
        else:
            print('WARNING! Unable to run test test_clm_clip_latlon because source file not found. '
                  'copy conus1_Grid_Centers_Short_Deg.format.sa into test_inputs/CONUS1_Inputs to enable test!')
            pass

    def test_clm_clip_latlon_padded(self):
        if os.environ.get('TRAVIS'):
            pass
        elif os.path.isfile(test_files.conus1_latlon):
            my_mask = SubsetMask(test_files.huc10190004.get('conus1_padded_mask'))
            bbox = my_mask.get_bbox()
            self.assertSequenceEqual((9, 9, 9, 9), bbox.get_padding())
            clm_clipper = ClmClipper(bbox)
            latlon_formatted, latlon_data = clm_clipper.clip_latlon(test_files.conus1_latlon)
            clm_clipper.write_lat_lon(latlon_formatted, 'WBDHU8_latlon_test.sa', x=latlon_data.shape[2],
                                      y=latlon_data.shape[1], z=latlon_data.shape[0])
            self.assertSequenceEqual((1, 48, 103), latlon_data.shape, 'check that clm clip when mask has padding works')
            os.remove('WBDHU8_latlon_test.sa')
        else:
            print('WARNING! Unable to run test test_clm_clip_latlon because source file not found. '
                  'copy conus1_Grid_Centers_Short_Deg.format.sa into test_inputs/CONUS1_Inputs to enable test!')
            pass

    def test_clm_clip_land_cover_padded(self):
        if os.environ.get('TRAVIS'):
            pass
        elif os.path.isfile(test_files.conus1_latlon):
            my_mask = SubsetMask(test_files.huc10190004.get('conus1_padded_mask'))
            bbox = my_mask.get_bbox()
            self.assertSequenceEqual((9, 9, 9, 9), bbox.get_padding())
            clm_clipper = ClmClipper(bbox)

            latlon_data, _ = clm_clipper.clip_latlon(test_files.conus1_latlon)
            land_cover_data, vegm_data = clm_clipper.clip_land_cover(lat_lon_array=latlon_data,
                                                                     land_cover_file=test_files.conus1_landcover)
            clm_clipper.write_land_cover(vegm_data, 'WBDHU8_vegm_test.dat')


            os.remove('WBDHU8_vegm_test.dat')
        else:
            print('WARNING! Unable to run test test_clm_clip_latlon because source file not found. '
                  'copy conus1_Grid_Centers_Short_Deg.format.sa into test_inputs/CONUS1_Inputs to enable test!')
            pass


if __name__ == '__main__':
    unittest.main()
