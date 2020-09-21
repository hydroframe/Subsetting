import os
import unittest
import parflow.subset.builders.solidfile as solidfile_generator
import tests.test_files as test_files
from parflow.subset.clipper import MaskClipper
from parflow.subset.mask import SubsetMask


class RegressionSolidFileTests(unittest.TestCase):
    """
    Regression tests1 to verify solidfile_generator creates solid file correctly from mask
    """

    asc_filenames = ['Back_Border.asc', 'Bottom_Border.asc', 'Front_Border.asc',
                     'Left_Border.asc', 'Right_Border.asc', 'Top_Border.asc']

    def tearDown(self):
        for f in self.asc_filenames:
            try:
                os.remove(f)
            except OSError:
                pass

    def test_create_solid_file_conus1(self):
        my_mask = SubsetMask(test_files.huc10190004.get('conus1_mask').as_posix())
        clipper = MaskClipper(subset_mask=my_mask, no_data_threshold=-1)
        batches = solidfile_generator.make_solid_file(clipped_mask=clipper.clipped_mask, out_name='conus1_solid')
        self.assertEqual(batches, '0 3 6 ')
        with open('conus1_solid.pfsol', 'r') as test_file:
            with open(test_files.huc10190004.get('conus1_sol'), 'r') as ref_file:
                self.assertEqual(test_file.read(), ref_file.read(),
                                 'Writing PFSOL file matches reference for conus1')
        """
        make sure the vtk files match as well. skip the first two lines as they contain a version and filename that
        may vary
        """
        with open('conus1_solid.vtk', 'r') as test_file:
            with open(test_files.huc10190004.get('conus1_vtk'), 'r') as ref_file:
                self.assertEqual(test_file.read().split('\n')[2:], ref_file.read().split('\n')[2:],
                                 'Writing vtk file matches reference for conus1')
        os.remove('conus1_solid.vtk')
        os.remove('conus1_solid.pfsol')

    def test_create_solid_file_conus2(self):
        my_mask = SubsetMask(test_files.huc10190004.get('conus2_mask').as_posix())
        clipper = MaskClipper(subset_mask=my_mask, no_data_threshold=-1)
        batches = solidfile_generator.make_solid_file(clipper.clipped_mask, 'conus2_solid')
        self.assertEqual(batches, '0 3 6 ')
        with open('conus2_solid.pfsol', 'r') as test_file:
            with open(test_files.huc10190004.get('conus2_sol'), 'r') as ref_file:
                self.assertEqual(test_file.read(), ref_file.read(),
                                 'Writing PFSOL file matches reference for conus2')
        """
        make sure the vtk files match as well. skip the first two lines as they contain a version and filename that
        may vary
        """
        with open('conus2_solid.vtk', 'r') as test_file:
            with open(test_files.huc10190004.get('conus2_vtk').as_posix(), 'r') as ref_file:
                self.assertEqual(test_file.read().split('\n')[2:], ref_file.read().split('\n')[2:],
                                 'Writing vtk file matches reference for conus2')
        for f in self.asc_filenames:
            os.remove(f)
        os.remove('conus2_solid.vtk')
        os.remove('conus2_solid.pfsol')


if __name__ == '__main__':
    unittest.main()
