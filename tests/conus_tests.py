import unittest
import pfsubset.subset.data as data
from pfsubset.subset.domain import Conus
import os
import tests.test_files as test_files

class ConusClassTests(unittest.TestCase):

    # Can't perform positive test case without all the required_files present
    def test_normal_startup(self):
        if os.getenv('TRAVIS', default=False):
            pass
        elif os.path.isdir(test_files.conus1_inputs): 
            conus1 = Conus(local_path=test_files.conus1_inputs)
            self.assertEqual(conus1.manifest_file, data.conus_manifest)
        else:
            print('CONUS1_inputs not found, skipping test.')
            pass


if __name__ == '__main__':
    unittest.main()
