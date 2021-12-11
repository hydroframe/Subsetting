import contextlib
import unittest
from parflow.subset.clipper import MaskClipper, BoxClipper
import parflow.subset.utils.io as file_io_tools
from parflow.subset.mask import SubsetMask
import numpy as np
import tests.test_files as test_files
import os




data_file = test_files.conus1_dem.as_posix()
my_mask = SubsetMask(test_files.huc10190004.get('conus1_mask').as_posix())

mask_file = test_files.huc10190004.get('conus1_mask').as_posix()
clipper = MaskClipper(mask_file=mask_file, no_data_threshold=-1)
mask_subset, _, _, bbox = clipper.subset(data_file, crop_inner=0)

data_array = file_io_tools.read_file(test_files.conus1_dem.as_posix())

box_clipper = BoxClipper(nx=10,ny=5,nz=5)
print("--------------")
print(data_file)
box_subset, _, _, _ = box_clipper.subset(data_file)
print("------thesubset---------")
print(box_subset)
print("abdf dfadf ")
