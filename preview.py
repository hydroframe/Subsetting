#!/usr/bin/env python3


import sys
import pfio
import numpy
import matplotlib.pyplot as plt


def preview(f):

    arr = pfio.pfread(f)

#    import pdb; pdb.set_trace()
    # set 0 to null
    arr[arr == 0.0] = numpy.nan

    # set non-zero to 1
    arr[~numpy.isnan(arr)] = 1

    imgplot = plt.imshow(arr[0, :, :])
    plt.show()


if __name__ == '__main__':

    args = sys.argv[1:]
    f = args[0]
    preview(f)
