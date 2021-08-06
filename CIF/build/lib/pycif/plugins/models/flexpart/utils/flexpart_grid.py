""" Class for reading FLEXPART binary output

"""

import numpy as np
try:
    import read_flexpart
except ImportError:
    raise Exception("Error importing module read_flexpart, make sure it is compiled")


class Flexpartgrid(object):
    """ Class for reading FLEXPART binary headers

    """

    def __init__(self):
        """Initialize a Flexpartheader object

        """

        self.grid = np.nan
        self.filedates = np.nan
        self.jd = np.nan
        self.numx = np.nan
        self.numy = np.nan
        self.maxngrid = np.nan
        self.ageclass = np.nan


    def read_grid(self, filename):
        """ Reads FLEXPART grid files using compiled fortran module read_flexpart

        """

        self.grid, self.filedates, self.jd, self.numx, self.numy, \
            self.xshift, self.ngrid, self. gtime, \
            self.ndgrid, self.trajdays = read_flexpart.read_grid(filename)
