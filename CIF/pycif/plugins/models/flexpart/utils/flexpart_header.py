""" Class for reading FLEXPART binary headers

"""

import datetime
import numpy as np


try:
    from . import mod_flexpart
except ImportError:
    raise Exception("Error importing module mod_flexpart, make sure it is compiled.\n"
                    "Normally this should be done by the setup script.\n"
                    "Try going to pycif/plugins/models/flexpart/utils and execute\n"
                    "'f2py -c mod_flexpart.f90 -m mod_flexpart'\n")


class Flexpartheader(object):
    """ Class for reading FLEXPART binary headers

    """

    def __init__(self):
        """Initialize a Flexpartheader object

        """

        self.numpoint = np.nan
        self.ibdate = np.nan
        self.ibtime = np.nan
        self.releases = np.nan
        self.numx = np.nan
        self.numy = np.nan
        self.bndx = np.nan
        self.bndy = np.nan
        self.delx = np.nan
        self.dely = np.nan
        self.xshift = np.nan
        self.ageclass = np.nan
        self.trajdays = np.nan
        self.ndgrid = np.nan
        self.outheight = np.nan
        self.maxngrid = np.nan


    def read_header(self, filename):
        """ Reads header using compiled fortran module mod_flexpart

        """

        self.numpoint, self.ibdate, self.ibtime, self.releases, \
        self.numx, self.numy, self.bndx, self.bndy, self.delx, self.dely, \
        self.xshift, self.ageclass, self.trajdays, self.ndgrid, self.outheight = \
        mod_flexpart.read_header(filename)

        self.maxngrid = self.trajdays*self.ndgrid + 2 


