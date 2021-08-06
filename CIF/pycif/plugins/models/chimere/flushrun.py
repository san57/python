from glob import glob

from pycif.utils import path


def flushrun(self, rundir, mode):
    """Cleaning the simulation directories to limit space usage"""

    list_subdirs = glob("{}/*/".format(rundir))

    # Removing big files in CHIMERE directories
    for subdir in list_subdirs:
        path.remove("{}/AEMISSIONS.nc".format(subdir))
        path.remove("{}/BOUN_CONCS.nc".format(subdir))
        path.remove("{}/INI_CONCS.0.nc".format(subdir))
        path.remove("{}/METEO.nc".format(subdir))
        path.remove("{}/par.nc".format(subdir))
        path.remove("{}/dep.nc".format(subdir))
        path.remove("{}/end.nc".format(subdir))
