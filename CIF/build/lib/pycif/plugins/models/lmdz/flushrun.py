from glob import glob

from pycif.utils import path


def flushrun(self, rundir, mode):
    """Cleaning the simulation directories to limit space usage"""

    list_subdirs = glob("{}/*/".format(rundir))

    # Removing all binary files used by LMDZ as interface with pycif
    for subdir in list_subdirs:
        if "chain" not in subdir:
            path.remove("{}/*.bin".format(subdir))

        else:
            path.remove("{}/restart*.nc".format(subdir))

    # Removing trajq binary files from previous forward simulations if any
    if mode == "adj" and getattr(self, "adj_refdir", False):
        path.remove("{}/chain/traj*.bin".format(self.adj_refdir))
