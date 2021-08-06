import os
from glob import glob


def flushrun(self, rundir, mode):
    """Cleaning the simulation directories to limit space usage"""

    list_subdirs = glob("{}/*/".format(rundir))

    # Removing all binary files used by LMDZ as interface with pycif
    for subdir in list_subdirs:
        if "chain" not in subdir:
            toremove = glob("{}/*.bin".format(subdir))

        else:
            toremove = glob("{}/restart*.nc".format(subdir))

        for rm_file in toremove:
            os.remove(rm_file)

    # Removing trajq binary files from previous forward simulations if any
    if mode == "adj" and getattr(self, "adj_refdir", False):
        for rm_file in glob("{}/traj*.bin".format(self.adj_refdir)):
            os.remove(rm_file)
