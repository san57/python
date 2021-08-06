import calendar
import datetime
import glob
import os
import shutil
import subprocess

import numpy as np

from pycif.utils import path
from pycif.utils.check import info
from pycif.utils.netcdf import save_nc


def run(self, runsubdir, mode, workdir, do_simu=True, **kwargs):
    """Run LMDZ model in forward or adjoint mode

    Args:
        runsubdir (str): working directory for the current run
        mode (str): forward or backward
        workdir (str): pycif working directory
        do_simu (bool): if False, considers that the simulation was
                        already run

    """

    if do_simu:
        # Cleaning the directory
        path.remove("{}/all_good".format(runsubdir))
        path.remove("{}/restart.nc".format(runsubdir))
        path.remove("{}/obs_out.bin".format(runsubdir))
        path.remove("{}/*_out.bin".format(runsubdir))

        if mode != "adj":
            path.remove("{}/traj*.bin".format(runsubdir))

        # Running the model
        info("Running sub-simulation in {}".format(runsubdir))
        with open("{}/dispersion.log".format(runsubdir), "w") as log:
            process = subprocess.Popen(
                ["mpirun", "{}/dispersion.e".format(runsubdir)],
                cwd=runsubdir,
                stdout=log,
                stderr=subprocess.PIPE,
            )
            _, err = process.communicate()

            if err != "" and not os.path.isfile(
                "{}/all_good".format(runsubdir)
            ):
                info("Exception in LMDZ")
                info(err)
                raise Exception(
                    "LMDZ did not run properly in {}".format(runsubdir)
                )

        # Dumps to NetCDF if asked
        if getattr(self, "dump", False):
            info(
                "Dumping LMDZ outputs to a NetCDF file in {}".format(runsubdir)
            )
            dump2nc(self, runsubdir)

        # Process here initial conditions for the next simulation
        date = datetime.datetime.strptime(
            os.path.basename(runsubdir), "%Y-%m-%d_%H-%M"
        )
        shutil.move(
            "{}/restart.nc".format(runsubdir),
            date.strftime(
                "{}/../chain/restart_%Y%m%d%H%M.nc".format(runsubdir)
            ),
        )

        if mode == "tl":
            shutil.move(
                "{}/restart_tl.bin".format(runsubdir),
                date.strftime(
                    "{}/../chain/restart_tl_%Y%m%d%H%M.bin".format(runsubdir)
                ),
            )

        if mode != "adj":
            # Saves trajectory files for later adjoint simulations
            list_file = glob.glob("{}/traj*".format(runsubdir))
            for traj_file in list_file:
                rad, ext = os.path.splitext(os.path.basename(traj_file))
                shutil.move(
                    "{}/{}{}".format(runsubdir, rad, ext),
                    date.strftime(
                        "{}/../chain/{}_%Y%m%d%H%M{}".format(
                            runsubdir, rad, ext
                        )
                    ),
                )

    if mode != "adj":
        # Keeps the running directory in memory for later adjoint simulations
        self.adj_refdir = "{}/../".format(runsubdir)


def dump2nc(self, runsubdir):
    """Dumps simulated concentration field to a netCDF file"""

    # Defines reference dates for NetCDF dates
    dref = datetime.datetime.strptime(
        os.path.basename(os.path.normpath(runsubdir)), "%Y-%m-%d_%H-%M"
    )

    # Grid dimension
    nlon = self.domain.nlon
    nlat = self.domain.nlat
    nlev = self.domain.nlev
    ndates = calendar.monthrange(dref.year, dref.month)[1] * 8 + 1

    # Dimensions
    dimnames = ["time", "lev", "lat", "lon"]
    dimlens = [None, nlev, nlat, nlon]
    units = [dref.strftime("hours since %Y-%m-%d %H:%M:%S")]
    vardims = [("time",)]
    varnames = ["time"]
    dtypes = ["d"]
    variables = [3 * np.arange(ndates)]

    # Loop over trajq files to extract concentrations
    for spec in self.chemistry.emis_species.attributes:
        traj_file = "{}/trajq_{}.bin".format(runsubdir, spec)

        field = np.fromfile(traj_file, dtype="float")
        field = field.reshape((-1, nlev, nlat, nlon), order="C")

        vardims.append(("time", "lev", "lat", "lon"))
        varnames.append(spec)
        variables.append(field)
        dtypes.append("d")
        units.append("kg/kg")

    # Add longitudes and latitudes
    vardims.extend([("lat"), ("lon")])
    varnames.extend(["lat", "lon"])
    variables.append(self.domain.zlat[:, 0])
    variables.append(self.domain.zlon[0, :])
    dtypes.extend(2 * ["d"])
    units.extend(["", ""])

    # Saves to NetCDF
    fileout = "{}/trajq.nc".format(runsubdir)
    save_nc(
        fileout, variables, varnames, vardims, dimnames, dimlens, units, dtypes
    )
