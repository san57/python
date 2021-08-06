from __future__ import print_function

import filecmp
import os
import shutil

import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

from pycif.utils import path


def make_fluxes(self, data, runsubdir, datei, mode):
    """Make AEMISSIONS.nc and BEMISSIONS.nc files for CHIMERE.
    Use chemical scheme to check which species is needed and either take it
    from the datastore (i.e. when defined in the control vector), or take it
    from prescribed emissions

    Args:
        self (pycif.utils.classes.Fluxes.fluxes): Flux plugin with all
                    attributes
        datastore (dict): information on flux species
        runsubdir (str): directory to the current run
        nho (int): number of hours in the run
        mode (str): running mode: 'fwd', 'tl' or 'adj'

    """

    datastore = {
        trid: data.datastore[trid]
        for trid in data.datastore
        if trid[0] in ["fluxes", "biofluxes"]
    }
    
    # List of dates for which emissions are needed
    list_dates = pd.date_range(datei, periods=self.nhours + 1, freq="H")

    # Getting the right emissions
    # Loop on all anthropogenic and biogenic species
    # If in datastore, take data, otherwise, link to original A/B EMISSIONS
    list_trid = [("fluxes", spec)
                 for spec in self.chemistry.emis_species.attributes] \
        + [("biofluxes", spec)
           for spec in self.chemistry.bio_species.attributes]
    
    for trid in list_trid:
        spec = trid[1]
        emis_type = trid[0]
        flx_plg = self.fluxes if emis_type == "anthro" else self.biofluxes
        if trid in datastore:
            pass
        
        # If spec not explicitly defined in datastore,
        # fetch general component information if available
        elif trid not in datastore and (emis_type, "") in datastore:
            trid = (emis_type, "")
        else:
            continue

        # Bio or anthro file
        if emis_type == "fluxes":
            file_emisout = "{}/AEMISSIONS.nc".format(runsubdir)
            file_emisincrout = "{}/AEMISSIONS.increment.nc".format(
                runsubdir)
        else:
            file_emisout = "{}/BEMISSIONS.nc".format(runsubdir)
            file_emisincrout = "{}/BEMISSIONS.increment.nc".format(runsubdir)

        tracer = datastore[trid]
        dirorig = tracer["dirorig"]
        fileorig = tracer["fileorig"]
        fileemis = datei.strftime("{}/{}".format(dirorig, fileorig))

        # If no data is provided, just copy from original file
        if "spec" not in tracer:
            linked = False
            # If does not exist, just link
            if not os.path.isfile(file_emisout):
                path.link(fileemis, file_emisout)
                linked = True

            # Otherwise, check for difference
            if not linked:
                if not filecmp.cmp(fileemis, file_emisout):
                    with Dataset(fileemis, "r") as fin:
                        emisin = fin.variables[spec][:]
                        emisin = xr.DataArray(
                            emisin,
                            coords={"time": list_dates},
                            dims=("time", "lev", "lat", "lon"),
                        )
                    flx_plg.write(spec, file_emisout, emisin)

            # Repeat operations for tangent linear
            if mode != "tl":
                continue

            if "spec" not in tracer:
                # If does not exist, just link
                if not os.path.isfile(file_emisincrout):
                    shutil.copy(fileemis, file_emisincrout)

                flx_incr = xr.DataArray(
                    np.zeros(
                        (
                            len(list_dates),
                            self.nlevemis if emis_type == "fluxes"
                            else self.nlevemis_bio,
                            self.domain.nlat,
                            self.domain.nlon,
                        )
                    ),
                    coords={"time": list_dates},
                    dims=("time", "lev", "lat", "lon"),
                )
                flx_plg.write(spec, file_emisincrout, flx_incr)

        else:
            # Replace existing link by copy of original file to modify it
            path.copyfromlink(file_emisout)

            # Put in dataset and write to input
            flx_fwd = datastore[trid]["spec"]
            flx_plg.write(spec, file_emisout, flx_fwd)

            if mode == "tl":
                path.copyfromlink(file_emisincrout)
                flx_tl = datastore[trid].get("incr", 0.0 * flx_fwd)
                flx_plg.write(spec, file_emisincrout, flx_tl)
