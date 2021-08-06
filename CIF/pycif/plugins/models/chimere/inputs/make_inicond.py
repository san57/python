from __future__ import print_function

import filecmp
import os
import shutil

import numpy as np
from netCDF4 import Dataset

from pycif.utils import path


def make_inicond(self, data, runsubdir, mode, datei):

    datastore = {
        trid: data.datastore[trid]
        for trid in data.datastore
        if trid[0] == "inicond"
    }

    # Fixed name for INI_CONCS files
    fileout = "{}/INI_CONCS.0.nc".format(runsubdir)
    fileoutincr = "{}/INI_CONCS.0.increment.nc".format(runsubdir)
    nho = self.nho

    # Getting the forward initial concentrations (needed even for adjoint)
    # if chained period
    if hasattr(self, "chain"):
        if mode in ["tl", "fwd"]:
            filein = self.chain.strftime(
                "{}/../chain/end.%Y%m%d%H.{}.nc".format(runsubdir, nho)
            )

        else:
            subsimu_dates = self.subsimu_dates
            date_index = np.where(subsimu_dates == datei)[0][0]
            refdir = self.adj_refdir
            filein = subsimu_dates[date_index - 1].strftime(
                "{}/chain/end.%Y%m%d%H.{}.nc".format(refdir, nho)
            )
        path.link(filein, fileout)

        # Adjoint needs to chain sensitivities as well
        if mode == "adj":
            filein = self.chain.strftime(
                "{}/../chain/aend.%Y%m%d%H.{}.nc".format(runsubdir, nho)
            )
            fileout = self.chain.strftime("{}/aini.nc".format(runsubdir))
            path.link(filein, fileout)

    # Exit if not first period
    if datei != self.datei:
        return

    # Loop on all active species
    # If in datastore, take data, otherwise, link to original INI_CONCS
    for spec in self.chemistry.acspecies.attributes:
        trid = ("inicond", spec)
        if trid in datastore:
            pass
        # If spec not explicitly defined in datastore,
        # fetch general component information if available
        elif trid not in datastore and ("inicond", "") in datastore:
            trid = ("inicond", "")
        else:
            continue

        tracer = datastore[trid]
        dirorig = tracer["dirorig"]
        fileorig = tracer["fileorig"]
        fileini = "{}/{}".format(dirorig, fileorig)

        # If no data is provided, just copy from original file
        if "spec" not in tracer:
            # If does not exist, just link
            linked = False
            if not os.path.isfile(fileout):
                path.link(fileini, fileout)
                linked = True

            # Otherwise, check for difference
            if not linked:
                if not filecmp.cmp(fileini, fileout):
                    raise Exception(
                        "I need to transfer inicond "
                        "from one file to the other one"
                    )

            # Repeat operations for tangent linear
            if mode != "tl":
                continue

            # If does not exist, just link
            if not os.path.isfile(fileoutincr):
                shutil.copy(fileini, fileoutincr)

            with Dataset(fileoutincr, "a") as fout:
                if spec in fout.variables:
                    fout.variables[spec][:] = 0.0

        else:
            # Replace existing link by copy
            # of original file to modify it
            path.copyfromlink(fileout)

            # Write initial conditions
            ini_fwd = datastore[trid]["spec"]
            self.inicond.write(spec, fileout, ini_fwd, comp_type="inicond")

            if mode == "tl":
                path.copyfromlink(fileoutincr)
                ini_tl = datastore[trid].get("incr", 0.0 * ini_fwd)
                self.inicond.write(
                    spec, fileoutincr, ini_tl, comp_type="inicond"
                )
