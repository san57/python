import os
import shutil

from netCDF4 import Dataset

from pycif.utils import path


def make_inicond(self, data_in, datei, datef, runsubdir, mode):

    ddi = min(datei, datef)
    ddf = max(datei, datef)

    # Reference initial condition file for looping sub-simulations
    if hasattr(self, "chain"):
        source = self.chain.strftime(
            "{}/../chain/restart_%Y%m%d%H%M.nc".format(runsubdir)
        )
        target = "{}/start.nc".format(runsubdir)
        path.link(source, target)

        if mode == "tl":
            source = self.chain.strftime(
                "{}/../chain/restart_tl_%Y%m%d%H%M.bin".format(runsubdir)
            )
            target = "{}/start_tl.bin".format(runsubdir)
            path.link(source, target)

        return

    # Generating reference initial conditions if first sub-simulation
    datastore = data_in.datastore

    for spec in self.chemistry.acspecies.attributes:
        restartID = getattr(self.chemistry.acspecies, spec).restart_id
        if ("inicond", spec) not in datastore:
            continue

        ini = datastore[("inicond", spec)]

        ref_dir = ini["dirorig"]
        ref_file = ini["fileorig"]

        # Copy if does not already exists
        source = ddi.strftime("{}/{}".format(ref_dir, ref_file))
        target = "{}/start.nc".format(runsubdir)
        if not os.path.isfile(target):
            shutil.copy(source, target)

        if mode == "tl":
            target_tl = "{}/start_tl.nc".format(runsubdir)
            if not os.path.isfile(target_tl):
                shutil.copy(source, target_tl)

        # Updates data
        if mode in ["fwd", "tl"] and ddi == self.datei:
            var = "q{:02d}".format(restartID)

            if "spec" in ini:
                data = ini["spec"]
                with Dataset(target, "a") as f:
                    if var in f.variables:
                        f.variables[var][:] = data.values
                    else:
                        self.inicond.write(var, target, data)

            if mode == "tl":
                if "incr" in ini:
                    data = ini["incr"]
                    with Dataset(target_tl, "a") as f:
                        if var in f.variables:
                            f.variables[var][:] = data.values
                        else:
                            self.inicond.write(var, target_tl, data)

                else:
                    with Dataset(target_tl, "a") as f:
                        if var in f.variables:
                            f.variables[var][:] = 0.0

        elif mode == "adj" and ddf == self.datef:
            for spec in self.chemistry.acspecies.attributes:
                restartID = getattr(self.chemistry.acspecies, spec).restart_id
                with Dataset(target, "a") as f:
                    var = "q{:02d}".format(restartID)
                    if var in f.variables:
                        f.variables[var][:] = 0.0
