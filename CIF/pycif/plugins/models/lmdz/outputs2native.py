import datetime
import glob
import os

import numpy as np
import xarray as xr

from pycif.utils.datastores.empty import init_empty


def outputs2native(
    self, data2dump, input_type, di, df, runsubdir, mode="fwd", dump=True
):
    """Reads outputs to pycif objects.

    If the mode is 'fwd' or 'tl', only observation-like outputs are extracted.
    For the 'adj' mode, all outputs relative to model sensitivity are extracted.

    Dumps to a NetCDF file with output concentrations if needed

    Args:
        self (pycif.utils.classes.models.Model): Model object
        runsubdir (str): current sub-sumilation directory
        mode (str): running mode; one of: 'fwd', 'tl', 'adj'; default is 'fwd'

        dump (bool): dumping outputs or not; default is True
    Return:
        dict

    """

    ddi = min(di, df)
    ddf = max(di, df)

    if mode in ["tl", "fwd"]:
        if not hasattr(self, "dataobs"):
            self.dataobs = init_empty()

        # Read simulated concentrations
        sim_file = "{}/obs_out.bin".format(runsubdir)
        if not os.path.isfile(sim_file):
            self.dataobs.loc[:, "sim"] = np.nan
            return

        sim = np.fromfile(sim_file, dtype="float").reshape((-1, 4), order="F")

        # Observations that were not extracted by LMDZ are set to NaN
        sim[sim == 0] = np.nan

        # Putting values to the local data store
        self.dataobs.loc[:, "sim"] = sim[:, 0]
        self.dataobs["pressure"] = sim[:, 2]
        self.dataobs["dp"] = sim[:, 3]

        if mode == "tl":
            self.dataobs.loc[:, "sim_tl"] = sim[:, 1]

        return self.dataobs

    elif mode == "adj":
        nlon = self.domain.nlon
        nlat = self.domain.nlat

        # Stores daily dates of the period for later aggregation
        dref = datetime.datetime.strptime(
            os.path.basename(os.path.normpath(runsubdir)), "%Y-%m-%d_%H-%M"
        )
        list_dates = self.input_dates[ddi]

        # Reading only output files related to given input_type
        ref_names = {
            "inicond": "init",
            "fluxes": "fluxes",
            "prescrconcs": "scale",
            "prodloss3d": "prodscale",
        }

        if input_type not in ref_names:
            return data2dump

        list_file = glob.glob(
            "{}/mod_{}_*_out.bin".format(runsubdir, ref_names[input_type])
        )
        specs2dump = [s[1] for s in data2dump]
        for out_file in list_file:
            spec = os.path.basename(out_file).split("_")[2]

            if spec not in specs2dump:
                continue

            with open(out_file, "rb") as f:
                data = np.fromfile(f, dtype=np.float)

            data = data.reshape((nlon, nlat, -1), order="F").transpose(
                (2, 1, 0)
            )

            if input_type == "inicond":
                data2dump[("inicond", spec)]["adj_out"] = xr.DataArray(
                    data[np.newaxis, ...],
                    coords={"time": np.array([dref])},
                    dims=("time", "lev", "lat", "lon"),
                )
                continue

            # Adding one time stamp to fit with input dates
            # including the first stamp of the next month
            data = np.concatenate((data, np.zeros((1, nlat, nlon))), axis=0)
            data = data[:, np.newaxis, ...]
            data2dump[(input_type, spec)]["adj_out"] = xr.DataArray(
                data,
                coords={"time": list_dates},
                dims=("time", "lev", "lat", "lon"),
            )

        return data2dump
