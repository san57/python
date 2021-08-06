import datetime
import glob
import os

import numpy as np
import pandas as pd
import xarray as xr
from builtins import zip
from netCDF4 import Dataset

from pycif.utils.check import info
from pycif.utils.datastores.empty import init_empty


def outputs2native(
    self, data2dump, input_type, di, df, runsubdir, mode="fwd", dump=True
):
    """Reads outputs to pyCIF objects.

    If the mode is 'fwd' or 'tl', only observation-like outputs are extracted.
    For the 'adj' mode, all outputs relative to model sensitivity are extracted.

    Dumps to a NetCDF file with output concentrations if needed"""

    ddi = min(di, df)

    if mode in ["fwd", "tl"]:
        if not hasattr(self, "dataobs"):
            self.dataobs = init_empty()

        # If no simulated concentration is available just pass
        sim_file = "{}/mod.txt".format(runsubdir)
        print("runsubdir", runsubdir)
        if os.stat(sim_file).st_size == 0:
            info(
                "CHIMERE ran without any observation "
                "to be compared with for sub-simu "
                "only CHIMERE's outputs are available"
            )
            self.dataobs.loc[:, "sim"] = np.nan
            return

        # Read simulated concentrations
        data = pd.read_csv(
            sim_file,
            delim_whitespace=True,
            header=None,
            usecols=range(6, 12),
            names=["sim", "pmid", "dp", "airm", "hlay", "simfwd"],
        )

        # Loop over observations in active species
        mask = (
            self.dataobs["parameter"]
            .str.upper()
            .isin(self.chemistry.acspecies.attributes)
        )

        # Putting values to the local data store
        # Assumes arithmetic averages upto now
        inds = [0] + list(np.cumsum(self.dataobs.loc[mask, "dtstep"][:-1]))

        column = "sim" if mode == "fwd" else "sim_tl"

        dataavg = pd.DataFrame(
            [
                data.iloc[k: k + dt].sum()
                for k, dt in zip(inds, self.dataobs.loc[mask, "dtstep"])
            ]
        )

        self.dataobs.loc[mask, column] = dataavg.loc[:, "sim"].values
        self.dataobs.loc[mask, "pressure"] = dataavg.loc[:, "pmid"].values
        self.dataobs.loc[mask, "dp"] = dataavg.loc[:, "dp"].values
        self.dataobs.loc[mask, "airm"] = dataavg.loc[:, "airm"].values
        self.dataobs.loc[mask, "hlay"] = dataavg.loc[:, "hlay"].values

        if column == "sim_tl":
            self.dataobs.loc[mask, "sim"] = dataavg.loc[:, "simfwd"].values

        return self.dataobs

    elif mode == "adj":
        # List of CHIMERE dates
        dref = datetime.datetime.strptime(
            os.path.basename(os.path.normpath(runsubdir)), "%Y-%m-%d_%H-%M"
        )
        list_dates = self.input_dates[ddi]

        # Reading only output files related to given input_type
        ref_names = {
            "inicond": "ini",
            "fluxes": "aemis",
            "biofluxes": "bemis",
            "latcond": "bc",
            "topcond": "bc",
        }

        if input_type not in ref_names:
            return data2dump

        list_file = glob.glob(
            "{}/aout.*{}*.nc".format(runsubdir, ref_names[input_type])
        )

        for out_file in list_file:
            with Dataset(out_file, "r") as f:
                # Load list of species and reformat it
                list_species = [
                    b"".join(s).strip().decode("ascii")
                    for s in f.variables["species"][:]
                ]

                # Restrict to species required in data2dump
                dump_species = [
                    t[1].lower() for t in data2dump if t[0] == input_type
                ]
                if "" not in dump_species:
                    list_species = [
                        s for s in list_species if s.lower() in dump_species
                    ]

                else:
                    for s in list_species:
                        if s not in data2dump:
                            data2dump[(input_type, s)] = {}

                # Different output structure between LBC and others
                if "ini" in out_file or "emis" in out_file:
                    data = {s: f.variables[s][:] for s in list_species}

                elif "bc" in out_file:
                    data_lat = {
                        s: f.variables["lat_conc"][..., k]
                        for k, s in enumerate(list_species)
                    }
                    data_top = {
                        s: f.variables["top_conc"][..., k]
                        for k, s in enumerate(list_species)
                    }

            if "ini" in out_file:
                for spec in data:
                    data2dump[("inicond", spec)]["adj_out"] = xr.DataArray(
                        data[spec][np.newaxis, ...],
                        coords={"time": np.array([dref])},
                        dims=("time", "lev", "lat", "lon"),
                    )

            elif "bc" in out_file:
                if input_type == "latcond":
                    for spec in data_lat:
                        data2dump[("latcond", spec)]["adj_out"] = xr.DataArray(
                            data_lat[spec][..., np.newaxis, :],
                            coords={"time": list_dates},
                            dims=("time", "lev", "lat", "lon"),
                        )
                if input_type == "topcond":
                    for spec in data_top:
                        data2dump[("topcond", spec)]["adj_out"] = xr.DataArray(
                            data_top[spec][..., np.newaxis, :],
                            coords={"time": list_dates},
                            dims=("time", "lev", "lat", "lon"),
                        )

            elif "aemis" in out_file or "bemis" in out_file:
                if "aemis" in out_file:
                    emis_type = "fluxes"
                else:
                    emis_type = "biofluxes"
                    
                for spec in data:
                    data2dump[(emis_type, spec)]["adj_out"] = xr.DataArray(
                        data[spec],
                        coords={"time": list_dates},
                        dims=("time", "lev", "lat", "lon"),
                    )

        return data2dump
