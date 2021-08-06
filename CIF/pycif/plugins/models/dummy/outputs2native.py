from pycif.utils.datastores.empty import init_empty


def outputs2native(
    self, data2dump, input_type, di, df, runsubdir, mode="fwd", dump=True
):
    """Reads outputs to pycif objects.

    If the mode is 'fwd' or 'tl', only onservation-like outputs are extracted.
    For the 'adj' mode, all outputs relative to model sensitivity are extracted.

    Dumps to a NetCDF file with output concentrations if needed"""

    if mode in ["tl", "fwd"]:
        if not hasattr(self, "dataobs"):
            self.dataobs = init_empty()

        # Read simulated concentrations
        # In classical model, this should correspond to reading output files
        # Here the observations are already stored in the model object
        datastore = self.dataobs

        # Re-aggregate observations spanning several time steps
        # Obsvect divides by number of tstep at higher level
        # (in case the observation spans several periods)
        ds = datastore.groupby([datastore.index, "station", "i", "j"]).sum()
        ds = ds[["sim", "sim_tl"]]
        ds.index = ds.index.get_level_values(0)

        self.dataobs = ds

        return self.dataobs

    elif mode == "adj":
        if input_type != "fluxes":
            return data2dump

        # Reads sensitivities
        # In the toy model's case, just take the data from the object itself
        datasensit = self.dflx

        # TODO: generalize with several species
        spec = self.chemistry.acspecies.attributes[0]
        data2dump[("fluxes", spec)]["adj_out"] = datasensit

        return data2dump
