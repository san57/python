import numpy as np


def make_fluxes(self, data, ddi, ddf):

    datastore = data.datastore
    for spec in datastore:
        ds = datastore[spec]

        flx = datastore[spec]["spec"]

        if "scale" in datastore[spec]:
            flx *= datastore[spec]["scale"]

        # Keeping only data in the simulation window
        mask = (flx.time >= np.datetime64(ddi)) & (
            flx.time <= np.datetime64(ddf)
        )
        mask = np.where(mask)[0]

        flx_fwd = flx[mask]
        flx_tl = 0.0 * flx_fwd
        if "incr" in datastore[spec]:
            flx_tl = datastore[spec]["incr"][mask]

        # Saving data in model object
        # For other models, intermediate NetCDF or binary files should be
        # saved for later use by numerical Fortran (or other) model
        self.dataflx = flx_fwd
        self.dataflx_tl = flx_tl
