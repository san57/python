#!/usr/bin/env python
# -*- coding: utf-8 -*-
import xarray as xr

from pycif.utils.check import info


def make_fluxes(self, data, ddi, ddf, runsubdir, mode):
    """Prepare a binary file per emitted species containing flux data.

    :param self:
    :param datastore:
    :param ddi:
    :param ddf:
    :param runsubdir:
    :param mode:
    :return:
    """

    datastore = data.datastore

    for spec in self.chemistry.emis_species.attributes:
        tracer = getattr(self.chemistry.emis_species, spec)

        if not ("fluxes", spec) in datastore:
            info("{} not available for being emitted in LMDZ".format(spec))
            continue

        info("LMDZ is generating flux inputs for {}".format(spec))

        data = datastore[("fluxes", spec)]

        # If not determined by the control vector
        if "spec" not in data:
            data["spec"] = self.fluxes.read(
                spec,
                data["dirorig"],
                data["fileorig"],
                data["varname"],
                self.input_dates[ddi],
            )

        # Adds empty increments if not available
        if "incr" not in data:
            data["incr"] = 0.0 * data["spec"]

        # Put in dataset for writing by 'write'
        ds = xr.Dataset({"fwd": data["spec"], "tl": data["incr"]})

        # Write to FORTRAN binary
        flx_file = "{}/mod_{}.bin".format(runsubdir, spec)
        self.emis_species.write(spec, flx_file, ds)
