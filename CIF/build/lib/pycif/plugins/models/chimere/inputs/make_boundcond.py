from __future__ import print_function

import os
import shutil

from netCDF4 import Dataset

from pycif.utils import path


def make_boundcond(self, data, runsubdir, ddi, mode, input_type):
    """
    Generates boundary conditions files for CHIMERE

    :param self:
    :param datastore:
    :type datastore: dict
    :param runsubdir:
    :type runsubdir: str
    :param sdc:
    :type sdc: str
    :param hour_dates:
    :param mode:
    :param input_type:
    :return:
    """

    datastore = {
        trid: data.datastore[trid]
        for trid in data.datastore
        if trid[0] == input_type
    }

    # Name of variable in netCDF
    nc_varname = "top_conc" if input_type == "topcond" else "lat_conc"

    # Fixed name for BC files
    fileout = "{}/BOUN_CONCS.nc".format(runsubdir)
    fileoutincr = "{}/BOUN_CONCS.increment.nc".format(runsubdir)

    # Loop on all active species
    # If in datastore, take data, otherwise, link to original INI_CONCS
    for spec in self.chemistry.acspecies.attributes:
        trid = (input_type, spec)
        if trid in datastore:
            pass

        # If spec not explicitly defined in datastore,
        # fetch general component information if available
        elif trid not in datastore and (input_type, "") in datastore:
            trid = (input_type, "")
        else:
            continue

        tracer = datastore[trid]
        dirorig = tracer["dirorig"]
        fileorig = tracer["fileorig"]
        fileini = ddi.strftime("{}/{}".format(dirorig, fileorig))

        # If no data is provided, just copy from original file
        if "spec" not in tracer:
            # If does not exist, just link
            if not os.path.isfile(fileout):
                path.link(fileini, fileout)

            # # Otherwise, check for difference
            # elif not filecmp.cmp(fileini, fileout):
            #
            #     print(__file__)
            #     import code
            #     code.interact(local=dict(locals(), **globals()))
            #     raise Exception('I need to transfer lbc '
            #                     'from one file to the other one')

        else:
            # Replace existing link by copy
            # of original file to modify it
            path.copyfromlink(fileout)

            # Write initial conditions
            lbc_fwd = datastore[trid]["spec"]
            self.latcond.write(spec, fileout, lbc_fwd, comp_type=input_type)

            if mode == "tl":
                path.copyfromlink(fileoutincr)
                lbc_tl = datastore[trid].get("incr", 0.0 * lbc_fwd)
                self.latcond.write(
                    spec, fileoutincr, lbc_tl, comp_type=input_type
                )

        # Repeat operations for tangent linear
        if mode != "tl":
            continue

        if "spec" not in tracer:
            # If does not exist, just link
            if not os.path.isfile(fileoutincr):
                shutil.copy(fileini, fileoutincr)

            with Dataset(fileoutincr, "a") as fout:
                nc_var_out = fout.variables[nc_varname][:]
                nc_names = [
                    str(b"".join(s).strip().lower(), "utf-8")
                    for s in fout.variables["species"][:]
                ]
                if spec.lower() not in nc_names:
                    continue

                # Apply to original data
                nc_var_out[..., nc_names.index(spec.lower())][:] = 0.0

                fout.variables[nc_varname][:] = nc_var_out
