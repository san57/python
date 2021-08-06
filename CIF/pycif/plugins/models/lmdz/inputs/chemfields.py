import numpy as np
from builtins import zip
from scipy.io import FortranFile

from pycif.utils import path


def make_chemfields(self, data_in, input_type, ddi, ddf, runsubdir, mode):

    # Deals only with species whose prodloss is activated
    if not hasattr(getattr(self, "chemistry", None), input_type):
        return

    # Binary and nc file name depending on input_type
    bin_name = "prodscale" if input_type == "prodloss3d" else "scale"
    nc_name = "prodloss" if input_type == "prodloss3d" else "prescr"

    datastore = data_in.datastore
    for spec in getattr(self.chemistry, input_type).attributes:
        trid = (input_type, spec)

        if trid not in datastore:
            continue

        data = datastore[trid]

        # Links reference netCDF files that are needed anyway by LMDZ
        try:
            dirorig = datastore[trid]["dirorig"]
            fileorig = datastore[trid]["fileorig"]

            if fileorig is None or dirorig is None:
                raise KeyError

        except KeyError:
            tracer = getattr(getattr(self.chemistry, input_type), spec)

            dirorig = tracer.dir
            fileorig = tracer.file

        origin = ddi.strftime("{}/{}".format(dirorig, fileorig))

        target = "{}/{}_{}.nc".format(runsubdir, nc_name, spec)
        path.link(origin, target)

        # Skip if that species is not in the control vector
        if "spec" not in data:
            continue

        # If tangent-linear mode, include tl increments
        if "incr" in data and mode == "tl":
            incr = data["incr"]
        else:
            incr = 0.0 * data["spec"]

        # Write to FORTRAN binary
        prod_file = "{}/mod_{}_{}.bin".format(runsubdir, bin_name, spec)
        with FortranFile(prod_file, "w") as f:
            # Looping over all values and writing to binary
            prod = data["spec"].values
            incr = data["incr"].values
            for d0, d1 in zip(np.ravel(prod), np.ravel(incr)):
                f.write_record(np.array([d0, d1], dtype=float))
