import numpy as np
import xarray as xr
from future import standard_library

from pycif.utils.path import init_dir
from .utils.scalemaps import scale2map

standard_library.install_aliases()
try:
    import cPickle as pickle
except ImportError:
    import pickle


def dump(self, cntrl_file, to_netcdf=False, dir_netcdf=None, **kwargs):
    """Dumps a control vector into a pickle file.
    Does not save large correlations.

    Args:
        self (pycif.utils.classes.controlvects.ControlVect):
                    the Control Vector to dump
        file (str): path to the file to dump as pickle
        to_netcdf (bool): save to netcdf files if True
        dir_netcdf (str): root path for the netcdf directory
    """

    # Saving recursive attributes from the Yaml
    exclude = ["transform", "domain", "datastore",
               "input_dates", "obsvect", "tracer",
               "tstep_dates", "tstep_all", "dataflx"]
    tosave = self.to_dict(self, exclude_patterns=exclude)

    # Adding control vector values
    if hasattr(self, "x"):
        tosave["x"] = self.x

    if hasattr(self, "dx"):
        tosave["dx"] = self.dx

    if hasattr(self, "xb"):
        tosave["xb"] = self.xb
    
    # Dumping the dictionary to a pickle
    with open(cntrl_file, "wb") as f:
        pickle.dump(tosave, f, pickle.HIGHEST_PROTOCOL)

    # Dumping to an ensemble of NetCDF files
    if not to_netcdf or dir_netcdf is None:
        return

    components = self.components
    for comp in components.attributes:
        component = getattr(components, comp)

        dir_comp = "{}/{}".format(dir_netcdf, comp)
        init_dir(dir_comp)

        for trcr in component.parameters.attributes:
            tracer = getattr(component.parameters, trcr)

            # Translating x and xb to maps
            x = np.reshape(
                self.x[tracer.xpointer: tracer.xpointer + tracer.dim],
                (tracer.ndates, -1),
            )
            x = scale2map(x, tracer, tracer.dates, self.domain)

            xb = np.reshape(
                self.xb[tracer.xpointer: tracer.xpointer + tracer.dim],
                (tracer.ndates, -1),
            )
            xb = scale2map(xb, tracer, tracer.dates, self.domain)

            std = np.reshape(
                self.std[tracer.xpointer: tracer.xpointer + tracer.dim],
                (tracer.ndates, -1),
            )
            std = scale2map(std, tracer, tracer.dates, self.domain)

            ds = xr.Dataset({"x": x, "xb": xb, "std": std})
            ds.to_netcdf("{}/statevect_{}_{}.nc".format(dir_comp, comp, trcr))


def load(self, cntrl_file, **kwargs):
    with open(cntrl_file, "r") as f:
        toread = pickle.load(f)

    out = self.from_dict(toread)

    if hasattr(out, "x"):
        self.x = out.x

    if hasattr(out, "dx"):
        self.dx = out.dx

    if hasattr(out, "xb"):
        self.xb = out.xb

    return out
