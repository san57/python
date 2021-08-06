import glob
import os
from types import MethodType

import numpy as np

from pycif.plugins.obsvects.standard.utils.check_monitor import check_monitor
from pycif.utils import path
from .init_rinvprod import init_rinvprod
from .init_y0 import init_y0
from .native2obsvect import native2obsvect
from .obsvect2native import obsvect2native
from .rinvprod import rinvprod

# It is necessary to have some measurements and some info about the meteo
# to initialize the observation vector
requirements = {
    "measurements": {"any": True, "empty": True},
    "model": {"any": True, "empty": False},
}


def ini_data(plugin, **kwargs):
    """Initializes the observation vector from information in the Yaml file

    Args:

    """

    # Set dump type if not defined; default is nc
    if not hasattr(plugin, "dump_type"):
        plugin.dump_type = "nc"

    # Set default file_obsvect
    file_default = "{}/obsvect/monitor.{}".format(
        plugin.workdir, plugin.dump_type
    )
    plugin.file_obsvect = getattr(plugin, "file_obsvect", file_default)

    # Keeping check_monitor as a class method
    plugin.check_monitor = MethodType(check_monitor, plugin)

    # Initializing y0
    measurements = plugin.measurements
    plugin.dump_type = getattr(plugin, "dump_type", "nc")
    init_y0(plugin, measurements, **kwargs)

    # Initialize R if any observation
    if plugin.datastore.size > 0:
        init_rinvprod(plugin, measurements, **kwargs)

    print("Link satellite files to the working directory")
    print(hasattr(plugin, "dir_satellites"), plugin.workdir)
    if hasattr(plugin, "dir_satellites"):
        path.init_dir("{}/obsvect/satellites/".format(plugin.workdir))
        for dd in plugin.model.subsimu_dates[:-1]:
            sat_files = glob.glob(
                dd.strftime("{}/infos_*%Y%m%d%H%M.nc").format(
                    plugin.dir_satellites
                )
            )

            for obs_file in sat_files:
                target = "{}/obsvect/satellites/{}".format(
                    plugin.workdir, os.path.basename(obs_file)
                )
                path.link(obs_file, target)

    plugin.has_satellites = False
    if np.where((plugin.datastore["level"] < 0))[0].size > 0:
        plugin.has_satellites = True
