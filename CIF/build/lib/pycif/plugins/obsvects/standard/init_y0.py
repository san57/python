import os
import shutil

from builtins import str

from pycif.plugins.obsvects.standard.utils.crop_monitor import crop_monitor
from pycif.plugins.obsvects.standard.utils.hcoord import hcoord
from pycif.plugins.obsvects.standard.utils.tstep import tstep
from pycif.plugins.obsvects.standard.utils.vcoord import vcoord
from pycif.utils import path
from pycif.utils.check import info
from pycif.utils.check.errclass import PluginError
from pycif.utils.datastores.dump import dump_datastore, read_datastore
from pycif.utils.datastores.empty import init_empty


def init_y0(obsvect, measurements, **kwargs):
    """Initializes the observation vector. In most cases the observation
    vector is similar to the measurement vector but there is no reason
    for it to be the same, especially when super-observations are used
    (e.g. daily or afternoon averages, gradients, etc.)

    Args:
        obsvect (Plugin): the observation vector with all its attributes
        measurements (Plugin): the pre-loaded measurements

    Returns:
        obsvect, updated with horizontal and vertical coordinates,
        as well as model time steps

    """

    # Saves initialized observation vector to workdir/obsvect/monitor
    # Try linking file_obsvect to workdir/obsvect if existing
    # Otherwise, define it and generate it later on
    path.init_dir("{}/obsvect".format(obsvect.workdir))
    file_monitor = "{}/obsvect/monitor.{}".format(
        obsvect.workdir, obsvect.dump_type
    )
    shutil.rmtree(file_monitor, ignore_errors=True)

    if os.path.isfile(obsvect.file_obsvect):
        path.link(obsvect.file_obsvect, file_monitor)

    # From here, no interaction with files outside the working directory
    obsvect.file_obsvect = file_monitor

    # Try reading file
    allcorrec, ok_hcoord, ok_vcoord, do_tstep = False, False, False, True
    exclude_stations = getattr(obsvect, "exclude_stat", [])
    try:
        # Reading a monitor file if available
        info("Try opening {}".format(file_monitor))
        obsvect.datastore = read_datastore(
            file_monitor, dump_type=obsvect.dump_type
        )

        # Check that the monitor is consistent with the simulation to be run
        # If True, returns directly the observation as is
        allcorrect, ok_hcoord, ok_vcoord, do_tstep = obsvect.check_monitor()
        if allcorrect and len(exclude_stations) == 0:
            return obsvect

        # Otherwise, recompute necessary items (i, j, lev, tstep)
        # of the observation error
        raise PluginError

    # If the monitor file could not be opened, start back from measurements
    except IOError:

        info("Couldn't open the observation file")
        info("Generating observation vector from measurements")

        # Copying attributes from measurements
        for key in measurements.attributes:
            setattr(obsvect, key, getattr(measurements, key))

        # Copying datastore from measurements if existing and not empty
        obsvect.datastore = getattr(measurements, "datastore", init_empty())

    except PluginError:
        info(
            "Warning! The prescribed monitor file is not consistent with "
            "the current simulation. Re-analysing it"
        )

    # Remove unwanted stations
    if len(exclude_stations) > 0:
        mask = ~(obsvect.datastore["station"].isin(exclude_stations))
        obsvect.datastore = obsvect.datastore.loc[mask]

    # Crops y0 to the simulation period
    obsvect.datastore = crop_monitor(
        obsvect.datastore, obsvect.datei, obsvect.datef, **kwargs
    )

    # Computes grid cells where observations fall
    if not ok_hcoord:
        obsvect = hcoord(obsvect, **kwargs)

    # Computes model layer where observations fall
    if not ok_vcoord:
        obsvect = vcoord(obsvect, **kwargs)

    # Compute time steps corresponding to the observations
    if do_tstep:
        obsvect = tstep(obsvect, **kwargs)

    # Dumps the datastore
    # Attributes to keep in the monitor if NetCDF format (recommanded!)
    nc_attributes = {
        "datei": obsvect.datei.strftime("%d-%m-%Y %H:%M:%S"),
        "datef": obsvect.datef.strftime("%d-%m-%Y %H:%M:%S"),
        "model name": obsvect.model.plugin.name,
        "model version": obsvect.model.plugin.version,
        "domain nlon": str(obsvect.model.domain.nlon),
        "domain nlat": str(obsvect.model.domain.nlat),
    }

    if getattr(obsvect, "dump", True):
        dump_datastore(
            obsvect.datastore,
            file_monit=obsvect.file_obsvect,
            dump_type=obsvect.dump_type,
            nc_attributes=nc_attributes,
            mode="w",
        )

    return obsvect
