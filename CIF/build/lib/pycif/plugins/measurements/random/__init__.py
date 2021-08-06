import numpy as np
import pandas as pd
from builtins import range

from pycif.utils.check import info
from pycif.utils.datastores.dump import dump_datastore, read_datastore
from pycif.utils.datastores.empty import init_empty

requirements = {"domain": {"any": True, "empty": False}}


def parse_tracers(self, datei, datef, file_monitor="", workdir="", **kwargs):
    """Generate random observations at random locations in the domain"""

    # If file_monitor is defined, tries reading it
    if hasattr(self, "file_monitor"):
        file_monitor = self.file_monitor

        try:
            info("Extracting measurements from {}".format(file_monitor))
            return read_datastore(
                file_monitor,
                dump_type=getattr(self, "dump_type", "nc"),
                **kwargs
            )

        except IOError as e:
            info(str(e))
            info("Could not find a monitor file, reading observations")

        except Exception as e:
            info(e)
            info("Could not read the specified monitor file: {}", file_monitor)
            raise e

    # Otherwise, create the monitor from observations
    else:
        if hasattr(self, "workdir"):
            workdir = self.workdir

        if file_monitor == "" or file_monitor is None:
            file_monitor = workdir + "/obs/monit_standard.nc"

    # Initialize empty datastore
    ds = init_empty()

    # Loop over species
    for trcr in self.species.attributes:
        tracer = getattr(self.species, trcr)

        # Make date range
        drange = pd.date_range(
            datei, datef, freq=getattr(tracer, "frequency", "1H")
        )
        ndates = drange.size

        # Pick random locations for x and y within the reference domain
        xmin = self.domain.zlon.min()
        xmax = self.domain.zlon.max()
        ymin = self.domain.zlat.min()
        ymax = self.domain.zlat.max()
        zmax = getattr(tracer, "zmax", 100)

        nstat = tracer.nstations

        statx = np.random.uniform(low=xmin, high=xmax, size=nstat)
        staty = np.random.uniform(low=ymin, high=ymax, size=nstat)
        statz = np.random.uniform(low=1, high=zmax, size=nstat)

        # Put locations into a monitor
        duration = pd.to_timedelta(getattr(tracer, "duration", "1H"))
        seconds_duration = duration.total_seconds() / 3600.0

        df = pd.DataFrame(
            {
                "alt": np.array(ndates * list(statz)),
                "lat": np.array(ndates * list(staty)),
                "lon": np.array(ndates * list(statx)),
                "station": ndates * list(range(nstat)),
                "parameter": trcr,
                "duration": seconds_duration,
            }
        )

        df.index = (
            np.array(nstat * list(drange))
            .reshape((ndates, nstat), order="F")
            .flatten()
        )
        if getattr(tracer, "random_subperiod_shift", False):
            df.index += np.random.uniform(0, 1, size=ndates * nstat) * duration

        # Appending datastore
        ds = pd.concat((ds, df), sort=False)

    # Dumping the datastore for later use by pycif
    dump_datastore(ds, file_monit=self.file_monitor, dump_type="nc")

    return ds
