from __future__ import print_function

import numpy as np

from pycif.utils.check import info
from pycif.utils.datastores.dump import dump_datastore


def execute(self, **kwargs):
    """Runs the model in forward mode

    Args:
        setup (Plugin): definition of the full set-up

    """

    # Working directory
    workdir = self.workdir

    # Control vector
    statevect = self.statevect

    # Observation operator
    obsoper = self.obsoperator

    # Simulation window
    datei = self.datei
    datef = self.datef

    # Some verbose
    info("Running a direct run")

    # Putting x at xb value if available
    if hasattr(statevect, "xb"):
        statevect.x = statevect.xb

    # Running the observation operator
    obsvect = obsoper.obsoper(
        statevect, "fwd", datei=datei, datef=datef, workdir=workdir, **kwargs
    )

    # Perturbs the output monitor if required in the Yaml
    if getattr(self, "perturb_obsvect", False):
        # Altering obsvect and save data
        obserror = self.obserror * obsvect.datastore["sim"].mean()

        obsvect.datastore["obs"] = (
            np.random.normal(
                loc=0, scale=obserror, size=obsvect.datastore.index.size
            )
            + obsvect.datastore["sim"]
        )
        obsvect.datastore["obserror"] = obserror

        # Dumping the datastore with reference data
        dump_datastore(
            obsvect.datastore,
            file_monit=obsvect.file_obsvect,
            dump_type="nc",
            mode="w",
        )

    return obsvect
