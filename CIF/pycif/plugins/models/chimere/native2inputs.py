from __future__ import absolute_import
from __future__ import print_function

import pandas as pd

from pycif.utils import path
from .inputs.make_boundcond import make_boundcond
from .inputs.make_fluxes import make_fluxes
from .inputs.make_inicond import make_inicond
from .inputs.make_meteo import make_meteo
from .inputs.make_obs import make_obs
from .inputs.params import make_nml


def native2inputs(
    self, datastore, input_type, datei, datef, runsubdir, mode="fwd"
):
    """Converts data at the model data resolution to model compatible input
    files.

    Args:
        self: the model Plugin
        input_type (str): one of 'fluxes'
        datastore: data to convert
            if input_type == 'fluxes',
        datei, datef: date interval of the sub-simulation
        mode (str): running mode: one of 'fwd', 'adj' and 'tl'
        runsubdir (str): sub-directory for the current simulation
        workdir (str): the directory of the whole pyCIF simulation

    Notes:
        - CHIMERE expects hourly inputs;

    """

    ddi = min(datei, datef)
    ddf = max(datei, datef)

    # Hour steps of the sub-run
    hour_dates = pd.date_range(ddi, ddf, freq="1H")

    if datastore is None:
        datastore = {}

    sdc = ddi.strftime("%Y%m%d%H")

    # Choose the right executable
    if input_type == "exe":
        if mode == "fwd":
            source = "fwdchimere.e"

        elif mode == "tl":
            source = "tlchimere.e"

        elif mode == "adj":
            source = "achimere.e"

        else:
            raise Exception(
                "Unknown mode '{}' for executable with CHIMERE".format(mode)
            )

        path.link(
            "{}/model/{}".format(self.workdir, source),
            "{}/chimere.e".format(runsubdir),
        )

    # Deals with the nml file
    if input_type == "nml":
        make_nml(self, runsubdir, sdc, mode)

    # Deals with fluxes
    # WARNING: Not specified emissions are set to zero??
    if input_type in ["fluxes", "biofluxes"]:
        make_fluxes(self, datastore, runsubdir, ddi, mode)

    if input_type in ["latcond", "topcond"]:
        make_boundcond(self, datastore, runsubdir, ddi, mode, input_type)

    if input_type == "inicond":
        make_inicond(self, datastore, runsubdir, mode, ddi)

    if input_type == "meteo":
        make_meteo(self, runsubdir, sdc)

    # Deals with observations
    if input_type == "obs":
        make_obs(self, datastore, runsubdir, mode)

    return
