from __future__ import absolute_import

from .inputs.chemfields import make_chemfields
from .inputs.fluxes import make_fluxes
from .inputs.inicond import make_inicond
from .inputs.obs import make_obs


def native2inputs(
    self, datastore, input_type, datei, datef, runsubdir, mode="fwd"
):
    """Converts data at the model data resolution to model compatible input
    files.

    Args:
        self: the model Plugin
        input_type (str): one of 'fluxes', 'obs'
        datastore: data to convert
            if input_type == 'fluxes', a dictionary with flux maps
            if input_type == 'obs', a pandas dataframe with the observations
        datei, datef: date interval of the sub-simulation
        mode (str): running mode: one of 'fwd', 'adj' and 'tl'
        runsubdir (str): sub-directory for the current simulation
        workdir (str): the directory of the whole pycif simulation

    Notes:
        - LMDZ expects daily inputs; if the periods in the control vector are
        longer than one day, period values are uniformly de-aggregated to the
        daily scale; this is done with pandas function 'asfreq' and the option
        'ffill' as 'forward-filling'
        See Pandas page for details:
        https://pandas.pydata.org/pandas-docs/stable/generated/pandas
        .DataFrame.asfreq.html

    """

    if datastore is None:
        datastore = {}

    # Switching datei and datef if adj
    ddi = min(datei, datef)
    ddf = max(datei, datef)

    # Deals with fluxes
    # Not specified emissions are set to zero by LMDZ
    if input_type == "fluxes":
        make_fluxes(self, datastore, ddi, ddf, runsubdir, mode)

    # Deals with 3D chemical production and loss
    elif input_type in ["prodloss3d", "prescrconcs"]:
        make_chemfields(self, datastore, input_type, ddi, ddf, runsubdir, mode)

    # Deals with initial conditions
    elif input_type == "inicond":
        make_inicond(self, datastore, datei, datef, runsubdir, mode)

    # Deals with observations
    elif input_type == "obs":
        make_obs(self, datastore, runsubdir, mode)

    else:
        self.make_input(datastore, input_type, datei, datef, runsubdir, mode)

    return
