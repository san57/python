from __future__ import absolute_import

from .params import make_rundef, make_species, make_totinput


def make_input(
    self,
    mod_input,
    di,
    df,
    mode,
    runsubdir,
    workdir,
    lmdz_yearref=1979,
    **kwargs
):
    """Prepares inputs for the current LMDZ simulation. It includes:
        - meteo files (defstoke.nc, fluxstoke.nc fluxstokev.nc and phystoke.nc)
        - def files (run.def, totinput, ACTIVE_SPECIES)
        - initial conditions (start.nc)
        - fields for SACS (inca.nc, incaMCF.nc)

    Args:
        - self (Plugin): LMDZ plugin
        - mod_input (str): one of: 'meteo', 'def', 'inicond', 'chem_fields'
        - di (datetime): beginning of the simulation window (along the fwd
        or bckwd axis)
        - df (datetime): end of simulation
        - workdir (str): path to the pycif simulation
        - runsubdir (str): path to the LMDZ sub-simulation
        - mode (str): running mode; fwd or adj
        - lmdz_yearref (int): reference year for LMDZ dates. Default is 1979

    """

    datei = min(di, df)

    return
