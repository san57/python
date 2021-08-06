from __future__ import division

import pandas as pd

from pycif.utils import path
from pycif.utils.check import info
from .compile import compile
from .flushrun import flushrun
from .ini_mapper import ini_mapper
from .ini_periods import ini_periods
from .native2inputs import native2inputs
from .outputs2native import outputs2native
from .run import run

requirements = {
    "domain": {
        "name": "CHIMERE",
        "version": "std",
        "empty": False,
        "any": False,
    },
    "chemistry": {
        "name": "CHIMERE",
        "version": "gasJtab",
        "empty": False,
        "any": False,
    },
    "fluxes": {
        "name": "CHIMERE",
        "version": "AEMISSIONS",
        "empty": True,
        "any": False,
        "subplug": True,
        "preftree": "statevect/components",
    },
    "biofluxes": {
        "name": "CHIMERE",
        "version": "AEMISSIONS",
        "type": "fluxes",
        "empty": True,
        "any": False,
        "subplug": True,
        "preftree": "statevect/components",
        "emis_type": "bio"
    },
    "meteo": {
        "name": "CHIMERE",
        "version": "std",
        "empty": True,
        "any": False,
        "subplug": True,
        "preftree": "statevect/components",
    },
    "latcond": {
        "name": "CHIMERE",
        "version": "icbc",
        "empty": True,
        "any": False,
        "type": "fields",
        "newplg": True,
    },
    "topcond": {
        "name": "CHIMERE",
        "version": "icbc",
        "empty": True,
        "any": False,
        "type": "fields",
        "newplg": True,
    },
    "inicond": {
        "name": "CHIMERE",
        "version": "icbc",
        "empty": True,
        "any": False,
        "type": "fields",
        "newplg": True,
    },
}

default_values = {
    #  Number of spin-up hours
    "ihoursu": 0,
    # Include biogenic emissions
    "optemisb": False,
    # Precision for output variables: double or float
    "dumpnctype": "float",
    # Dump outputs into a NetCDF file
    "dumpncoutput": True,
    # Dump parameters into a NetCDF file
    "dumpncpar": False,
    # ! Save cumulated deposition every ... hours
    "nsavedepos": 4,
    # Use chemistry
    "usechemistry": 0,
    # emissions
    "useemissions": 1,
    # transport and mixing
    "usetransmix": 1,
    # wet deposition
    "usewetdepos": 0,
    # Use dry deposition
    "usedepos": 0,
    # Use dry deposition
    "dryairout": 0,
    # Number of Gauss-Seidel iterations in the TWOSTEP solver.
    # 1 for model testing, 2 for higher accuracy
    "nitgs": 1,
    # Same but during spin-up
    "nitgssu": 1,
    # clipping of small (in absolute value) concentrations
    "useabsclipconc": 0,
    # Clipping value for the TWOSTEP algorithm
    "clipconc": "1d0",
    # Max number of reaction types
    "ntyperate": 50,
    # Number of vegetation types
    "nvegtype": 16,
    # Max number of landuse classes
    "nlduse": "09",
    # Max number of output parameters
    "nparammax": 30,
    # hour of the pulse
    "hpulse": 0,
    # Auto-compile before running
    "auto-recompile": False,
    # Compile mode: DEBUG or PROD
    "compile-mode": "PROD",
    # Clean compile
    "compile-clean": True,
    # Number of levels in BEMISSIONS
    "nlevemis_bio": 1,
}

# Replacement component if not defined in state vector
backup_comps = {"latcond": "background",
                "topcond": "background"}

# Required inputs for running a CHIMERE simulation
required_inputs = [
    "exe",
    "nml",
    "fluxes",
    "biofluxes",
    "meteo",
    "inicond",
    "latcond",
    "topcond",
]


def ini_data(plugin, **kwargs):
    """Initializes CHIMERE

    Args:
        plugin (dict): dictionary defining the plugin
        **kwargs (dictionary): possible extra parameters

    Returns:
        loaded plugin and directory with executable

    """

    info("Initializing the model")

    workdir = getattr(plugin, "workdir", "./")

    # Initializes the directory
    path.init_dir("{}/model".format(workdir))

    # Default values:
    # period: '1D'
    plugin.periods = getattr(plugin, "periods", "1D")

    # Number of hours per period
    plugin.nhours = int(
        pd.to_timedelta(plugin.periods).total_seconds() // 3600
    )
    plugin.nho = "{:.0f}".format(plugin.nhours)

    # Replacing nsaveconcs if not specified
    # Forces the end.nc file to contain concentration every N hours
    # By default, saves only at the end
    if not hasattr(plugin, "nsaveconcs"):
        plugin.nsaveconcs = plugin.nhours

    # Replace name for METEO files
    plugin.meteo.file = plugin.meteo.file.format(nho=plugin.nho)

    # Replace name for AEMISSION files and BEMISSIONS files
    plugin.fluxes.file = plugin.fluxes.file.format(nho=plugin.nho)
    plugin.fluxes.nlevemis = plugin.nlevemis
    
    plugin.biofluxes.file = plugin.biofluxes.file.format(nho=plugin.nho)
    plugin.biofluxes.nlevemis = plugin.nlevemis_bio
    
    # Replace name for BOUN_CONCS files
    plugin.latcond.file = plugin.latcond.file.format(nho=plugin.nho)
    plugin.topcond.file = plugin.topcond.file.format(nho=plugin.nho)

    return plugin
