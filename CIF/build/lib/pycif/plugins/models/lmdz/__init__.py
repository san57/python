from __future__ import absolute_import

import os
import shutil

from pycif.utils import path
from pycif.utils.check import info
from pycif.utils.classes.setup import Setup
from .flushrun import flushrun
from .ini_mapper import ini_mapper
from .ini_periods import ini_periods
from .inputs import make_input
from .inputs.fluxes import make_fluxes
from .inputs.meteo import make_meteo
from .native2inputs import native2inputs
from .outputs2native import outputs2native
from .run import run

requirements = {
    "domain": {"name": "LMDZ", "version": "std", "empty": False},
    "fluxes": {
        "name": "LMDZ",
        "version": "sflx",
        "empty": True,
        "newplg": True,
    },
    "chemistry": {"name": "CHIMERE", "version": "gasJtab", "empty": False},
    "emis_species": {
        "name": "LMDZ",
        "version": "bin",
        "type": "fluxes",
        "empty": True,
        "newplg": True,
    },
    "meteo": {
        "name": "LMDZ",
        "version": "mass-fluxes",
        "empty": True,
        "any": False,
        "subplug": True,
        "preftree": "statevect/components",
    },
    "inicond": {
        "name": "LMDZ",
        "version": "ic",
        "empty": True,
        "any": False,
        "type": "fields",
        "newplg": True,
    },
    "prescrconcs": {
        "name": "LMDZ",
        "version": "prescrconcs",
        "empty": True,
        "any": False,
        "type": "fields",
        "newplg": True,
    },
}

# Required inputs for running a LMDz simulations
required_inputs = [
    "fluxes",
    "meteo",
    "inicond",
    "def",
    "chem_fields",
    "prescrconcs",
    "prodloss3d",
    "traj",
    "exe",
]


def ini_data(plugin, **kwargs):
    """Initializes LMDZ

    Args:
        plugin (Plugin): the model plugin to initialize
        **kwargs (dictionary): possible extra parameters

    Returns:
        loaded plugin and directory with executable

    """

    info("Initializing the model")

    workdir = getattr(plugin, "workdir", "./")

    # Cleaning the model working directory
    shutil.rmtree("{}/model/".format(workdir), ignore_errors=True)

    # Initializes the directory
    path.init_dir("{}/model".format(workdir))

    # copying the executable
    target = "{}/model/".format(workdir) + os.path.basename(plugin.fileexec)
    source = plugin.fileexec
    shutil.copy(source, target)

    # copying the definition file
    target = "{}/model/run.def".format(workdir)
    source = plugin.filedef
    shutil.copy(source, target)

    # LMDZ has a fixed integration time step
    plugin.tstep = 0

    # Initializes default values
    # Period of sub-simulations: default = 1 month
    if not hasattr(plugin, "periods"):
        plugin.periods = "1MS"

    # Convection scheme: default = TK = Tiedke
    if not hasattr(plugin, "conv_scheme"):
        plugin.conv_scheme = "TK"

    # Loading input fluxes if specified, otherwise, loads default inputs
    for spec in plugin.chemistry.emis_species.attributes:
        tracer = getattr(plugin.chemistry.emis_species, spec)

        if hasattr(tracer, "provider") and hasattr(tracer, "format"):
            name = tracer.provider
            version = tracer.format

        else:
            name = "LMDZ"
            version = "sflx"

        tracer = Setup.load_registered(
            name, version, "fluxes", plg_orig=tracer
        )
        Setup.load_setup(
            Setup.from_dict({"fluxes": tracer}), "fluxes", level=1, **kwargs
        )

        setattr(plugin.chemistry.emis_species, spec, tracer)

    return plugin
