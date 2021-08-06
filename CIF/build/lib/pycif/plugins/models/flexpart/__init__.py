import shutil
from types import MethodType
import pandas as pd
from pycif.utils import path
from pycif.utils.check import verbose
from .ini_periods import ini_periods
#from .inputs import make_input
from .native2inputs import native2inputs
from .outputs2native import outputs2native
from .run import run
from .ini_mapper import ini_mapper

requirements = {'domain': {'name': 'FLEXPART', 'version': 'std',
                           'empty': False, 'any': False},
                'fluxes': {'name': 'FLEXPART', 'version': 'nc',
                           'empty': True, 'any': False}}

# Required inputs for running a FLEXPART simulation
required_inputs = ['fluxes']


def ini_data(plugin, **kwargs):
    """Initializes FLEXPART

    Args:
        plugin (dict): dictionary defining the plugin
        **kwargs (dictionary): possible extra parameters

    Returns:
        loaded plugin and directory with executable

    """

    # Initializes default values
    # Period of sub-simulations: default = 1 month
    if not hasattr(plugin.plugin, 'periods'):
        plugin.periods = '1MS'
    else:
        plugin.periods = plugin.plugin.periods
    
    return plugin

