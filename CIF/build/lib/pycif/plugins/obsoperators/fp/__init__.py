from types import MethodType
from .obsoper import obsoper

from pycif.utils import path
from pycif.utils.check import verbose

requirements = {'model': {'any': True, 'empty': False},
                'obsvect': {'any': True, 'empty': True,
                            'name': 'standard', 'version': 'std'},
                'controlvect': {'any': True, 'empty': True,
                                'name': 'standard', 'version': 'std'},
                }


def ini_data(plugin, **kwargs):
    """Initializes the observation operator

    Args:
        plugin (dict): dictionary defining the plugin
        **kwargs (dictionary): possible extra parameters

    """
    
    verbose("Initializing the observation operator")
    
    workdir = plugin.workdir
    
    # Initializes the directory
    path.init_dir('{}/obsoperator'.format(workdir))
    
    return plugin
