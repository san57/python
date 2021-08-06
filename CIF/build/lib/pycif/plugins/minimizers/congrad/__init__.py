"""
Python version of the congrad minimization algorithm

Mike Fisher (ECMWF), April 2002
Frederic Chevallier (LSCE), April 2004, for the Python adaptation
"""

from types import MethodType

from .check import check_options
from .congrad import congrad
from .minimize import minimize

requirements = {
    "simulator": {
        "any": True,
        "empty": True,
        "name": "gausscost",
        "version": "std",
    },
}


def ini_data(plugin, **kwargs):
    """Initializes congrad.

    Args:
        plugin (Plugin): options for the variational inversion

    Returns:
        updated plugin

    """

    # Function to check the consistency of options and arguments
    plugin.check_options = MethodType(check_options, plugin)

    # M1QN3 itself
    plugin.congrad = MethodType(congrad, plugin)

    return plugin
