from types import MethodType

from .aux import mlis0
from .check import check_options
from .minimize import minimize
from .opt import m1qn3

requirements = {
    "simulator": {
        "any": True,
        "empty": True,
        "name": "gausscost",
        "version": "std",
    },
}


def ini_data(plugin, **kwargs):
    """Initializes M1QN3.

    Args:
        plugin (Plugin): options for the variational inversion

    Returns:
        updated plugin

    """

    # Function to check the consistency of options and arguments
    plugin.check_options = MethodType(check_options, plugin)

    # M1QN3 itself
    plugin.m1qn3 = MethodType(m1qn3, plugin)
    plugin.mlis0 = MethodType(mlis0, plugin)

    return plugin
