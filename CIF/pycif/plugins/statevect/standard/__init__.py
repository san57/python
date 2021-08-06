from pycif.utils.path import init_dir
from .control2native import control2native
from .dump import dump, load
from .init_bprod import init_bprod
from .init_components import init_components
from .init_xb import init_xb
from .native2control import native2control
from .sqrtbprod import sqrtbprod, sqrtbprod_ad

# It is necessary to initialize a domain, fluxes and the model itself
requirements = {
    "domain": {"any": True, "empty": False},
    "model": {"any": True, "empty": False},
    "components": {"any": True, "empty": True, "type": "fields"},
    "obsvect": {
        "any": True,
        "empty": True,
        "name": "standard",
        "version": "std",
    },
}


def ini_data(plugin, **kwargs):
    """Initializes the state vector from information in the Yaml file

    Args:
        plugin (pycif.classes.plugins): the plugin to initialize

    """

    # Initializes reference directories if needed
    init_dir("{}/statevect/".format(plugin.workdir))

    # Saves reference directories and file formats if not prescribed
    # Look for directory by order of priority:
    # 1) directly in tracer definition
    # 2) in component definition if specified
    # 3) in model fluxes if any
    # Getting the right emissions
    init_components(plugin)

    # Initializes the control vector
    plugin.init_xb(plugin, **kwargs)

    # TODO: fix general data structure to enable compatibility
    #  between control vector and transforms
    plugin.datastore = {}

    # Initializing the product of chi by B^(1/2), only if components specified
    if hasattr(plugin, "components"):
        plugin.init_bprod(plugin, **kwargs)
