from pycif.utils import path
from pycif.utils.check import info
from .obsoper import obsoper
from .transforms import init_transform
from .transforms.do_transforms import do_transforms

requirements = {
    "model": {"any": True, "empty": False},
    "obsvect": {
        "any": True,
        "empty": True,
        "name": "standard",
        "version": "std",
    },
    "statevect": {
        "any": True,
        "empty": True,
        "name": "standard",
        "version": "std",
    },
    "platform": {"any": True, "empty": True},
}


def ini_data(plugin, **kwargs):
    """Initializes the observation operator

    Args:
        plugin (dict): dictionary defining the plugin
        **kwargs (dictionary): possible extra parameters

    """

    info("Initializing the observation operator")

    workdir = plugin.workdir

    # Initializes the directory
    path.init_dir("{}/obsoperator".format(workdir))

    # Initializes transforms
    init_transform(plugin, plugin.statevect)
    init_transform(plugin, plugin.obsvect, transform_type="obs")

    # Re-compile model if necessary
    if hasattr(plugin.model, "compile"):
        plugin.model.compile()

    return plugin
