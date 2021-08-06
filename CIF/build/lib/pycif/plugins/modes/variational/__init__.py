from .execute import execute

requirements = {
    "obsvect": {
        "any": True,
        "empty": False,
        "name": "standard",
        "version": "std",
    },
    "statevect": {
        "any": True,
        "empty": False,
        "name": "standard",
        "version": "std",
    },
    "obsoperator": {
        "any": True,
        "empty": True,
        "name": "standard",
        "version": "std",
    },
    "minimizer": {
        "any": True,
        "empty": True,
        "name": "m1qn3",
        "version": "std",
    },
    "simulator": {
        "any": True,
        "empty": True,
        "name": "gausscost",
        "version": "std",
    },
}


def ini_data(plugin, **kwargs):
    """Initializes the variational inversions. It means that all sub-modules
    and plugins required to run the inversion are built at this step.

    Args:
        plugin (pycif.classes.plugins): the plugin to initialize

    Returns:
        dict with all options and sub-modules

    """

    if not hasattr(plugin, "minimizer"):
        raise Exception(
            "A variational inversion was asked with no " "optimization method"
        )

    # Initializing the simulator from the minimizer
    if not hasattr(plugin.minimizer, "simulator"):
        raise Exception("You are trying to run a minimizer without simulator")

    plugin.simulator = plugin.minimizer.simulator

    return plugin
