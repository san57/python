import os

from pycif.utils import path
from .simul import simul
from .svd import svd_init

requirements = {
    "statevect": {"any": True, "empty": False},
    "obsvect": {"any": True, "empty": False},
    "obsoperator": {
        "any": True,
        "empty": True,
        "name": "standard",
        "version": "std",
    },
}


def ini_data(plugin, **kwargs):

    # Initializes the directory
    workdir = getattr(plugin, "workdir", "./")
    path.init_dir("{}/simulator".format(workdir))

    os.system("rm -f {}/simulator/cost.txt".format(workdir))
    os.system("rm -f {}/simulator/gradcost.txt".format(workdir))

    # Initializing Singular Value Decomposition if required
    plugin.do_svd = getattr(plugin, "do_svd", False)
    if plugin.do_svd:
        plugin.svd_vectors = svd_init(plugin.obsvect.datastore)
