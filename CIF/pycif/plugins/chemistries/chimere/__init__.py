from __future__ import absolute_import

import shutil
from distutils.dir_util import copy_tree

from pycif.utils.check import info
from pycif.utils.path import init_dir
from .make_chemistry import create_chemicalscheme
from .read_chemistry import read_chemicalscheme

default_values = {
    #  Max number of reactants/reaction
    "nreactamax": 4,
    #  Number of tabul. temperatures for stoichio.
    "ntemps": 4,
    # Max number of rate constants
    "ntabmax": 22,
    # Max number of tabulated photolysis levels
    "nlevphotmax": 50,
    # Max number of tabulated zenith angles
    "ntabuzenmax": 20,
    # Max number of photolysis reactions
    "nphotmax": 50,
}


def ini_data(self, **kwargs):
    """Initializes the chemistry depending on the model used
    for the inversion.

    Args:
        plugin (ChemistryPlugin): chemistry definition

    Returns:
        Updates on the fly the chemistry
    """

    info("Initializing the Chemistry")

    # Copying the chemical scheme to the working directory
    workdir = self.workdir
    dirchem_ref = "{}/chemical_scheme/{}/".format(workdir, self.schemeid)
    self.dirchem_ref = dirchem_ref

    shutil.rmtree(dirchem_ref, ignore_errors=True)
    init_dir(dirchem_ref)

    # If pre-computed scheme is specified
    if hasattr(self, "dir_precomp"):
        copy_tree(
            "{}/{}/".format(self.dir_precomp, self.schemeid), dirchem_ref
        )

        # Read chemistry
        self.read_chemicalscheme(**kwargs)

    # Otherwise, initialize files from the yaml
    else:
        self.create_chemicalscheme()
