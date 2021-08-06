from __future__ import division

import numpy as np
from builtins import str

from pycif.utils.check import info


def execute(self, **kwargs):
    """Performs a variational inversion given a minimizer method and a
    simulator (i.e. a function to minimize and its gradient)

    Args:
        self (Plugin): definition of the mode set-up

    """

    # Working directory
    workdir = self.workdir

    # Control vector
    statevect = self.statevect

    # Observation operator
    obsoper = self.obsoperator

    # Simulation window
    datei = self.datei
    datef = self.datef

    # Minimizer
    minimizer = self.minimizer

    # Simulator
    simulator = self.simulator

    # Some verbose
    towrite = """
        Running a variational inversion with the following modules:
            Minimizer: {}
            Simulator: {}
        """.format(
        minimizer.plugin.name, simulator.plugin.name
    )
    info(towrite)

    # Initial run of the simulator as a starting point for M1QN3
    costinit, gradinit = simulator.simul(statevect.chi, run_id=-1, **kwargs)

    zgnorm = np.sqrt(np.sum(gradinit ** 2))
    info("Nb of observations: " + str(len(obsoper.obsvect.datastore)))
    info("Initial cost: " + str(costinit))
    info("Initial gradient norm: " + str(zgnorm))

    # Runs the minimizer
    chiopt = minimizer.minimize(costinit, gradinit, statevect.chi, **kwargs)

    # Last call to the simulator for final diagnostics
    costend, gradend = simulator.simul(chiopt, run_id=-2, **kwargs)

    zgnorm = np.sqrt(np.dot(gradend, gradend)) / zgnorm
    info("Final cost: " + str(costend))
    info("Ratio final/initial gradient norm: " + str(zgnorm))

    # Save results
    statevect.dump(
        "{}/statevect_final.pickle".format(workdir),
        to_netcdf=getattr(self, "save_out_netcdf", False),
        dir_netcdf="{}/statevect/".format(workdir),
    )

    return chiopt
