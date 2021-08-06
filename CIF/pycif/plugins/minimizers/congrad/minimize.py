import copy

import numpy as np
from builtins import str

from pycif.utils.check import info


def minimize(self, finit, gradinit, chi0, **kwargs):
    # x, f, g, auxil, io, niter, nsim, iz, df1, m=5, dxmin=1.e-20,
    #  epsg=1.e-20, impres=1, mode=0, **kwargs
    """Entry point for CONGRAD algorithm.

    Args:
        finit (float): initial value for the function to minimize
        gradinit (np.array): gradient at the starting point
        chi (np.array): initial state for the unknown to optimize
        simulator (module): simulator module to evaluate the function and
                               its gradient
        minimizer (module): minimizer module, used to define minimizer options

    Returns:
        (np.array, float): a tuple with the optimized vector and the
                           corresponding function maximum

    """

    # Initializing options (and filling missing values with default)
    self = self.check_options(chi0, **kwargs)

    # Running CONGRAD
    lanczvect0 = copy.deepcopy(gradinit)
    xopt, gradopt, preduc, pevecs, iiter = self.congrad(
        chi0, gradinit, lanczvect0, **kwargs
    )

    # Final verbose and output
    towrite = """
        CONGRAD:
            number of iterations: {}
            achieved relative reduction of the gradient: {}
        """.format(
        iiter, preduc
    )

    info(towrite)

    r1 = np.sqrt(np.dot(xopt, xopt))
    r2 = np.sqrt(np.dot(gradopt, gradopt))

    info("norm of x = " + str(r1))
    info("norm of g = " + str(r2))

    return xopt
