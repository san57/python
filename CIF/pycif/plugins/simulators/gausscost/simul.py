from __future__ import absolute_import

import numpy as np

from pycif.utils.check import info
from .svd import svd_cost


def simul(self, chi, grad=True, run_id=-1, **kwargs):
    """Computes the cost function J (and its gradient) based on the Gaussian
    formulation of the inversion framework:
    J(x) = chi^T chi + (Hx-y)^T R^-1 (Hx-y)
         = j_b         + j_r

    gradJ(x) = 2 * chi + 2 * H^T R^(-1) (Hx-y)


    Args:
        chi (np.array): a flat vector defining the current state of the control
                        vector
        grad (bool, optional): if True, returns both the function value and
                               its gradient
        run_id (int): ID for the current run (determines the folder names)

    Returns:
        J(x), gradJ(x)
    """

    # Various variables
    datei = self.datei
    datef = self.datef
    workdir = self.workdir
    reload_results = getattr(self, "reload_from_previous", False)

    # Get the observation operator from extra arguments
    if not hasattr(self, "obsoperator"):
        raise Exception(
            "Observation operator is missing to compute the "
            "simulator. Please check your setup files"
        )

    obsoper = self.obsoperator
    statevect = self.statevect

    # Saving chi to the control vector for later
    statevect.chi = chi

    # Control space contribution of the cost function
    j_b = 0.5 * np.dot(chi, chi)

    # Computes forward simulations
    statevect.x = statevect.sqrtbprod(chi, **kwargs)

    obsvect = obsoper.obsoper(
        statevect,
        "fwd",
        datei=self.datei,
        datef=self.datef,
        workdir=self.workdir,
        run_id=run_id,
        reload_results=reload_results,
        **kwargs
    )

    # Computes the observation term of the cost function
    departures = obsvect.datastore["sim"] - obsvect.datastore["obs"]
    j_o = 0.5 * (departures * obsvect.rinvprod(departures)).sum()

    # Computes SVD-based cost function
    if self.do_svd:
        j_o = svd_cost(self, obsvect.datastore, j_o)

    zcost = j_b + j_o

    towrite = (
        "In Simulator {}:\n" "    Jb = {}\n" "    Jo = {}\n" "    Total = {}"
    ).format(run_id, j_b, j_o, zcost)
    info(towrite)

    # Saves cost function value
    with open("{}/simulator/cost.txt".format(workdir), "a") as f:
        f.write("{},{},{},{}\n".format(run_id, j_b, j_o, zcost))

    # Return only the cost function if grad = False
    if not grad:
        return zcost

    # Runs the adjoint to get the gradients
    if not self.do_svd:
        obsvect.datastore["obs_incr"] = obsvect.rinvprod(departures)

    statevect = obsoper.obsoper(
        obsvect,
        "adj",
        datei=datei,
        datef=datef,
        workdir=workdir,
        run_id=run_id,
        reload_results=reload_results,
        **kwargs
    )

    # Observational part of the gradient
    zgrad_b = chi

    # Control part of the gradient
    zgrad_o = statevect.sqrtbprod_ad(statevect.dx, **kwargs)

    zgrad = zgrad_b + zgrad_o

    # Verbose the norms
    znorm_grad_b = np.dot(zgrad_b, zgrad_b) ** 0.5
    znorm_grad_o = np.dot(zgrad_o, zgrad_o) ** 0.5
    info(
        "In Simulator {}:\n"
        "    grad(Jb) = {}\n"
        "    grad(Jo) = {}\n".format(run_id, znorm_grad_b, znorm_grad_o)
    )

    # Saves cost function gradient
    with open("{}/simulator/gradcost.txt".format(workdir), "a") as f:
        f.write("{},{},{}\n".format(run_id, znorm_grad_b, znorm_grad_o))

    return zcost, zgrad
