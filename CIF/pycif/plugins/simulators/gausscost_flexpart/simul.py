from __future__ import absolute_import
import numpy as np
import pandas as pd

import pycif.utils.check as check
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
    reload_results = getattr(self, 'reload_from_previous', False)
    
    # Get the observation operator from extra arguments
    if not hasattr(self, 'obsoperator'):
        raise Exception("Observation operator is missing to compute the "
                        "simulator. Please check your setup files")
    
    obsoper = self.obsoperator
    controlvect = self.controlvect
    
    # Saving chi to the control vector for later
    controlvect.chi = chi
    
    # Control space contribution of the cost function
    j_b = 0.5 * np.dot(chi, chi)
    
    # Computes forward simulations
    controlvect.x = controlvect.sqrtbprod(chi, **kwargs)

    
    obsvect = obsoper.obsoper(controlvect, 'fwd',
                              datei=self.datei, datef=self.datef,
                              workdir=self.workdir, run_id=run_id,
                              reload_results=reload_results,
                              **kwargs)
    
    # Computes the observation term of the cost function
    # At the moment, this is valid for diagonal observation matrices only
    # TODO: Include these lines into rinvprod.py, with option for
    # non-diagonal matrices, eventually
    departures = obsvect.datastore['sim'] - obsvect.datastore['obs']
    j_o = 0.5 * (departures ** 2 / (obsvect.datastore['obserror']**2 + obsvect.datastore['obs_bkgerr']** 2) ).sum()
    
    zcost = j_b + j_o
    
    towrite = (
        "In Simulator {}:\n"
        "    Jb = {:.4E}\n"
        "    Jo = {:.4E}\n"
        "    Total = {:.4E}") \
        .format(run_id, j_b, j_o, zcost)
    check.verbose(towrite)
    
    # Saves cost function value
    with open('{}/cost.txt'.format(workdir), 'a') as f:
        f.write('{},{:.4E},{:.4E},{:.4E}\n'.format(run_id, j_b, j_o, zcost))
    
    # Return only the cost function if grad = False
    if not grad:
        return zcost

    controlvect.dx = obsvect.dx

    # Observational part of the gradient
    zgrad_o = controlvect.sqrtbprod_ad(controlvect.dx, **kwargs)

    # Control part of the gradient
    zgrad_b = chi
    
    zgrad = zgrad_b + zgrad_o
    
    # Verbose the norms
    znorm_grad_b = np.dot(zgrad_b, zgrad_b) ** 0.5
    znorm_grad_o = np.dot(zgrad_o, zgrad_o) ** 0.5
    check.verbose("In Simulator {}:\n"
                  "    grad(Jb) = {:.4E}\n"
                  "    grad(Jo) = {:.4E}\n"
                  .format(run_id, znorm_grad_b, znorm_grad_o))
    return zcost, zgrad
