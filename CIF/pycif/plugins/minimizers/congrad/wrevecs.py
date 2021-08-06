from __future__ import division

import numpy as np
from builtins import range
from builtins import str

from pycif.utils.check import info


def wrevecs(
    kdprob, knits, kmaxit, kverbose, pevals, pbnds, pblim, pv, plancv, logfile
):
    """
  Python version of the congrad minimization algorithm

  Mike Fisher (ECMWF), April 2002
  Frederic Chevallier (LSCE), April 2004, for the Python adaptation
  verbose from the Fortran subroutine:

     WREVECS - Called from CONGRAD. Computes and saves eigenvectors.

     Purpose.
     --

     Interface.
     -
        CALL WREVECS (kdprob,knits,kmaxit,pevals,pbnds,pblim,pv,plancv,&
                     &kngood,ptheta)

        Explicit arguments:

        Inputs: kdprob  -- dimension of the problem.
                knits   -- Number of Lanczos vectors.
                kmaxit  -- First dimension of 'pv'.
                kulout  -- unit number for information messages (e.g.
                standard output)
                kverbose -- verbosity (0 => no messages, 1 => taciturn,
                2 => verbose)
                pevals  -- Approximate eigenvalues (Ritz values).
                pbnds   -- Error bounds on eigenvectors.
                pblim   -- Error bound below which eigenvector is 'good'.
                pv      -- Eigenvector matrix of tridiagonal problem.
                plancv  -- Lanczos vectors
        Outputs: kngood -- Number of eigenvectors writen.
                 ptheta -- Eigenvalues .
                 pvcglev -- Eigenvectors
                 prcglev -- Eigenvalues again (originally, this was a file)

     Externals.
     -

     Reference.
     -
         None yet!

     Author.
     -
         Mike Fisher  *ECMWF*

     modifications.
     --
         Original 20/5/94
         M. Fisher 01-01-96 : Optionally save vectors in memory
         M. Fisher 26-11-96 : USE YOM_DISTRIBUTED_VECTORS
         M. Fisher 31-03-99 : Removed the last vestiges of the MIO package
         F. Chevallier 04/05: Translate into python

    """

    ptheta = np.zeros(kmaxit)
    prcglev = np.zeros(kmaxit)
    pvcglev = np.zeros((kmaxit, kdprob))

    # Count and collect converged eigenvalues

    kngood = 0
    for jk in range(knits, -1, -1):
        # WARNING : test on convergence replaced - F Chevallier 20 March 2007 -
        #   if pbnds[jk] <= pblim:
        if pevals[jk] >= 1.0:
            ptheta[kngood] = pevals[jk]
            kngood += 1

    if kverbose > 0:
        info(" ", logfile)
        info(
            "Calculating eigenvectors for the following approximate "
            "eigenvalues:",
            logfile,
        )
        info(str(ptheta[0:kngood]), logfile)
        info(" ", logfile)

    prcglev[0:kngood] = ptheta[0:kngood]

    # Calculate converged eigenvectors

    ij = -1
    for jm in range(knits, -1, -1):

        # WARNING : test on convergence replaced - F Chevallier 20 March 2007 -
        #   if pbnds[jm] <= pblim:
        if pevals[jm] >= 1.0:
            ij += 1
            pvcglev[ij, :] = 0.0
            for jk in range(knits + 1):
                pvcglev[ij, :] += plancv[jk, :] * pv[jk, jm]

        for jk in range(ij):
            dla = np.dot(pvcglev[jk, :], pvcglev[ij, :])
            pvcglev[ij, :] -= dla * pvcglev[jk, :]

        dla = np.dot(pvcglev[ij, :], pvcglev[ij, :])
        pvcglev[ij, :] = pvcglev[ij, :] / np.sqrt(dla)

    return ptheta, pvcglev, prcglev, kngood
