import copy
import sys

import numpy as np

from pycif.utils.check import info
from .wrevecs import wrevecs


def congrad(self, px, pgrad, planc1, **kwargs):
    """
    """

    # Simulator
    simulator = self.simulator

    # Initializes variables from the minimizer
    zreqrd = self.zreduc
    pevbnd = self.pevbnd
    kvadim = self.kvadim
    kmaxit = self.maxiter
    knevecout = self.knevecout
    kverbose = self.kverbose
    ldsolve = self.ldsolve
    logfile = self.logfile

    if kverbose > 0:
        info("- Lanczos Solver -")

    # Other local parameters
    idprob = kvadim
    itheta1 = 0
    ztheta1 = 0.0

    # Declare local arrays
    zqg0 = np.zeros(kmaxit + 1)
    zbeta = np.zeros(kmaxit + 1)
    zbnds = np.zeros(kmaxit)
    zdelta = np.zeros(kmaxit)
    zv = np.zeros((kmaxit + 1, kmaxit + 1))
    zw = np.zeros(kvadim)
    zetheta = np.zeros(kmaxit)
    zcglwk = np.zeros((kmaxit + 1, kvadim))

    # 'zeta' is an upper bound on the relative error of the gradient.
    zeta = 1e-4

    preduc = 1.0

    # (Note: the factor 0.5 below is because simul defines the cost
    #  function to be J = transpose(chi)*(chi) + 2*Jo, which makes
    #  J''b = 2I, whereas CONGRAD defines J to be half of this value,
    #  so that J''b = I.)
    pgrad /= 2.0
    zgnorm = np.sqrt(np.dot(pgrad, pgrad))

    znorm2l1 = np.dot(planc1, planc1)

    zcglwk[0] = 0.5 * planc1
    zcglwk[0] /= np.sqrt(np.dot(zcglwk[0], zcglwk[0]))

    # Save initial control vector and gradient
    zgrad0 = copy.deepcopy(pgrad)
    zx0 = copy.deepcopy(px)

    zqg0[0] = np.dot(zcglwk[0], zgrad0)

    # Lanczos iteration starts here
    ingood = 0
    iiter = 0
    while 1:  # Lanczos_loop

        # Evaluate the Hessian applied to the latest Lanczos vector
        zw = zx0 + zcglwk[iiter, :]

        zdummy, pgrad = simulator.simul(zw, run_id=iiter, **kwargs)

        pgrad *= 0.5
        pgrad -= zgrad0

        # Calculate zdelta
        zdelta[iiter] = np.dot(zcglwk[iiter, :], pgrad)

        if zdelta[iiter] <= 0.0:
            iiter -= 1
            info("CONGRAD: Hessian is not positive definite")
            info("Stopping after {} iterations".format(iiter))
            sys.exit()

        # Calculate the new Lanczos vector (This is the Lanczos recurrence)
        pgrad -= zdelta[iiter] * zcglwk[iiter]

        if iiter > 0:
            pgrad -= zbeta[iiter] * zcglwk[iiter - 1]

        # Orthonormalize gradient against previous gradients
        for jm in range(iiter, -1, -1):
            dla = np.dot(pgrad, zcglwk[jm, :])
            pgrad -= dla * zcglwk[jm, :]

        zbeta[iiter + 1] = np.sqrt(np.dot(pgrad, pgrad))
        zcglwk[iiter + 1, :] = pgrad / zbeta[iiter + 1]

        zqg0[iiter + 1] = np.dot(zcglwk[iiter + 1, :], zgrad0)

        # Calculate the reduction in the gradient norm
        if ldsolve:
            zwork1 = zdelta[0: iiter + 1]
            zwork2 = zbeta[1: iiter + 1]
            zwork3 = np.zeros(iiter + 1)
            zwork3[0: iiter + 1] = -zqg0[0: iiter + 1]

            zmat = np.zeros((iiter + 1, iiter + 1))
            for i in range(iiter + 1):
                zmat[i, i] = zwork1[i]

            for i in range(iiter):
                zmat[i, i + 1], zmat[i + 1, i] = zwork2[i], zwork2[i]

            zwork3 = np.linalg.solve(zmat, zwork3)

            zw[:] = (
                zgrad0[:]
                + zbeta[iiter + 1] * zcglwk[iiter + 1, :] * zwork3[iiter]
            )

            for j in range(iiter + 1):
                zw[:] = zw[:] - zcglwk[j, :] * zqg0[j]

            preduc_prev = copy.copy(preduc)
            preduc = np.sqrt(np.dot(zw, zw)) / zgnorm
            if kverbose > 1:
                info("Reduction in norm of gradient is: " + str(preduc))

        else:
            preduc = 1.0

        # Determine eigenvalues and eigenvectors of the tri-diagonal problem
        zwork4 = zdelta[0: iiter + 1]
        zwork = zbeta[1: iiter + 1]

        if iiter != 0:
            zmat = np.zeros((iiter + 1, iiter + 1))
            for i in range(iiter + 1):
                zmat[i, i] = zwork4[i]

            for i in range(iiter):
                zmat[i, i + 1], zmat[i + 1, i] = zwork[i], zwork[i]

            zritz, zv = np.linalg.eigh(zmat)

        else:
            zv[0, 0] = 1.0
            zritz = zwork4

        if kverbose > 1:
            info(
                "congrad: ritz values are: {}".format(zritz[0: iiter + 1]),
                logfile,
            )

        # Estimate error bounds
        zbndlm = zeta * zritz[iiter]

        zbnds[0: iiter + 1] = abs(zbeta[iiter + 1] * zv[iiter, 0: iiter + 1])
        if kverbose > 1:
            info(
                "congrad: error bounds are: " + str(zbnds[0: iiter + 1]),
                logfile,
            )

        # Check for exploding or negative Ritz values
        if min(zritz[0: iiter + 1]) < 0.0:
            if kverbose > 0:
                info("congrad: stopping: negative ritz value")
            preduc = preduc_prev
            iiter = iiter - 1
            zwork4 = zdelta[0: iiter + 1]
            zwork = zbeta[1: iiter + 1]

            if iiter > 0:
                zmat = np.zeros((iiter + 1, iiter + 1))
                for i in range(iiter + 1):
                    zmat[i, i] = zwork4[i]

                for i in range(iiter):
                    zmat[i, i + 1], zmat[i + 1, i] = zwork[i], zwork[i]

                zritz, zv = np.linalg.eigh(zmat)

            else:
                zv[0, 0] = 1.0

            zbnds[0: iiter + 1] = abs(
                zbeta[iiter + 1] * zv[iiter, 0: iiter + 1]
            )
            break  # Lanczos loop

        if ingood > 0:
            if zritz[itheta1] > 1.01 * ztheta1 and kverbose > 0:
                info("CONGRAD: warning -- ritz values exploding", logfile)
                info("leading ritz value=" + str(zritz[itheta1]), logfile)
                info("leading converged eigenvalue=" + str(ztheta1), logfile)

        # Count the converged eigenvectors
        ing = 0
        for jm in range(iiter + 1):
            if zbnds[jm] <= zbndlm:
                ing = ing + 1
                if kverbose > 1:
                    info(
                        "leading converged eigenvalue=" + str(zritz[jm]),
                        logfile,
                    )

        # Deal with newly converged eigenvector and recompute eigenvectors
        ingood = ing

        # Save leading converged eigenvalue for explosion test
        if ingood > 0:
            for jm in range(iiter, -1, -1):
                if zbnds[jm] <= zbndlm:
                    ztheta1 = zritz[jm]
                    itheta1 = jm
                    break

        if kverbose > 1:
            info("congrad: End of iteration: " + str(iiter + 1), logfile)

        if iiter >= (kmaxit - 1) or preduc <= zreqrd:
            break  # Lanczos loop

        iiter = iiter + 1
        if ingood > 0:
            itheta1 = itheta1 + 1

    # End of Lanczos iteration
    if kverbose > 0:
        info("Summary of Lanczos iteration:")
        info("   Number of iterations performed: " + str(iiter + 1), logfile)
        if ldsolve:
            info(
                "   Maximum allowed number of iterations: " + str(kmaxit),
                logfile,
            )
            info(
                "   Required reduction in norm of gradient: " + str(zreqrd),
                logfile,
            )
            info(
                "   Achieved reduction in norm of gradient: " + str(preduc),
                logfile,
            )
    if preduc > zreqrd:
        info(
            "   *** Failed to achieve required reduction in gradient ***",
            logfile,
        )

    # Calculate sufficiently converged eigenvectors of the
    # preconditioned Hessian
    if ldsolve:
        zbndlm = pevbnd
        zbnds[0: iiter + 1] = zbnds[0: iiter + 1] / zritz[0: iiter + 1]

        ztheta, zvcglev, zrcglev, ingood = wrevecs(
            idprob, iiter, kmaxit, kverbose, zritz, zbnds, zbndlm, zv, zcglwk,
            logfile
        )

        if ingood == 0 and kverbose > 0:
            info("Warning: congrad found no eigenpairs")

        if kverbose > 0:
            info(
                "number of eigenpairs converged to requested accuracy="
                + str(ingood)
            )

    ###########
    # pevecs = np.zeros((knevecout,idprob))
    # return px, pgrad, preduc, pevecs, iter
    ###########
    if not ldsolve:
        # Determine bounds on the quadratic form: v'*inv(J'')*v
        # (where v=zgrad0) using theorem 5.3 of Golub and Meurant 1994.
        zbjm1 = 1.0 / zdelta[0]
        zdjm1 = zdelta[0]
        zcjm1 = 1.0
        zdhatjm1 = zdelta[0] - 1.0
        zdbarjm1 = zdelta[0] - ztheta1
        pgolubl_inv = np.zeros(iiter + 1)
        pgolubu_inv = np.zeros(iiter + 1)

        for j in range(1, iiter + 1):
            # zbj is the Gauss rule lower bound:
            zbj = zbjm1 + zbeta[j] * zbeta[j] * zcjm1 * zcjm1 / (
                zdjm1 * (zdelta[j] * zdjm1 - zbeta[j] * zbeta[j])
            )

            zdj = zdelta[j] - zbeta[j] * zbeta[j] / zdjm1
            zcj = zcjm1 * zbeta[j] / zdjm1
            zdhatj = zdelta[j] - 1.0 - zbeta[j] * zbeta[j] / zdhatjm1
            zdbarj = zdelta[j] - ztheta1 - zbeta[j] * zbeta[j] / zdbarjm1
            zohatj = 1.0 + zbeta[j + 1] * zbeta[j + 1] / zdhatj
            zobarj = ztheta1 + zbeta[j + 1] * zbeta[j + 1] / zdbarj

            # zbhatj is the Gauss-Radau rule upper bound:
            # zbbarj is the Gauss-Radau rule lower bound:
            zbhatj = zbj + zbeta[j + 1] * zbeta[j + 1] * zcj * zcj / (
                zdj * (zohatj * zdj - zbeta[j + 1] * zbeta[j + 1])
            )
            zbbarj = zbj + zbeta[j + 1] * zbeta[j + 1] * zcj * zcj / (
                zdj * (zobarj * zdj - zbeta[j + 1] * zbeta[j + 1])
            )

            # zbuj is the Gauss-Lobatto rule upper bound:
            zouj = (
                zdhatj
                * zdbarj
                * (ztheta1 / zdhatj - 1.0 / zdbarj)
                / (zdbarj - zdhatj)
            )
            zgamuj2 = zdhatj * zdbarj * (ztheta1 - 1.0) / (zdbarj - zdhatj)

            zbuj = zbj + zgamuj2 * zcj * zcj / (zdj * (zouj * zdj - zgamuj2))

            zbjm1 = zbj
            zdjm1 = zdj
            zcjm1 = zcj
            zdhatjm1 = zdhatj
            zdbarjm1 = zdbarj
            pgolubl_inv[j] = max(zbbarj, zbj)
            pgolubu_inv[j] = min(zbhatj, zbuj)

        # Calculate bounds on (planc1)' log[(J'')^-1] planc1
        # note that the former computation is an alternate version of
        #   the one just above
        # The lines for the inv have been commented in this alternate algorithm
        pgolubu_inv = np.zeros(iiter + 1)
        pgolubu_log = np.zeros(iiter + 1)
        pgolubu_inv[0] = zv[0, 0] * zv[0, 0] / zritz[0]
        pgolubu_log[0] = zv[0, 0] * zv[0, 0] * np.log(zritz[0])
        for j in range(1, iiter + 1):
            pgolubu_inv[j] = (
                pgolubu_inv[j - 1] + zv[0, j] * zv[0, j] / zritz[j]
            )
            pgolubu_log[j] = pgolubu_log[j - 1] + zv[0, j] * zv[0, j] * np.log(
                zritz[j]
            )

        za = 1.0  # za must be less than the smallest eigenvalue
        #   za = 1.e-8

        zwork1 = zdelta[0: iiter + 1] - za
        zwork2 = zbeta[1: iiter + 1]
        zwork3 = np.zeros(iiter + 1)
        zwork3[iiter] = zbeta[iiter + 1] * zbeta[iiter + 1]

        zmat = np.zeros((iiter + 1, iiter + 1))
        for i in range(iiter + 1):
            zmat[i, i] = zwork1[i]

        for i in range(iiter):
            zmat[i, i + 1], zmat[i + 1, i] = zwork2[i], zwork2[i]

        zwork3 = np.linalg.solve(zmat, zwork3)

        zwork4 = np.zeros(iiter + 2)
        zwork4[0: iiter + 1] = zdelta[0: iiter + 1]
        zwork4[iiter + 1] = za + zwork3[iiter]
        zwork = np.zeros(iiter + 1)
        zwork[0: iiter + 1] = zbeta[1: iiter + 2]
        iiter += 1
        zmat = np.zeros((iiter + 1, iiter + 1))
        for i in range(iiter + 1):
            zmat[i, i] = zwork4[i]

        for i in range(iiter):
            zmat[i, i + 1], zmat[i + 1, i] = zwork[i], zwork[i]

        zritz, zv = np.linalg.eigh(zmat)
        iiter -= 1

        pgolubl_inv = np.zeros(iiter + 1)
        pgolubl_log = np.zeros(iiter + 1)
        pgolubl_inv[0] = zv[0, 0] * zv[0, 0] / zwork4[0]
        pgolubl_log[0] = zv[0, 0] * zv[0, 0] * np.log(zwork4[0])
        for j in range(1, iiter + 1):
            pgolubl_inv[j] = (
                pgolubl_inv[j - 1] + zv[0, j] * zv[0, j] / zwork4[j]
            )
            pgolubl_log[j] = pgolubl_log[j - 1] + zv[0, j] * zv[0, j] * np.log(
                zwork4[j]
            )

        # Convert bounds to degrees of freedom for signal and entropy
        pgolubu_inv = znorm2l1 * (1.0 - pgolubu_inv)
        pgolubl_inv = znorm2l1 * (1.0 - pgolubl_inv)
        pgolubu_log = 0.5 * znorm2l1 * pgolubu_log / np.log(2.0)
        pgolubl_log = 0.5 * znorm2l1 * pgolubl_log / np.log(2.0)

    # Calculate the solution vector and gradient
    if ldsolve:

        zwork1 = zdelta[0: iiter + 1]
        zwork2 = zbeta[1: iiter + 1]
        zwork3 = np.zeros(iiter + 1)
        zwork3[0: iiter + 1] = -zqg0[0: iiter + 1]

        zmat = np.zeros((iiter + 1, iiter + 1))
        for i in range(iiter + 1):
            zmat[i, i] = zwork1[i]

        for i in range(iiter):
            zmat[i, i + 1], zmat[i + 1, i] = zwork2[i], zwork2[i]

        zwork3 = np.linalg.solve(zmat, zwork3)

        for j in range(iiter + 1):
            px[:] = px[:] + zcglwk[j, :] * zwork3[j]
            pgrad -= zcglwk[j, :] * zqg0[j]

    else:
        pgrad = copy.deepcopy(zgrad0)

    # Determine eigenpairs of the un-preconditioned Hessian in 'chi'
    # space, transform to model variable space and write out.
    pevecs = np.zeros((knevecout, idprob))
    if ldsolve:
        intot = ingood
        if intot > 0:
            # subroutine xformev removed - F Chevallier 4 April 2012
            knevecout = min(intot, knevecout)
            for jk in range(knevecout):
                pevecs[jk, :] = zvcglev[jk, :] * np.sqrt(
                    1.0 - 1.0 / zrcglev[jk]
                )

    # Transform control variable and gradient back to space with scalar
    # product SCAAS
    pgrad *= 2.0

    # Return the number of iterations actually performed
    iiter = iiter + 1

    if ldsolve:
        return px, pgrad, preduc, pevecs, iiter
    else:
        return pgolubl_inv, pgolubu_inv, pgolubl_log, pgolubu_log, iiter
