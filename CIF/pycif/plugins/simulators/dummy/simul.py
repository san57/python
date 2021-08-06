import numpy as np

from pycif.utils.check import info


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

    zcost = np.sum((chi - np.arange(len(chi))) ** 2)
    zgrad = 2 * (chi - np.arange(len(chi)))

    # Verbose the norms
    znorm_grad_b = np.dot(zgrad, zgrad) ** 0.5
    info(
        "In Simulator:\n"
        "    grad(Jb) = {}\n"
        "         Jb  = {}\n".format(znorm_grad_b, zcost)
    )

    return zcost, zgrad
