from __future__ import division

import copy

import numpy as np

from pycif.utils.check import info


def execute(self, **kwargs):
    """Runs the model in forward mode

    Args:
        self (Plugin): the mode Plugin with all set-up arguments

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

    # Increments in x
    increments = getattr(self, "increments", 0.0)
    incrmode = getattr(self, "incrmode", "cst")
    testspace = getattr(self, "testspace", "control")

    if testspace == "control":
        if incrmode == "cst":
            statevect.dx = increments * statevect.std

        elif incrmode == "rand":
            statevect.dx = (
                np.random.normal(0, increments, statevect.dim) * statevect.std
            )

    elif testspace == "chi":
        if incrmode == "cst":
            statevect.chi[:] = increments

        elif incrmode == "rand":
            statevect.chi = np.random.normal(0, increments, statevect.chi_dim)

        statevect.dx = (
            statevect.sqrtbprod(statevect.chi, **kwargs) - statevect.xb
        )

    if testspace not in ["control", "chi"] or incrmode not in ["cst", "rand"]:
        info(
            "The increment mode you specified can't be parsed: {}".format(
                incrmode
            )
        )
        info(
            "Please check the definition of the running mode "
            "in your Yaml file"
        )
        raise Exception

    statevect.x = copy.deepcopy(statevect.xb)
    statevect.xb += statevect.dx

    # Some verbose
    info("Computing the test of the adjoint")

    # Get the machine accuracy
    accuracy = np.finfo(np.float64).eps

    # Running the tangent linear code of the model
    obsvect = obsoper.obsoper(
        statevect,
        "tl",
        datei=datei,
        datef=datef,
        workdir=workdir,
        reload_results=getattr(self, "reload_results", False),
        **kwargs
    )

    # Computing < H.dx, H.dx >
    scaleprod1 = obsvect.datastore["sim_tl"].pow(2.0).sum()

    # Putting increments in the observation vector
    obsvect.datastore["obs_incr"] = obsvect.datastore["sim_tl"]
    obsvect.datastore.loc[
        np.isnan(obsvect.datastore["obs_incr"]), "obs_incr"
    ] = 0

    # Running the observation operator
    statevect = obsoper.obsoper(
        obsvect,
        "adj",
        datei=datei,
        datef=datef,
        workdir=workdir,
        reload_results=getattr(self, "reload_results", False),
        **kwargs
    )

    # Computing < dx, H*(H.dx) >
    if testspace == "control":
        scaleprod2 = ((statevect.xb - statevect.x) * statevect.dx).sum()

    elif testspace == "chi":
        scaleprod2 = (
            statevect.sqrtbprod_ad(statevect.dx, **kwargs) * statevect.chi
        ).sum()

    else:
        scaleprod2 = np.nan

    # Final verbose
    info("Machine accuracy: {}".format(accuracy))
    info("< H.dx, H.dx >   = {:.17E}".format(scaleprod1))
    info("< dx, H*(H.dx) > = {:.17E}".format(scaleprod2))
    info(
        "The difference is {:.1E} times the machine accuracy".format(
            np.floor(np.abs(scaleprod2 / scaleprod1 - 1) / accuracy)
        )
    )
    
    # Return test
    return np.floor(np.abs(scaleprod2 / scaleprod1 - 1) / accuracy)
    