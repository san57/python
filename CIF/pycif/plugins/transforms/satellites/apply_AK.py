import numpy as np


def apply_ak(sim_ak, dpavgs, aks, nbformula, qa0lus, chosenlevel):
    """Apply the corresponding AK values"""
    if nbformula == 1:
        y0 = (sim_ak * dpavgs * aks).sum(axis=0) / (dpavgs * aks).sum(axis=0)

    # OMIQA4ECV
    elif nbformula == 2:
        y0 = (sim_ak * aks).sum(axis=0)

    # MOPITT
    elif nbformula == 3:
        y0 = np.log10(qa0lus[chosenlevel]) + (
            aks * (np.log10(sim_ak) - np.log10(qa0lus))
        ).sum(axis=0)
        y0 = np.power(10, y0)

    else:
        raise Exception("Don't know formula")

    return y0


def apply_ak_tl(
    sim_ak_tl, dpavgs, aks, nbformula, qa0lus, chosenlevel, sim_ak
):
    if nbformula == 1 or nbformula == 2:
        y0 = apply_ak(sim_ak_tl, dpavgs, aks, nbformula, qa0lus, chosenlevel)

    elif nbformula == 3:
        tmp = sim_ak_tl / (sim_ak * np.log(10))
        y0 = apply_ak(sim_ak, dpavgs, aks, nbformula, qa0lus, chosenlevel)
        y0 = y0 * np.log(10) * (aks * tmp).sum(axis=0)

    else:
        raise Exception("Don't know formula")

    return y0


def apply_ak_ad(sim_ak, dpavgs, aks, nbformula, qa0lus, chosenlevel, obs_incr):
    if nbformula == 1:
        obs_incr = (
            dpavgs
            * aks.values
            / (dpavgs * aks.values).sum(axis=0)
            * obs_incr.values
        )

    elif nbformula == 2:
        obs_incr = aks.values * obs_incr.values

    elif nbformula == 3:
        y0 = apply_ak(sim_ak, dpavgs, aks, nbformula, qa0lus, chosenlevel)
        y0 *= np.log(10) * obs_incr.values
        y0 *= aks.values

        obs_incr = y0 / (sim_ak * np.log(10))

    else:
        raise Exception("Don't know formula")

    return obs_incr
