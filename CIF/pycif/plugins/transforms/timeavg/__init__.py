import numpy as np
import pandas as pd

requirements = {"model": {"any": True, "empty": False}}


def obsvect2native(
    transf,
    y0,
    mapper,
    mod_input,
    ddi,
    ddf,
    mode,
    runsubdir,
    workdir,
    trans_mode,
    **kwargs
):
    # Increments are scaled according to dtstep in the adjoint
    # This assumes that obs values are averages from sub-step values
    # TODO: more complex operations might be considered in the future:
    #  e.g., interpolation, gradients, etc.
    #  These should anyway be dealt with
    #  with another transform and not here
    if mode == "adj":
        y0.loc[:, "obs_incr"] /= y0["dtstep"]

    # Observations overlapping two simulation sub-periods are dealt with
    model = transf.model

    subsimu_dates = model.subsimu_dates
    tstep_dates = model.tstep_dates
    tstep_all = model.tstep_all

    # Saving global tstep and dtstep to the transform for later use
    transf.global_dtsteps = y0.loc[:, ["tstep", "dtstep"]]

    # Cropping observations starting before the sub-simulation
    mask = y0.index < ddi
    y0.loc[mask, "tstep"] = 0
    y0.loc[mask, "dtstep"] -= (
        np.argmax(tstep_all == ddi) - y0.loc[mask, "tstep_glo"]
    )

    # Cropping observations ending after the sub-simulation
    if y0.size > 0:
        mask = y0.index + pd.to_timedelta(y0["duration"], unit="h") > ddf
        y0.loc[mask, "dtstep"] = (
            np.argmax(tstep_all == ddf) - y0.loc[mask, "tstep_glo"]
        )

        # For observations from several sub-periods back in time,
        # cropping to the full sub-period extend
        y0.loc[mask, "dtstep"] = np.minimum(
            y0.loc[mask, "dtstep"], len(tstep_dates[ddi]) - 1
        )

    # Change type to integer
    y0.loc[~np.isnan(y0["tstep"]), "tstep"] = y0["tstep"][
        ~np.isnan(y0["tstep"])
    ].astype(int)
    y0.loc[~np.isnan(y0["dtstep"]), "dtstep"] = y0["dtstep"][
        ~np.isnan(y0["dtstep"])
    ].astype(int)

    return y0


def native2obsvect(
    transf,
    y0,
    mapper,
    mod_input,
    ddi,
    ddf,
    mode,
    runsubdir,
    workdir,
    trans_mode,
    **kwargs
):
    # TODO: add in a general way extra columns
    columns = ["sim", "sim_tl", "pressure", "dp", "airm", "hlay"]
    for col in columns:
        if col in y0:
            y0.loc[:, col] /= transf.global_dtsteps.loc[:, "dtstep"]

    return y0


def ini_mapper(
    transf, mapper,
):
    """Initialize mapper for time cropping.
    Does not need to change the mapper as cropping time keep similar
    names and dimensions"""

    return mapper
