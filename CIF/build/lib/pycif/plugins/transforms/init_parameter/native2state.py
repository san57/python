import numpy as np

from pycif.utils.check import info
from .utils.dates import dateslice
from .utils.scalemaps import map2scale, vmap2vaggreg


def native2state(
    transform,
    data,
    mapper,
    mod_input,
    di,
    df,
    mode,
    runsubdir,
    workdir,
    trans_mode,
    **kwargs
):
    """Translates real-size data as extracted from the model outputs to
     the control space. This includes mainly temporal and spatial aggregation.
     This routine is used in adjoint mode, thus computes operations on
     increments.

     Args:
         self (Plugin): the control vect
         datastore (dict): the data at the model resolution to be converted
                           to the control space
         di (datetime): starting date of the simulation window
         df (datetime): ending date of the simulation window
         workdir (str): pycif working directory

     """

    ddi = min(di, df)
    datastore = data.datastore

    tracer_ids = mapper["outputs"]

    # If dx is not yet defined, initialize it
    if not hasattr(data, "dx"):
        info("Setting dx to zero in the control vector")
        data.dx = 0.0 * data.x

    # Loop over model sensitivities
    for tracer_id in tracer_ids:
        mod_input = tracer_id[0]
        trcr = tracer_id[1]

        # If this type of input is not considered in the control vector,
        # ignoring the model sensitivity
        component = getattr(getattr(data, "components", None), mod_input, None)
        parameters = getattr(component, "parameters", None)

        if parameters is None:
            info(
                "{} is sensitive to {} but your control vector doesn't "
                "include it as a component".format(
                    transform.model.plugin.name, mod_input
                )
            )
            continue

        # Skip tracers not in the control space
        if not hasattr(parameters, trcr) or not getattr(
            getattr(parameters, trcr, None), "iscontrol", False
        ):
            info(
                "{} is sensitive to {} as a {} "
                "but your control vector doesn't "
                "include it as a component".format(
                    transform.model.plugin.name, trcr, mod_input
                )
            )
            continue

        # Ignore initial conditions if not the first sub-period
        if mod_input == "inicond" and data.datei != min(di, df):
            info(
                "Initial conditions for {} are not considered "
                "for sub-period {} compared to general start date {}".format(
                    trcr, min(di, df), transform.datei
                )
            )
            continue

        # Check that the model provides information about sensitivities
        if "adj_out" not in datastore[tracer_id]:
            info(
                "Couldn't get any model sensitivity. "
                "Assuming zero sensitivity"
            )
            continue

        # Process other input types:
        # - re-project map sensitivities to control space
        # - sum date slices in the sensitivities to control space periods
        tracer = getattr(parameters, trcr)

        dslice = dateslice(tracer, di, df)

        # Loop over control space periods for temporal aggregation
        # Make vertical aggregation per temporal slice
        sensit = datastore[tracer_id]["adj_out"]
        sensit_data = sensit.data
        data_dates = sensit.time.data

        for idd, ds in enumerate(dslice):
            dd0 = tracer.dates[ds]
            dd1 = tracer.dates[min(ds + 1, tracer.ndates - 1)]

            # Either take the corresponding slice of time,
            # or take the exact date
            # if the control variable is on a time stamp
            if dd0 < dd1:
                mask = (data_dates >= np.datetime64(dd0)) & (
                    data_dates < np.datetime64(dd1)
                )
            else:
                mask = data_dates == np.datetime64(dd0)

            # For variables stored as a scaling factor,
            # scaling by the original value
            phys = np.ones((len(dslice)))
            if getattr(tracer, "type", "scalar") == "scalar":
                phys = tracer.read(
                    trcr,
                    tracer.dir,
                    tracer.file,
                    tracer.varname,
                    transform.model.input_dates[ddi][mask],
                    comp_type=mod_input,
                    tracer=tracer,
                    model=transform.model,
                    **kwargs
                ).data

            # Vertical aggregation
            vdata = np.sum(sensit_data[mask] * phys, axis=0)
            vaggreg = vmap2vaggreg(vdata[np.newaxis], tracer, tracer.domain)

            # 2d maps to control vector slices
            data.dx[tracer.xpointer: tracer.xpointer + tracer.dim][
                ds
                * tracer.hresoldim
                * tracer.vresoldim: (ds + 1)
                * tracer.hresoldim
                * tracer.vresoldim
            ] += map2scale(vaggreg, tracer, tracer.domain).flatten()

    return data

