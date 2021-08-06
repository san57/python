import numpy as np
import xarray as xr

from .utils.dates import dateslice
from pycif.utils.dataarrays.reindex import reindex
from .utils.scalemaps import scale2map


def state2native(
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

    ddi = min(di, df)

    if data == {}:
        xout = transform.from_dict({})
        xout.datastore = {}
        xmod = {}

    else:
        xout = data
        xmod = data.datastore

    # Basic objects
    tracer_ids = mapper["outputs"]

    for trid in tracer_ids:
        tracer = tracer_ids[trid]

        comp = trid[0]
        trcr = trid[1]

        force_read = tracer.get("force_read", False)
        in_files = tracer["input_files_parameters"][ddi]
        in_dates = tracer["input_dates_parameters"][ddi]

        component = transform.orig_component_plg
        tracer = transform.orig_parameter_plg

        # Saving reference directories if specified
        xmod[trid] = {"tracer": tracer}
        xmod[trid]["fileorig"] = getattr(tracer, "file", None)
        xmod[trid]["dirorig"] = getattr(tracer, "dir", None)
        xmod[trid]["varname"] = getattr(tracer, "varname", None)

        # Skip parameters not in the control space
        if not getattr(tracer, "iscontrol", False):
            if force_read:
                inputs = tracer.read(
                    trcr,
                    "",
                    in_files,
                    tracer.varname,
                    in_dates,
                    comp_type=comp,
                    tracer=tracer,
                    model=transform.model,
                    **kwargs
                )
                xmod[trid]["spec"] = inputs

            xout.datastore = xmod
            continue

        # Otherwise, reformat from x

        # Translates control vector, and increments if tangent-linear
        variables = {"scale": data.x}
        if mode == "tl":
            variables["incr"] = data.dx

        # Deal with control vect dates and cut if controlvect period spans
        # outside the sub-simulation period
        dslice = dateslice(tracer, di, df)
        cdates = tracer.dates[dslice]
        if cdates[0] < ddi:
            cdates[0] = ddi

        # Translates only control variables corresponding to the
        # simulation period
        for x in variables:
            tmp = np.reshape(
                variables[x][tracer.xpointer: tracer.xpointer + tracer.dim],
                (tracer.ndates, tracer.nlev, -1),
            )[dslice]

            # Deals with different resolutions
            xmod[trid][x] = scale2map(tmp, tracer, cdates, tracer.domain)

        # Now deals with scalars and physical variables
        if getattr(tracer, "type", "scalar") == "scalar":
            # Read the tracer array and apply the present control vector
            # scaling factor
            inputs = tracer.read(
                trcr,
                tracer.dir,
                tracer.file,
                tracer.varname,
                transform.model.input_dates[ddi],
                comp_type=comp,
                tracer=tracer,
                model=transform.model,
                **kwargs
            )
            
            scale = reindex(
                xmod[trid]["scale"],
                levels={"time": inputs.time, "lev": inputs.lev},
            )
            xmod[trid]["spec"] = inputs * scale

            if mode == "tl":
                incr = reindex(
                    xmod[trid]["incr"],
                    levels={"time": inputs.time, "lev": inputs.lev},
                )
                xmod[trid]["incr"] = incr * inputs

        # Data already contains the correct info for physical control variables
        # WARNING: so far, assumes that the vertical resolution is already
        # correct
        elif getattr(tracer, "type", "scalar") == "physical":
            spec = reindex(
                xmod[trid]["scale"],
                levels={
                    "time": xr.DataArray(transform.model.input_dates[ddi]),
                    "lev": xmod[trid]["scale"].lev,
                },
            )
            xmod[trid]["spec"] = spec

            if mode == "tl":
                incr = reindex(
                    xmod[trid]["incr"],
                    levels={
                        "time": xr.DataArray(transform.model.input_dates[ddi]),
                        "lev": xmod[trid]["incr"].lev,
                    },
                )
                xmod[trid]["incr"] = incr

        # Removing the scaling factor as all information is stored in
        # 'spec' and 'incr' now
        xmod[trid].pop("scale")

        xout.datastore = xmod

    return xout
