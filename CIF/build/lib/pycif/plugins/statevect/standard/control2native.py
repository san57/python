from __future__ import absolute_import

import datetime

import numpy as np

from .utils.dates import dateslice
from .utils.reindex import reindex
from .utils.scalemaps import scale2map


def control2native(
    self, mod_input, di, df, mode, runsubdir, workdir, transforms, **kwargs
):
    """Translates information in the control vector to real size-data. This
    is in preparation to the generation of proper model inputs from the
    real-size data

    Args:
        self (Plugin): the control vect to use to generate inputs
        mod_input (str): type of inputs to prepare
        di (datetime.datetime): starting date of the simulation window
        df (datetime.datetime): ending date of the simulation window
        mode (str): running mode: one of 'fwd', 'adj' and 'tl'
        runsubdir (str): sub-directory for the current simulation
        workdir (str): the directory of the whole pyCIF simulation

    Returns:
        xarray Dataset with native resolution control variables

    """

    ddi = min(di, df)
    ddf = min(di, df)

    # Initializing the output with necessary information
    xout = self.from_dict({})
    xout.datastore = {}

    # If no mapper is available, just return an empty data set
    # Pass file information from the state vector to xout if available
    mapper = transforms.mapper[mod_input]

    if mapper == {} or mapper["input_parameter"] == []:
        if hasattr(self.components, mod_input):
            component = getattr(self.components, mod_input)
            if hasattr(component, "dir"):
                xout.dir = component.dir
            if hasattr(component, "file"):
                xout.file = component.file
        return xout

    # Loop over inputs to translate to the model resolution
    xmod = {}
    for trcr, comp_type, force_read, in_dates, in_files in zip(
        mapper["input_parameter"],
        mapper["input_component"],
        mapper["force_read"],
        mapper["input_dates_parameters"],
        mapper["input_files_parameters"],
    ):
        component = getattr(self.components, comp_type)
        tracer = getattr(component.parameters, trcr)
        trid = (trcr, comp_type)

        # Saving reference directories if specified
        xmod[trid] = {"tracer": tracer}
        xmod[trid]["fileorig"] = getattr(tracer, "file", None)
        xmod[trid]["dirorig"] = getattr(tracer, "dir", None)
        xmod[trid]["varname"] = getattr(tracer, "varname", None)

        # Saving component directories if available
        if hasattr(component, "dir"):
            xout.dir = component.dir
        if hasattr(component, "file"):
            xout.file = component.file

        # Skip parameters not in the control space
        if not tracer.iscontrol:
            if force_read:
                inputs = tracer.read(
                    trcr,
                    "",
                    in_files[ddi],
                    tracer.varname,
                    in_dates[ddi],
                    comp_type=mod_input,
                    tracer=tracer,
                    model=self.model,
                    **kwargs
                )
                xmod[trid]["spec"] = inputs

            continue

        # Translates control vector, and increments if tangent-linear
        variables = {"scale": self.x}
        if mode == "tl":
            variables["incr"] = self.dx

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
                self.model.input_dates[ddi],
                comp_type=mod_input,
                tracer=tracer,
                model=self.model,
                **kwargs
            )

            scale = reindex(
                xmod[trid],
                "scale",
                levels={"time": inputs.time, "lev": inputs.lev},
            )
            xmod[trid]["spec"] = inputs * scale

            if mode == "tl":
                incr = reindex(
                    xmod[trid],
                    "incr",
                    levels={"time": inputs.time, "lev": inputs.lev},
                )
                xmod[trid]["incr"] = incr * inputs

        # Data already contains the correct info for physical control variables
        # WARNING: so far, assumes that the vertical resolution is already
        # correct
        elif getattr(tracer, "type", "scalar") == "physical":
            spec = reindex(
                xmod[trid],
                "scale",
                levels={
                    "time": self.model.input_dates[ddi],
                    "lev": xmod[trid]["scale"].lev,
                },
            )
            xmod[trid]["spec"] = spec

            if mode == "tl":
                incr = reindex(
                    xmod[trid],
                    "incr",
                    levels={
                        "time": self.model.input_dates[ddi],
                        "lev": xmod[trid]["incr"].lev,
                    },
                )
                xmod[trid]["incr"] = incr

        # Removing the scaling factor as all information is stored in
        # 'spec' and 'incr' now
        xmod[trid].pop("scale")

    xout.datastore = xmod

    return xout
