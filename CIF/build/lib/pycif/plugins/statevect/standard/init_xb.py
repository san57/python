import numpy as np
from scipy import ndimage

from pycif.utils import dates
from .utils.dimensions import hresol2dim, vresol2dim
from .utils.scalemaps import map2scale, vmap2vaggreg


def init_xb(cntrlv, **kwargs):
    """Initializes the prior control vector. Loops over all components and
    tracers and process temporal and horizontal resolution.

    Args:
        cntrlv (Plugin): definition of the control vector.
        datei (datetime): initial date of the inversion window
        datei (datetime): end date of the inversion window

    """

    datei = getattr(cntrlv, "datei")
    datef = getattr(cntrlv, "datef")

    cntrlv.dim = 0
    cntrlv.xb = np.ones(0)
    cntrlv.std = np.ones(0)

    # If no definition is specified for the control vector in the Yaml,
    # return empty control vector
    if not hasattr(cntrlv, "components"):
        return cntrlv

    # Else, carry on initializing
    components = cntrlv.components
    for comp in components.attributes:
        component = getattr(components, comp)

        # Skip if component does not have parameters
        if not hasattr(component, "parameters"):
            continue

        for trcr in component.parameters.attributes:
            tracer = getattr(component.parameters, trcr)

            # Skip tracers that are not control variables
            tracer.iscontrol = hasattr(tracer, "hresol")
            if not tracer.iscontrol:
                continue

            # Deals with temporal resolution of the control component
            # A negative period is equivalent to no period
            tracer.dates = dates.date_range(
                datei,
                datef,
                getattr(tracer, "period", ""),
                subperiod=getattr(tracer, "subperiod", ""),
            )
            tracer.ndates = len(tracer.dates)

            # Keeping a pointer to the correct location in the whole control
            tracer.xpointer = cntrlv.dim

            # Updating dimension
            tracer.hresoldim = hresol2dim(tracer, tracer.domain, **kwargs)
            tracer.vresoldim = vresol2dim(tracer, tracer.domain, **kwargs)
            tracer.dim = tracer.ndates * tracer.hresoldim * tracer.vresoldim
            cntrlv.dim += tracer.dim

            # Scale xb if prescribed in Yaml
            xb_scale = getattr(tracer, "xb_scale", 1.0)

            # Filling with defined value if prescribed
            xb_value = getattr(tracer, "xb_value", 0.0)

            # Filling Xb and uncertainties
            # Most types are scalar factors
            if (
                getattr(tracer, "type", "scalar") == "scalar"
                or getattr(tracer, "hresol", "") != "hpixels"
                or getattr(tracer, "vresol", "") != "vpixels"
            ):
                tracer.type = "scalar"
                xb = np.ones(tracer.dim) * xb_scale + xb_value
                std = tracer.err * xb

            # Physical variables stores uncertainties in the physical space
            # Valid only for control variables at the pixel resolution
            else:
                tracer.type = "physical"
                flxall = (
                    tracer.read(
                        trcr,
                        tracer.dir,
                        tracer.file,
                        tracer.varname,
                        tracer.dates,
                        comp_type=comp,
                        tracer=tracer,
                        **kwargs
                    ).data
                    * xb_scale
                    + xb_value
                )

                # Vertical aggregation depending on vertical stacks
                vaggreg = vmap2vaggreg(flxall, tracer, tracer.domain)
                xstack = map2scale(vaggreg, tracer, cntrlv.domain)

                # Putting flatten values into xb
                xb = xstack.flatten()

                # Filling uncertainties depending on uncertainty type
                # If 'max', takes the maximum of neighbouring cells as
                # reference flux (spatially and temporally)
                if getattr(tracer, "errtype", "") == "max":
                    xstack = map2scale(
                        ndimage.maximum_filter(np.abs(vaggreg), size=3),
                        tracer,
                        cntrlv.domain,
                    )

                # If errtype is 'avg', prescribes uniform uncertainties
                elif getattr(tracer, "errtype", "") == "avg":
                    xstack = map2scale(
                        0.0 * np.abs(vaggreg) + np.mean(np.abs(vaggreg)),
                        tracer,
                        cntrlv.domain,
                    )

                std = tracer.err * np.abs(xstack).flatten()

            # Appending local xb to the control vector
            cntrlv.xb = np.append(cntrlv.xb, xb)
            cntrlv.std = np.append(cntrlv.std, std)

    return cntrlv
