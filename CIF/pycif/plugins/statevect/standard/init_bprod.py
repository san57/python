import numpy as np

from .build_hcorr import build_hcorrelations
from .build_tcorr import build_tcorrelations


def init_bprod(cntrlv, options={}, **kwargs):
    """Initilializes the product of chi by sqrt-B. It allows translating
    information from the minimization space to the control space. This first
    needs to initialize correlation matrices

    Args:
        cntrlv (dict): definition of the control vector

    Returns:
        updated control vector:

    """

    # Get the grid for horizontal correlations
    grid = cntrlv.domain

    # Initializing sqrt-B
    cntrlv.hcorrelations = {}
    cntrlv.tcorrelations = {}
    cntrlv.chi_dim = 0
    components = cntrlv.components
    for comp in components.attributes:
        component = getattr(components, comp)

        # Skip if component does not have parameters
        if not hasattr(component, "parameters"):
            continue

        for trcr in component.parameters.attributes:
            tracer = getattr(component.parameters, trcr)

            # Skip tracers that are not control variables
            if not hasattr(tracer, "hresol"):
                continue

            # Deals with horizontal correlations if any
            if hasattr(tracer, "hcorrelations") and tracer.hresol == "hpixels":
                corr = tracer.hcorrelations
                dump_hcorr = getattr(corr, "dump_hcorr", False)
                dircorrel = getattr(corr, "dircorrel", "")
                evalmin = getattr(corr, "evalmin", 0.5)
                projection = getattr(grid, "projection", "gps")

                # Two possible options: - uniform correlation length,
                # or separated land and sea, along a land-sea mask
                # Default is no separation
                lsm = getattr(corr, "landsea", False)
                if lsm:
                    sigma_land = getattr(corr, "sigma_land", -1)
                    sigma_sea = getattr(corr, "sigma_sea", -1)
                    file_lsm = getattr(corr, "filelsm")

                else:
                    sigma = getattr(corr, "sigma", -1)
                    sigma_land = sigma
                    sigma_sea = -999
                    file_lsm = None

                # Load or compute the horizontal correlations
                # Checks before whether they were already loaded
                if (sigma_land, sigma_sea) in cntrlv.hcorrelations:
                    sqrt_evalues = cntrlv.hcorrelations[
                        (sigma_land, sigma_sea)
                    ]["sqrt_evalues"]
                    evectors = cntrlv.hcorrelations[(sigma_land, sigma_sea)][
                        "evectors"
                    ]

                else:
                    sqrt_evalues, evectors = build_hcorrelations(
                        grid.zlat,
                        grid.zlon,
                        lsm,
                        sigma_land,
                        sigma_sea,
                        file_lsm=file_lsm,
                        evalmin=evalmin,
                        dump=dump_hcorr,
                        dir_dump=dircorrel,
                        projection=projection,
                        **kwargs
                    )

                    # Storing computed correlations for use by other components
                    cntrlv.hcorrelations[(sigma_land, sigma_sea)] = {
                        "evectors": evectors,
                        "sqrt_evalues": sqrt_evalues,
                    }

                corr.sqrt_evalues = sqrt_evalues
                corr.evectors = evectors

                tracer.chi_hresoldim = sqrt_evalues.size

            else:
                tracer.chi_hresoldim = tracer.hresoldim

            # Initializes temporal correlations if any
            if hasattr(tracer, "tcorrelations"):
                corr = tracer.tcorrelations
                sigma_t = getattr(corr, "sigma_t", -1)
                dump_tcorr = getattr(corr, "dump_tcorr", False)
                dircorrel = getattr(corr, "dircorrel", "")
                period = getattr(tracer, "period")
                dates = getattr(tracer, "dates")
                evalmin = getattr(corr, "evalmin", 1e-5)

                # Load or compute the temporal correlations
                # Checks before whether they were already loaded
                if (sigma_t, period) in cntrlv.tcorrelations:
                    sqrt_evalues = cntrlv.tcorrelations[(sigma_t, period)][
                        "sqrt_evalues"
                    ]
                    evectors = cntrlv.tcorrelations[(sigma_t, period)][
                        "evectors"
                    ]

                else:
                    # Computes correlations for dates[:-1]
                    # as the last dates indicates the end of the period
                    sqrt_evalues, evectors = build_tcorrelations(
                        period,
                        dates,
                        sigma_t,
                        evalmin=evalmin,
                        dump=dump_tcorr,
                        dir_dump=dircorrel,
                        **kwargs
                    )

                    # Storing computed correlations for use by other components
                    cntrlv.tcorrelations[(sigma_t, period)] = {
                        "evectors": evectors,
                        "sqrt_evalues": sqrt_evalues,
                    }

                corr.sqrt_evalues = sqrt_evalues
                corr.evectors = evectors

                tracer.chi_tresoldim = sqrt_evalues.size

            else:
                tracer.chi_tresoldim = tracer.ndates
            
            # TODO: dealing with vertical resolution
            tracer.chi_vresoldim = tracer.vresoldim
            
            # Update the dimension of the full control vector
            tracer.chi_pointer = cntrlv.chi_dim
            tracer.chi_dim = \
                tracer.chi_hresoldim \
                * tracer.chi_tresoldim \
                * tracer.chi_vresoldim
            cntrlv.chi_dim += tracer.chi_dim

    # Defining chi from the total dimension
    cntrlv.chi = np.zeros((cntrlv.chi_dim,))

    return cntrlv
