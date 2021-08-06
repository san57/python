import numpy as np
from netCDF4 import Dataset


def hresol2dim(tracer, dom, **kwargs):
    """Computes the horizontal size of a control vector from its resolution

    Args:
        tracer (Plugin): definition of the tracer, including the resolution and
                       additional information on the resolution
        domain (dict): the domain grid

    Returns
        int: the size of the control vector for this component

    """

    if tracer.hresol == "hpixels":
        return (
            dom.zlon.size
            if not getattr(tracer, "is_lbc", False)
            else dom.zlon_side.size
        )

    elif tracer.hresol == "bands":
        # Check that the lists are sorted
        if (
            not sorted(tracer.bands_lat) == tracer.bands_lat
            or not sorted(tracer.bands_lon) == tracer.bands_lon
        ):
            raise Exception(
                "Bands in the control vector are not sorted!\n"
                "Please check your configuration file."
            )

        # Compute the total number of dimensions in the grid
        if not hasattr(tracer, "nbands"):
            tracer.nbands = (len(tracer.bands_lat) - 1) * (
                len(tracer.bands_lon) - 1
            )
        return tracer.nbands

    elif tracer.hresol == "ibands":
        # Check that the lists are sorted
        if (
            not sorted(tracer.bands_i) == tracer.bands_i
            or not sorted(tracer.bands_j) == tracer.bands_j
        ):
            raise Exception(
                "Bands in the control vector are not sorted!\n"
                "Please check your configuration file."
            )

        # Compute the total number of dimensions in the grid
        if not hasattr(tracer, "nbands"):
            tracer.nbands = (len(tracer.bands_i) - 1) * (
                len(tracer.bands_j) - 1
            )
        return tracer.nbands

    elif tracer.hresol == "regions":
        if not hasattr(tracer, "regions"):
            with Dataset(tracer.fileregions, "r") as f:
                tracer.regions = f.variables["regions"][:]
            tracer.nregions = len(np.unique(tracer.regions))

            # Check that regions have the correct dimensions
            if tracer.regions.shape != (dom.nlat, dom.nlon):
                raise Exception(
                    "Regions were not correctly defined in {}".format(
                        tracer.fileregions
                    )
                )

        return tracer.nregions

    elif tracer.hresol == "global":
        return 1


def vresol2dim(tracer, dom, **kwargs):
    """Computes the horizontal size of a control vector from its resolution

    Args:
        tracer (Plugin): definition of the tracer, including the resolution and
                       additional information on the resolution
        domain (dict): the domain grid

    Returns
        int: the size of the control vector for this component

    """

    # Default vertical resolution is integrated columns
    tracer.vresol = getattr(tracer, "vresol", "column")
    tracer.nlev = getattr(tracer, "nlev", 1)

    # Loop over possible vertical resolution
    if tracer.vresol == "vpixels":
        tracer.levels = np.arange(tracer.nlev)
        return tracer.nlev

    elif tracer.vresol == "kbands":
        # Compute the total number of dimensions in the grid
        if not hasattr(tracer, "nbands"):
            tracer.nbands = len(tracer.kbands) - 1
        tracer.levels = np.array(tracer.kbands[:-1])
        return tracer.nbands

    elif tracer.vresol == "column":
        tracer.levels = np.array([0])
        return 1
