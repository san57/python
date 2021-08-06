import numpy as np
import xarray as xr
from builtins import zip


def scale2map(x, tracer, dates, dom):
    """Translates a control vector slice to a map

    Args:
        x: the vector slice
        tracer: tracer Class with corresponding attributes
        dates: dates corresponding to the slice
        dom: the model domain

    Returns:
        xarray.Dataset with the map at native resolution
    """

    ndates = len(x)

    if tracer.hresol == "hpixels":
        xmap = x.reshape(-1, tracer.vresoldim, dom.nlat, dom.nlon)

    elif tracer.hresol == "bands":
        xmap = np.zeros((ndates, tracer.vresoldim, dom.nlat, dom.nlon))
        iband = 0

        for lat1, lat2 in zip(tracer.bands_lat[:-1], tracer.bands_lat[1:]):
            for lon1, lon2 in zip(tracer.bands_lon[:-1], tracer.bands_lon[1:]):
                reg = (
                    (dom.zlon >= lon1)
                    & (dom.zlon < lon2)
                    & (dom.zlat >= lat1)
                    & (dom.zlat < lat2)
                )
                xmap[..., reg] = x[..., iband, np.newaxis]
                iband += 1

    elif tracer.hresol == "ibands":
        xmap = np.zeros((ndates, tracer.vresoldim, dom.nlat, dom.nlon))
        iband = 0
        for i1, i2 in zip(tracer.bands_i[:-1], tracer.bands_i[1:]):
            for j1, j2 in zip(tracer.bands_j[:-1], tracer.bands_j[1:]):
                xmap[..., i1:i2, j1:j2] = x[..., iband, np.newaxis]
                iband += 1

    elif tracer.hresol == "regions":
        xmap = np.zeros((ndates, tracer.vresoldim, dom.nlat, dom.nlon))
        for regID, reg in enumerate(np.unique(tracer.regions)):
            xmap[..., tracer.regions == reg] = x[..., regID, np.newaxis]

    elif tracer.hresol == "global":
        xmap = x[..., np.newaxis] * np.ones((ndates, dom.nlat, dom.nlon))

    else:
        raise Exception("{} is not recognized".format(tracer.hresol))

    # Putting in xarray dataset
    xmap = xr.DataArray(
        xmap,
        coords={"time": dates, "lev": tracer.levels},
        dims=("time", "lev", "lat", "lon"),
    )

    return xmap


def map2scale(xmap, tracer, dom):
    """Translates a map to a control vector slice

    Args:
        xmap: the map
        tracer: tracer Class with corresponding attributes
        dom: the model domain

    Returns:
        the control vector vertical slice
    """

    # Checking that xmap has the correct dimension
    if not len(xmap.shape) == 4:
        raise Exception(
            "Warning! map2scale expects inputs data of dimension:"
            "(time, levels, lat, lon). Got only {} dimensions "
            "instead".format(len(xmap.shape))
        )

    ndates = xmap.shape[0]

    # Dealing different resolution types
    if tracer.hresol == "hpixels":
        x = xmap.reshape((ndates, -1))

    elif tracer.hresol == "bands":
        x = np.zeros((ndates, tracer.vresoldim, tracer.nbands))
        iband = 0
        for lat1, lat2 in zip(tracer.bands_lat[:-1], tracer.bands_lat[1:]):
            for lon1, lon2 in zip(tracer.bands_lon[:-1], tracer.bands_lon[1:]):
                reg = (
                    (dom.zlon >= lon1)
                    & (dom.zlon < lon2)
                    & (dom.zlat >= lat1)
                    & (dom.zlat < lat2)
                )
                x[..., iband] = np.sum(xmap[..., reg], axis=2)
                iband += 1

    elif tracer.hresol == "ibands":
        x = np.zeros((ndates, tracer.vresoldim, tracer.nbands))
        iband = 0
        for lat1, lat2 in zip(tracer.bands_i[:-1], tracer.bands_i[1:]):
            for lon1, lon2 in zip(tracer.bands_j[:-1], tracer.bands_j[1:]):
                reg = (
                    (dom.zlon >= lon1)
                    & (dom.zlon < lon2)
                    & (dom.zlat >= lat1)
                    & (dom.zlat < lat2)
                )
                x[..., iband] = np.sum(xmap[..., reg], axis=2)
                iband += 1

    elif tracer.hresol == 'regions':
        x = np.zeros((ndates, tracer.vresoldim, tracer.nregions))
        if getattr(tracer, "region_scale_area", False):
            for regID, reg in enumerate(np.unique(tracer.regions)):
                x[..., regID] = \
                    np.sum(xmap[..., tracer.regions == reg]
                           * tracer.domain.areas[..., tracer.regions == reg],
                           axis=2)
                x[..., regID] /= tracer.region_areas[regID]
        
        elif getattr(tracer, "region_max_val", False):
            for regID, reg in enumerate(np.unique(tracer.regions)):
                x[..., regID] = np.max(xmap[..., tracer.regions == reg])
        else:
            for regID, reg in enumerate(np.unique(tracer.regions)):
                x[..., regID] = np.sum(xmap[..., tracer.regions == reg], axis=2)

    elif tracer.hresol == "global":
        x = xmap.sum(axis=(2, 3))[..., np.newaxis]

    else:
        raise Exception("{} is not recognized".format(tracer.hresol))

    return x


def vmap2vaggreg(data, tracer, dom):
    """
    Aggregate full 3-dimensional field to vertical slices in the control vector

    :param data:
    :param tracer:
    :param dom:
    :return:
    """

    # Checking that xmap has the correct dimension
    if not len(data.shape) == 4:
        raise Exception(
            "Warning! vmap2vaggreg expects inputs data "
            "of dimension: (time, levels, lat, lon). "
            "Got only {} dimensions instead.".format(len(data.shape))
        )

    # Flat fields are directly returned
    if tracer.vresol == "column":
        return data.sum(axis=1)[:, np.newaxis, ...]

    elif tracer.vresol == "vpixels":
        return data

    elif tracer.vresol == "kbands":
        outshape = list(data.shape)
        outshape[1] = tracer.nbands
        outdata = np.zeros(outshape)

        iband = 0
        for k1, k2 in zip(tracer.kbands[:-1], tracer.kbands[1:]):
            outdata[:, iband, ...] = data[:, k1:k2, ...].sum(axis=1)
            iband += 1

        return outdata

    else:
        raise Exception("{} is not recognized".format(tracer.vresol))
