import numpy as np


def create_domain(domain, **kwargs):
    """Creates a grid if needed

    Args:
        domain (dictionary): dictionary defining the domain.

    Returns:
         initialized domain

    """

    nlon = domain.nlon
    nlat = domain.nlat

    lonc = np.linspace(domain.xmin, domain.xmax, nlon + 1)
    latc = np.linspace(domain.ymin, domain.ymax, nlat + 1)

    lon = 0.5 * (lonc[1:] + lonc[:-1])
    lat = 0.5 * (latc[1:] + latc[:-1])

    # Meshgrids
    zlon, zlat = np.meshgrid(lon, lat)
    zlonc, zlatc = np.meshgrid(lonc, latc)

    domain.zlon = zlon
    domain.zlat = zlat
    domain.zlonc = zlonc
    domain.zlatc = zlatc
