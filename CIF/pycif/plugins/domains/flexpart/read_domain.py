import numpy as np
import pycif.utils.check as check


def read_grid(domain,
              **kwargs):
    """Reads a grid from an existing file

    Args:
        domain (Plugin): dictionary defining the domain. Should include
        ficgrid to be able to read the grid from a file

    Return:
        Grid dictionary with meshgrids for center lon/lat and corner lon/lat

    Notes:
        For now, this function is not used, so create_domain() will be called.


    """

    zlon = np.loadtxt(domain.ficlon)
    zlat = np.loadtxt(domain.ficlat)

    nlon = zlon.size
    nlat = zlat.size
    
    # Corner coordinates
    dlon = np.ptp(zlon) / (nlon - 1) / 2.
    zlonc = zlon - dlon
    zlonc = np.append(zlonc, zlonc[-1] + 2 * dlon)
    
    dlat = np.ptp(zlat) / (nlat - 1) / 2.
    zlatc = zlat - dlat
    zlatc = np.append(zlatc, zlatc[-1] + 2 * dlat)
    
    # Meshgrids
    zlon, zlat = np.meshgrid(zlon, zlat)
    zlonc, zlatc = np.meshgrid(zlonc, zlatc)
    
    # Saving information to domain attributes
    domain.nlon = nlon
    domain.nlat = nlat
    domain.zlon = zlon
    domain.zlat = zlat
    domain.zlonc = zlonc
    domain.zlatc = zlatc
