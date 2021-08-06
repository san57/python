import numpy as np
from builtins import map
from builtins import range


def read_grid(domain, **kwargs):
    """Reads a grid from an existing file

    Args:
        domain (Plugin): dictionary defining the domain. Should include
        filegrid to be able to read the grid from a file

    Return:
        Grid dictionary with meshgrids for center lon/lat and corner lon/lat

    Notes:
        In LMDZ the grid cell at -180 degW is repeated

    """

    with open(domain.filegrid, "r") as fgrid:
        # Loading geometry
        nlon = int(fgrid.readline()) + 1
        zlon = [float(fgrid.readline()) for _ in range(nlon - 1)]

        nlat = int(fgrid.readline())
        zlat = [float(fgrid.readline()) for _ in range(nlat)]

        zlon, zlat = np.meshgrid(zlon, zlat)

        # Loading dynamic and physics time steps wrt on-archive
        # One of them have to divide the other
        split = list(map(int, fgrid.readline().split()))

    # Corner coordinates
    zlonc = np.concatenate((zlon, zlon[:, np.newaxis, 1]), axis=1)
    zlonc = np.concatenate((zlonc, zlonc[-1, np.newaxis, :]), axis=0)
    zlonc -= 360.0 / (nlon - 1) / 2.0

    zlatc = zlat + 180.0 / (nlat - 1) / 2.0
    zlatc = np.concatenate((zlatc, -90.0 * np.ones((1, nlon - 1))), axis=0)
    zlatc[0, :] = 90
    zlatc = np.concatenate((zlatc, zlatc[:, -1, np.newaxis]), axis=1)

    # Saving information to domain attributes
    domain.zlon = zlon
    domain.zlat = zlat
    domain.zlonc = zlonc
    domain.zlatc = zlatc

    domain.dsplit = split[0]
    domain.psplit = split[1]
