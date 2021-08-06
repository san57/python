import numpy as np


def read_grid(domain, **kwargs):
    """Reads an horizontal and a vertical grid from existing files

    Args:
        domain (dictionary): dictionary defining the domain. Should include
        filegrid to be able to read the grid from a file

    Return:
        Grid dictionary with meshgrids for center lon/lat and corner lon/lat

    Notes:

    """

    # Loading geometry from the list of existing domains
    with open("{}/domainlist.nml".format(domain.repgrid), "r") as fgrid:
        ln = fgrid.readlines()
        for l in ln[1:]:
            s = l.split()
            if s[0] == domain.domid:
                nzo = int(s[1])
                nme = int(s[2])
                break

    # Saving information to domain attributes
    domain.nlat = nme
    domain.nlon = nzo

    # Reading lat/lon and latc/lonc
    file_hcoord = "{}/HCOORD/COORD_{}".format(domain.repgrid, domain.domid)
    data = np.genfromtxt(file_hcoord)
    lon = data[:, 0].reshape(nme, nzo)
    lat = data[:, 1].reshape(nme, nzo)

    file_hcoord = "{}/HCOORD/COORDcorner_{}".format(
        domain.repgrid, domain.domid
    )
    data = np.genfromtxt(file_hcoord)
    lonc = data[:, 0].reshape(nme + 1, nzo + 1)
    latc = data[:, 1].reshape(nme + 1, nzo + 1)

    # Putting the data into the domain
    domain.zlonc = lonc
    domain.zlatc = latc
    domain.zlon = lon
    domain.zlat = lat
    domain.corners = "{}/HCOORD/COORDcorner_{}".format(
        domain.repgrid, domain.domid
    )

    # Reading vertical coordinates
    file_vcoord = "{}/VCOORD/VCOORD_{}_{}_{}".format(
        domain.repgrid, domain.nlev, domain.p1, domain.pmax
    )
    data = np.genfromtxt(file_vcoord)
    domain.sigma_a = data[:, 0] * 1e5
    domain.sigma_b = data[:, 1] * 1e2

    domain.nlev = len(domain.sigma_a)
