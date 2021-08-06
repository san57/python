from __future__ import division

from osgeo import ogr, osr

import numpy as np
from builtins import range


def calc_areas(domain, **kwargs):
    """Compute grid cells surfaces

    Args:
        domain (Plugin): a domain dictionary with pre-loaded zlonc and zlatc
        **kwargs (dictionary): any extra arguments

    Returns:
        numpy.array with all areas in m2

    """

    # Reference GPS projection
    srsQuery = osr.SpatialReference()
    srsQuery.ImportFromEPSG(4326)

    # Grid corners
    zlonc = domain.zlonc
    zlatc = domain.zlatc
    nmerid = domain.nlat
    nzonal = domain.nlon

    # Loop over all cells
    areas = np.array(
        [
            [
                calc_cellarea(i, j, zlonc, zlatc, srsQuery)
                for j in range(nzonal)
            ]
            for i in range(nmerid)
        ]
    )

    domain.areas = areas


def calc_cellarea(i, j, zlonc, zlatc, srsQuery):
    # Create ring
    ring = ogr.Geometry(ogr.wkbLinearRing)
    _ = ring.AddPoint(zlonc[i, j], zlatc[i, j])
    _ = ring.AddPoint(zlonc[i + 1, j], zlatc[i + 1, j])
    _ = ring.AddPoint(zlonc[i + 1, j + 1], zlatc[i + 1, j + 1])
    _ = ring.AddPoint(zlonc[i, j + 1], zlatc[i, j + 1])
    _ = ring.AddPoint(zlonc[i, j], zlatc[i, j])

    # Define local projection
    srsArea = osr.SpatialReference()
    srsArea.ImportFromProj4(
        "+proj=laea +lat_0=80 +lon_0={} +x_0=0 +y_0=0 "
        "+ellps=WGS84 +units=m +no_defs ".format(zlatc[i, j], zlonc[i, j])
    )

    transf = osr.CoordinateTransformation(srsQuery, srsArea)

    ring.Transform(transf)

    # Returns area in m2
    return ring.GetArea()
