from __future__ import division
from __future__ import print_function

from osgeo import ogr

import numpy as np
from builtins import range
from builtins import zip
from future.moves.itertools import zip_longest
from pyproj import Proj, transform

from pycif.utils.check import info


def hcoord(obsvect, chunksize=1e3, **kwargs):
    """Finds out in which domain grid cells the observations are

    Args:
        chunksize (int): size of the datastore chunks used to compute
            coordinates; this allows the code to be faster and limiting memory
            issues when feeding the entire datastore to other methods; defaults
            is 1e3

    """

    info("Finding model grid cells corresponding to observations")

    # Don't do anything if the datastore is empty
    ds = obsvect.datastore
    if len(ds) == 0:
        return obsvect

    # Corner coordinates
    zlonc = obsvect.model.domain.zlonc
    zlatc = obsvect.model.domain.zlatc

    # Initialize i, j coordinates
    ds.loc[:, "i"] = np.nan
    ds.loc[:, "j"] = np.nan
    
    # Check domain regularity
    isregular = (
        np.sum(zlonc[0, np.newaxis] - zlonc) == 0
        and np.sum(zlatc[:, 0, np.newaxis] - zlatc) == 0
    )
    
    is_unstructured = getattr(obsvect.model.domain, "unstructured_domain",
                              False)
    
    # Makes simplified operations if regular
    if isregular or is_unstructured:
        lon = ds["lon"]
        lat = ds["lat"]

        listi = []
        listj = []

        # Loops over chunks to make the code faster
        for lon_chunk, lat_chunk in zip(
            zip_longest(*[iter(lon)] * int(chunksize), fillvalue=-999),
            zip_longest(*[iter(lat)] * int(chunksize), fillvalue=-999),
        ):
            lon_chunk = np.array(lon_chunk)[np.array(lon_chunk) != -999]
            lat_chunk = np.array(lat_chunk)[np.array(lat_chunk) != -999]

            discont = (
                180
                if getattr(obsvect.model.domain, "projection", "gps") == "gps"
                else np.ptp(zlonc)
            )
            
            if isregular:
                i, j = find_gridcell(
                    lon_chunk,
                    lat_chunk,
                    zlonc,
                    zlatc,
                    isregular=isregular,
                    discont=discont,
                )
            
            else:
                # For unstructured domains, take closest point
                # TODO: generalize
                dist = (obsvect.model.domain.zlon
                        - lon_chunk) ** 2 \
                    + (obsvect.model.domain.zlat
                     - lat_chunk) ** 2
                
                i = np.argmin(dist, axis=0)
                j = np.zeros(len(lon_chunk))
                
            listi.extend(list(i))
            listj.extend(list(j))
        
        ds.loc[:, "i"] = listi
        ds.loc[:, "j"] = listj
    
    else:
        # Remove duplicates to reduce the number of stations to locate
        locations = ds[["lon", "lat"]].drop_duplicates(subset=["lon", "lat"])

        # Loop over (lon, lat) tuples
        k = 0
        nlocs = np.floor(locations.size / 10) + 1
        for lon, lat in zip(locations["lon"], locations["lat"]):
            i, j = find_gridcell(lon, lat, zlonc, zlatc,
                                 isregular=isregular,
                                 is_unstructured=is_unstructured)

            ds.loc[(ds["lon"] == lon) * (ds["lat"] == lat),
                   ["i", "j"]] = [i, j]

            k += 1
            if k % nlocs == 0:
                info("{}%".format(k * 10.0 / nlocs))

    # Cropping observations outside the model domain
    mask = ~np.isnan(ds["i"]) & ~np.isnan(ds["j"])
    obsvect.datastore = ds.loc[mask]

    return obsvect


def find_gridcell(
    lon,
    lat,
    zlonc,
    zlatc,
    orig_proj=Proj(init="epsg:4326"),
    isregular=False,
    is_unstructured=False,
    discont=180,
):
    """Finds the grid cell corresponding to a coordinate

    Args:
        lon (np.array): longitude of the point to find
        lat (np.array): latitude of the point to find
        zlonc (np.array): longitudes of the corners of the grid
        zlatc (np.array): latitudes of the corners of the grid
        orig_proj (projection): the projection for reporting coordinates.
                                Default is WSG84
        isregular (Boolean): if True, simplified operations are computed.
                             Default is False

    Returns:
        i, j the grid cell ID

    Notes: For very regular grids, this script could be made more effileient
    and shorter, but the objective is to be able to deal with any domain

    """

    # If regular domain make simplified operations
    if isregular and not is_unstructured:
        try:
            xlon = 1.0 / np.unwrap(
                zlonc[0, np.newaxis] - lon[:, np.newaxis], discont=discont
            )
            xlat = 1.0 / (lat[:, np.newaxis] - zlatc[:, 0])
            return xlat.argmin(axis=1), xlon.argmin(axis=1)

        # If only one (lon, lat) tuple was provided
        except TypeError as e:
            info(e)
            xlon = 1.0 / np.unwrap(zlonc[0] - lon, discont=discont)
            xlat = 1.0 / (lat - zlatc[:, 0])
            return xlat.argmin(), xlon.argmin()

    # Grid shape
    nlat, nlon = zlonc.shape

    # Define measurement point geometry in local coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(0, 0)

    # Define local projection
    target_proj = Proj(
        "+proj=laea +lat_0={} +lon_0={} +x_0=0 +y_0=0 "
        "+ellps=WGS84 +units=m +no_defs ".format(lat, lon)
    )

    zxc, zyc = transform(
        orig_proj, target_proj, (zlonc + 180) % 360 - 180, zlatc
    )
    dist = zxc ** 2 + zyc ** 2

    imin, jmin = np.unravel_index(dist.argmin(), dist.shape)

    imin = imin % nlat
    jmin = jmin % nlon

    for i in range(imin - 3, imin + 2):
        for j in range(jmin - 3, jmin + 2):
            # Looping if the minimum is close to the domain side
            i = max(min(i, nlat - 2), 0)
            j = max(min(j, nlon - 2), 0)

            # Create ring
            ring = ogr.Geometry(ogr.wkbLinearRing)
            _ = ring.AddPoint(zxc[i, j], zyc[i, j])
            _ = ring.AddPoint(zxc[i + 1, j], zyc[i + 1, j])
            _ = ring.AddPoint(zxc[i + 1, j + 1], zyc[i + 1, j + 1])
            _ = ring.AddPoint(zxc[i, j + 1], zyc[i, j + 1])
            _ = ring.AddPoint(zxc[i, j], zyc[i, j])

            # Create polygon
            poly = ogr.Geometry(ogr.wkbPolygon)
            _ = poly.AddGeometry(ring)

            # Buffering the polygon to include points on the edge
            buffer = 1.0
            poly = poly.Buffer(buffer)

            if point.Within(poly):
                return i, j

    # If no grid cell was fund, raise exception
    #
    # raise IndexError(
    #     "No index was found for measurement at {}, {}".format(lon, lat))
    # TODO: Decide whether we raise exception
    #   when the observation is outside the domain
    return np.nan, np.nan
