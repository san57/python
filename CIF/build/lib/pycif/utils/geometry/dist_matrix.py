from __future__ import division

import numpy as np


def dist_matrix(zlat, zlon, projection="gps"):
    """Computes the distance matrix for two arrays of longitudes and
    latitudes

    Args:
        zlat (np.array): numpy array of latitudes
        zlon (np.array): numpy array of longitudes
        projection (str): the projection used for the longitudes and latitudes

    Returns:
        np.array((zlat.size, zlat.size)): matrix of distances

    """

    if zlat.size != zlon.size:
        raise ValueError(
            "Warning: longitudes and latitudes do not have the "
            "same dimension. Cannot compute the distance matrix"
        )

    # Dealing differently with gps and xy projections
    if projection == "gps":
        # Earth radius
        rearth = 6371.03

        # Flatten lon/lat
        radlat = np.radians(zlat).flatten()
        radlon = np.radians(zlon).flatten()

        # Compute the distance
        val = 0.5 * (
            1
            - np.sin(radlat[:, np.newaxis]) * np.sin(radlat[np.newaxis, :])
            - np.cos(radlat[:, np.newaxis])
            * np.cos(radlat[np.newaxis, :])
            * np.cos(radlon[:, np.newaxis] - radlon[np.newaxis, :])
        )
        val[val < 0] = 0

        dx = rearth * 2 * np.arcsin(val ** 0.5)

    elif projection == "xy":
        lat = zlat.flatten()
        lon = zlon.flatten()

        dlat = (lat[:, np.newaxis] - lat[np.newaxis, :]) ** 2
        dlon = (lon[:, np.newaxis] - lon[np.newaxis, :]) ** 2
        dx = dlat + dlon
        dx = np.sqrt(dx)

    else:
        raise ValueError("Projection {} is not recognized".format(projection))

    return dx
