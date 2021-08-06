from __future__ import division

import os

import numpy as np

from pycif.utils.check import info
from pycif.utils.geometry.dist_matrix import dist_matrix
from pycif.utils.netcdf import readnc
from pycif.utils.path import init_dir


def build_hcorrelations(
    zlat,
    zlon,
    lsm,
    sigma_land,
    sigma_sea,
    file_lsm=None,
    evalmin=0.5,
    dump=False,
    dir_dump="",
    projection="gps",
    **kwargs
):
    """Build horizontal correlation matrix based on distance between grid
    cells.
    For cells i and j, the corresponding correlation is:
    c(i,j) = exp(-dist(i, j) / sigma)
    sigma depends on the land-sea mask: land and sea cells are assumed
    un-correlated

    Args:
        zlat (np.array): 2D array of latitudes
        zlon (np.array): 2D array of longitudes
        file_lsm (str): path to NetCDF file with land-sea mask (grid must be
        consistent with LMDZ grid); the land-sea mask is assumed to be stored
        in the varible 'lsm'
        sigma_land (float): decay distance for correlation between land cells
        sigma_sea (float): idem for sea
        evalmin (float): flag out all eigenvalues below this value. Default
        is 0.5
        dump (bool): dumps computed correlations if True
        dir_dump (str): directory where correlation matrices are stored
        projection (str): the projection used for the longitudes and latitudes

    Return:
        tuple with:
            - square roots of eigenvalues
            - eigenvectors

    """

    # Define domain dimensions
    nlon, nlat = zlat.shape

    # Try reading existing file
    try:
        evalues, evectors = read_hcorr(
            nlon, nlat, sigma_sea, sigma_land, dir_dump
        )

    # Else build correlations
    except IOError:
        info("Computing hcorr")
        # No correlation between land and sea if lsm = True
        if lsm:
            landseamask = readnc(file_lsm, ["lsm"]).flatten()
            sigma = sigma_land * (landseamask[:, np.newaxis] >= 0.5) * (
                landseamask[np.newaxis, :] >= 0.5
            ) + sigma_sea * (landseamask[:, np.newaxis] < 0.5) * (
                landseamask[np.newaxis, :] < 0.5
            )  # Otherwise, isotropic correlation
        # Takes sigma_land
        else:
            sigma = sigma_land

        # Compute matrix of distance
        dx = dist_matrix(zlat, zlon, projection)

        # Compute the correlation matrix itself
        corr = np.exp(-dx / sigma)
        corr[sigma <= 0] = 0

        # Component analysis
        evalues, evectors = np.linalg.eigh(corr)

        # Re-ordering values
        # (not necessary in principle in recent numpy versions)
        index = np.argsort(evalues)[::-1]

        evalues = evalues[index]
        evectors = evectors[:, index]

        # Dumping to a txt file
        if dump:
            dump_hcorr(
                nlon, nlat, sigma_sea, sigma_land, evalues, evectors, dir_dump
            )

    except Exception as e:
        raise e

    # Truncating values < evalmin
    mask = evalues >= evalmin

    return evalues[mask] ** 0.5, evectors[:, mask]


def dump_hcorr(nlon, nlat, sigma_sea, sigma_land, evalues, evectors, dir_dump):
    """Dumps eigenvalues and vectors to a txt file. The default file format
    is:
    '{}/horcor{}x{}cs{}cl{}_lmdz5.txt'.format(
        dir_dump, nlon, nlat, sigma_sea, sigma_land)

    """

    ncell = evalues.size

    file_dump = "{}/horcor_{}x{}_cs{}_cl{}_lmdz.bin".format(
        dir_dump, nlon, nlat, sigma_sea, sigma_land
    )

    if os.path.isfile(file_dump):
        raise IOError(
            "Warning: {} already exists. "
            "I don't want to overwrite it".format(file_dump)
        )

    datasave = np.concatenate((evalues[np.newaxis, :], evectors), axis=0)

    # Creating path if does not exist
    if not os.path.isdir(os.path.dirname(file_dump)):
        init_dir(os.path.dirname(file_dump))

    # Saving data
    np.array(datasave).tofile(file_dump)


def read_hcorr(nlon, nlat, sigma_sea, sigma_land, dir_dump):
    """Reads horizontal correlations from existing text file

    Args:
        nlon, nlat (ints): dimensions of the domain
        sigma_land, sigma_sea (floats): horizontal correlation distances
        dir_dump (str): where the horizontal correlations have been stored

    """

    file_dump = "{}/horcor_{}x{}_cs{}_cl{}_lmdz.bin".format(
        dir_dump, nlon, nlat, sigma_sea, sigma_land
    )

    if not os.path.isfile(file_dump):
        raise IOError(
            "{} does not exist. "
            "Please compute correlations from scratch".format(file_dump)
        )

    data = np.fromfile(file_dump).reshape((-1, nlon * nlat))

    evalues = data[0]
    evectors = data[1:]

    evalues[evalues < 0] = 0.0

    return evalues, evectors
