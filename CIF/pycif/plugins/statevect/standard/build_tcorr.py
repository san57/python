import os

import numpy as np
import pandas as pd

from pycif.utils.check import info


def build_tcorrelations(
    period, dates, sigma_t, evalmin=0.5, dump=False, dir_dump="", **kwargs
):
    """Build temporal correlation matrix based on timedelta between periods.
    For period i and j, the corresponding correlation is:
    c(i,j) = exp(-timedelta(i, j) / sigma)

    Args:
        period (int): period duration
        dates (np.array): dates sub-dividing the control vector periods
        sigma_t (float): decay distance for correlation between periods
                         (in days)
        evalmin (float): flag out all eigenvalues below this value.
                         Default is 0.5
        dump (bool): dumps computed correlations if True
        dir_dump (str): directory where correlation matrices are stored

    Return:
        tuple with:
            - square roots of eigenvalues
            - eigenvectors

    """

    # Try reading existing file
    try:
        evalues, evectors = read_tcorr(period, dates, sigma_t, dir_dump)

    # Else build correlations
    except IOError:
        info("Computing temporal correlations")

        # Compute matrix of distance
        dt = (
            pd.DatetimeIndex(dates).values[:, np.newaxis]
            - pd.DatetimeIndex(dates).values[np.newaxis, :]
        ) / np.timedelta64(sigma_t, "h")

        # Compute the correlation matrix itself
        corr = np.exp(-(dt ** 2))

        # Component analysis
        evalues, evectors = np.linalg.eigh(corr)

        # Re-ordering values
        # (not necessary in principle in recent numpy versions)
        index = np.argsort(evalues)[::-1]

        evalues = evalues[index]
        evectors = evectors[:, index]

        # Dumping to a txt file
        if dump:
            dump_tcorr(period, dates, sigma_t, evalues, evectors, dir_dump)

    except Exception as e:
        raise e

    # Truncating values < evalmin
    mask = evalues >= evalmin

    return evalues[mask] ** 0.5, evectors[:, mask]


def dump_tcorr(period, dates, sigma_t, evalues, evectors, dir_dump):
    """Dumps eigenvalues and vectors to a bin file. The default file format
    is:
    '{}/tempcor_{}_{}_per{}_ct{}_lmdz.bin'.format(
        dir_dump, datei, datef, period, sigma_t)

    Args:
        period (int): period duration
        dates (np.array): dates sub-dividing the control vector periods
        sigma_t (float): decay distance for correlation between periods
                         (in days)

    """

    datei = dates[0]
    datef = dates[-1]

    file_dump = "{}/tempcor_{}_{}_per{}_ct{}_lmdz.bin".format(
        dir_dump,
        datei.strftime("%Y%m%d%H%M"),
        datef.strftime("%Y%m%d%H%M"),
        period,
        sigma_t,
    )

    if os.path.isfile(file_dump):
        raise IOError(
            "Warning: {} already exists. "
            "I don't want to overwrite it".format(file_dump)
        )

    datasave = np.concatenate((evalues[np.newaxis, :], evectors), axis=0)
    datasave.tofile(file_dump)


def read_tcorr(period, dates, sigma_t, dir_dump):
    """Reads temporal correlations from existing bin file

    Args:
        period (int): period duration
        dates (np.array): dates sub-dividing the control vector periods
        sigma_t (float): decay distance for correlation between periods
                         (in days)
        dir_dump (str): where the horizontal correlations have been stored

    """

    datei = dates[0]
    datef = dates[-1]

    file_dump = "{}/tempcor_{}_{}_per{}_ct{}_lmdz.bin".format(
        dir_dump,
        datei.strftime("%Y%m%d%H%M"),
        datef.strftime("%Y%m%d%H%M"),
        period,
        sigma_t,
    )

    if not os.path.isfile(file_dump):
        raise IOError(
            "{} does not exist. "
            "Please compute correlations from scratch".format(file_dump)
        )

    data = np.fromfile(file_dump).reshape((-1, len(dates)))

    evalues = data[0]
    evectors = data[1:]

    evalues[evalues < 0] = 0.0

    return evalues, evectors
