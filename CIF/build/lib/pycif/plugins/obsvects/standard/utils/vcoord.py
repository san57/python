import numpy as np
from builtins import zip

from pycif.utils.check import info


def vcoord(obsvect, **kwargs):
    """Computes the vertical layer in which fall the observations
    """

    info("Finding model levels corresponding to observations")

    # Don't do anything if the datastore is empty
    if len(obsvect.datastore) == 0:
        return obsvect

    # If a file with fixed vertical coordinates is specified, use it
    if hasattr(obsvect, "file_statlev"):
        file_statlev = obsvect.file_statlev
        info(
            "Using pre-defined vertical coordinates for stations: {}".format(
                file_statlev
            )
        )
        obsvect.datastore = vcoordfromfile(
            obsvect.datastore, file_statlev, **kwargs
        )

    # Else compute vertical coordinates from meteo files
    # To be coded from pyCIF-CHIMERE
    # To be generalized by using meteo plugin
    else:
        obsvect.datastore = vcoordfrommeteo(
            obsvect.workdir, obsvect.datastore, **kwargs
        )

    return obsvect


def vcoordfromfile(datastore, file_lev, **kwargs):
    """Computes vertical levels from pre-defined file

    Args:
        datastore (dict): dictionary of pd.Dataframes containing measurement
                          informations
        file_lev (str): path to the file defining station levels

    """

    lev_stats = np.genfromtxt(file_lev, skip_header=1, usecols=0, dtype=str)
    lev_infos = np.genfromtxt(file_lev, skip_header=1)[:, 1:]

    datastore.loc[:, "level"] = np.nan

    for s, linfo in zip(lev_stats, lev_infos):
        mask = datastore["station"] == s.lower()
        datastore.loc[mask, "level"] = linfo[2]

    return datastore


def vcoordfrommeteo(workdir, datastore, **kwargs):
    meteodir = workdir + "/meteo/"

    datastore["level"] = np.nan

    # ds = xr.open_dataset(meteodir + 'fluxstoke.an2012.m01.nc')
    # alt = ds['phi'] / 9.81
    # alt = alt.mean(dim=['time_counter', 'lat', 'lon'])

    alt = np.array(
        [
            419.947858,
            491.729716,
            590.03352,
            731.746909,
            933.003129,
            1209.481853,
            1577.072424,
            2052.229597,
            2652.034343,
            3392.864731,
            4287.454615,
            5338.23227,
            6527.101582,
            7802.693814,
            9078.312502,
            10258.393071,
            11298.825452,
            12251.020931,
            13215.014871,
            14244.018004,
            15347.801437,
            16530.826455,
            17803.910146,
            19179.506405,
            20667.282816,
            22274.191707,
            24012.661236,
            25898.57798,
            27949.867586,
            30191.221244,
            32660.342496,
            35402.092974,
            38490.531506,
            41997.093232,
            46044.013939,
            50511.472241,
            55800.118528,
            61192.883408,
            72735.110678,
        ]
    )

    for index, row in datastore.iterrows():
        idx = (np.abs(alt - row["alt"])).argmin()
        datastore.at[index, "level"] = idx + 1

    return datastore
