import numpy as np
import pandas as pd


def make(meteo, meteo_file, **kwargs):
    """Creates random meteo time series for toy Gaussian model

    Args:
        meteo (pycif.utils.classes.Meteo): dictionary defining the domain.

    Returns:
         initialized meteo

    """

    # Default temporal resolution of 1H
    if not hasattr(meteo, "resolution"):
        meteo.resolution = "1H"

    dates = pd.date_range(meteo.datei, meteo.datef, freq=meteo.resolution)
    ndates = len(dates)

    windspeed = np.random.uniform(0, 5, size=ndates)
    winddir = np.random.uniform(0, 360, size=ndates)
    stabclass = np.random.choice(["A", "B", "C", "D", "E", "F"], ndates)

    df = pd.DataFrame(
        data={
            "windspeed": windspeed,
            "winddir": winddir,
            "stabclass": stabclass,
        },
        index=dates,
    )
    df.index.rename("date", inplace=True)

    # Save to file
    df.to_csv(meteo_file)
