import numpy as np

from pycif.utils.classes.setup import Setup
from .utils import grib_file_reader, find_valid_file


def get_domain(ref_dir, ref_file, input_dates, target_dir, tracer=None):

    date_ref = input_dates.values()[0][0]

    dir_dd = date_ref.strftime(ref_dir)
    files_3d, dates_3d = find_valid_file(dir_dd, tracer.file_3d, date_ref)
    lon, lat = grib_file_reader(files_3d[0], ["longitude", "latitude"])
    pv = grib_file_reader(files_3d[0], [], attribute="pv")

    lon_min = lon.min() - (lon[1] - lon[0]) / 2
    lon_max = lon.max() + (lon[-1] - lon[-2]) / 2
    lat_min = lat.min() - (lat[1] - lat[0]) / 2
    lat_max = lat.max() + (lat[-1] - lat[-2]) / 2

    nlon = lon.size
    nlat = lat.size

    # Reconstruct alpha and beta
    ecm_nlevs = len(pv) / 2 - 1
    sigma_a = np.empty(ecm_nlevs)
    sigma_b = np.empty(ecm_nlevs)
    for ii in range(ecm_nlevs - 1):
        sigma_a[ii] = (pv[ecm_nlevs - ii] + pv[ecm_nlevs - ii - 1]) / 2
        sigma_b[ii] = (
            100
            * (pv[1 + 2 * ecm_nlevs - ii] + pv[1 + 2 * ecm_nlevs - ii - 1])
            / 2
        )

    # Initializes domain
    setup = Setup.from_dict(
        {
            "domain": {
                "plugin": {
                    "name": "dummy",
                    "version": "std",
                    "type": "domain",
                },
                "xmin": lon_min,
                "xmax": lon_max,
                "ymin": lat_min,
                "ymax": lat_max,
                "nlon": nlon,
                "nlat": nlat,
                "nlev": ecm_nlevs,
                "sigma_a": sigma_a,
                "sigma_b": sigma_b,
            }
        }
    )

    Setup.load_setup(setup, level=1)

    return setup.domain
