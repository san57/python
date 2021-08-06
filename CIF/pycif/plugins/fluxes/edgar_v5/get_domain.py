import datetime
import glob
import os

import xarray as xr

from pycif.utils.classes.setup import Setup


def get_domain(ref_dir, ref_file, input_dates, target_dir, tracer=None):

    # Looking for a reference file to read lon/lat in
    list_file = glob.glob("{}/*nc".format(ref_dir))

    domain_file = None
    # Either a file is specified in the Yaml
    if ref_file in list_file:
        domain_file = "{}/{}".format(ref_dir, ref_file)

    # Or loop over available file regarding file pattern
    else:
        for flx_file in list_file:
            try:
                date = datetime.datetime.strptime(
                    os.path.basename(flx_file), ref_file
                )
                domain_file = flx_file
                break
            except ValueError:
                continue

    if domain_file is None:
        raise Exception(
            "EDGARv5 domain could not be initialized as no file was found"
        )

    # Read lon/lat in
    nc = xr.open_dataset(domain_file, decode_times=False)

    lon = nc["lon"]
    lat = nc["lat"]

    lon_min = lon.min() - (lon[1] - lon[0]) / 2
    lon_max = lon.max() + (lon[-1] - lon[-2]) / 2
    lat_min = lat.min() - (lat[1] - lat[0]) / 2
    lat_max = lat.max() + (lat[-1] - lat[-2]) / 2

    nlon = lon.size
    nlat = lat.size

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
            }
        }
    )

    Setup.load_setup(setup, level=1)

    return setup.domain
