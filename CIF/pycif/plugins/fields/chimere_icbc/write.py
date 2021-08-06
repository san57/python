import os

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from pycif.utils.netcdf import save_nc


def write(self, name, lbc_file, data, mode="a", comp_type=None):
    """Write flux to INICOND or BOUN_CONC CHIMERE compatible files.
    """

    # If mode is 'a' but file does not exit, switch to mode 'w'
    if mode == "a" and not os.path.isfile(lbc_file):
        mode = "w"

    # Loading longitudes and latitudes
    lon = self.domain.zlon
    lat = self.domain.zlat
    lon_side = self.domain.zlon_side
    lat_side = self.domain.zlat_side
    nlev = self.domain.nlev

    # Write INI_CONCS
    if comp_type == "inicond":
        write_iniconcs(name, lbc_file, data, lon, lat, mode)

    else:
        write_bounconcs(
            name,
            lbc_file,
            data,
            lon,
            lat,
            lon_side,
            lat_side,
            nlev,
            mode,
            comp_type,
        )


def write_iniconcs(name, lbc_file, data, lon, lat, mode="a"):
    nhours, nlev, nmerid, nzonal = np.shape(data)

    # Dimensions
    spstrlen = 23
    dimnames = [
        "Time",
        "south_north",
        "west_east",
        "bottom_top",
        "SpStrLen",
        "DateStrLen",
        "Species",
    ]
    dimlens = [None, nmerid, nzonal, nlev, spstrlen, 19, 1]

    # Variables names, dimension and attributes
    varnames = ["Times", "species", "lon", "lat", name]
    vardims = [
        ("Time", "DateStrLen"),
        ("Species", "SpStrLen"),
        ("south_north", "west_east"),
        ("south_north", "west_east"),
        ("Time", "bottom_top", "south_north", "west_east"),
    ]
    dtypes = ["c", "c", "f", "f", "d"]
    units = ["", "", "degrees_east", "degrees_north", "ppb"]
    attributes = [
        {},
        {},
        {"long_name": "Longitude"},
        {"long_name": "Latitude"},
        {"long_name": "{}".format(name)},
    ]

    # Variables to save
    times = [
        list(pd.to_datetime(d).strftime("%Y-%m-%d_%H:00:00"))
        for d in data.time.values
    ]
    specs = [list(name.ljust(spstrlen))]

    variables = [times, specs, lon, lat, data]

    save_nc(
        lbc_file,
        variables,
        varnames,
        vardims,
        dimnames,
        dimlens,
        units,
        dtypes,
        mode=mode,
        attributes=attributes,
        format="NETCDF3_CLASSIC",
    )


def write_bounconcs(
    name,
    lbc_file,
    data,
    lon,
    lat,
    lon_side,
    lat_side,
    nlev,
    mode="a",
    comp_type=None,
):

    nhours = len(data)
    nsides = len(lon_side)
    nmerid, nzonal = np.shape(lon)

    # Dimensions
    spstrlen = 23
    datestrlen = 19
    dimnames = [
        "Time",
        "south_north",
        "west_east",
        "bottom_top",
        "h_boundary",
        "SpStrLen",
        "DateStrLen",
        "Species",
    ]
    dimlens = [None, nmerid, nzonal, nlev, nsides, spstrlen, datestrlen, 1]

    # Variables names, dimension and attributes
    varnames = ["Times", "species", "lon", "lat", "top_conc", "lat_conc"]
    vardims = [
        ("Time", "DateStrLen"),
        ("Species", "SpStrLen"),
        ("south_north", "west_east"),
        ("south_north", "west_east"),
        ("Time", "south_north", "west_east", "Species"),
        ("Time", "bottom_top", "h_boundary", "Species"),
    ]
    dtypes = ["c", "c", "f", "f", "d", "d"]
    units = ["", "", "degrees_east", "degrees_north", "ppb", "ppb"]
    attributes = [
        {},
        {},
        {"long_name": "Longitude"},
        {"long_name": "Latitude"},
        {"long_name": "{}".format(varnames[-2])},
        {"long_name": "{}".format(varnames[-1])},
    ]

    # Variables to save
    times = [
        list(pd.to_datetime(d).strftime("%Y-%m-%d_%H:00:00"))
        for d in data.time.values
    ]

    # Output data
    top_in = (
        data[:, 0]
        if comp_type == "topcond"
        else np.zeros((nhours, nmerid, nzonal, 1))
    )
    lat_in = (
        data.data[..., 0, :, np.newaxis]
        if comp_type == "latcond"
        else np.zeros((nhours, nlev, nsides, 1))
    )
    ljust_specs_in = [list(name.ljust(spstrlen))]
    specs_in = [name]

    # Append species to existing file
    if os.path.isfile(lbc_file) and mode == "a":
        with Dataset(lbc_file, "a") as f:
            ljust_specs_in = f.variables["species"][:]
            specs_in = ["".join(p).strip() for p in ljust_specs_in]
            top_in = f.variables["top_conc"][:]
            lat_in = f.variables["lat_conc"][:]

        if name in specs_in:
            ispec = specs_in.index(name)
            if comp_type == "topcond":
                top_in[..., ispec] = data[:, 0]
            else:
                lat_in[..., ispec] = data[..., 0, :]

        else:
            specs_in = specs_in + [name]
            ljust_specs_in = [list(s.ljust(spstrlen)) for s in specs_in]
            if comp_type == "topcond":
                top_in = np.concatenate(top_in, data, axis=3)
                lat_in = np.concatenate(
                    lat_in, np.zeros((nhours, nlev, nsides, 1)), axis=3
                )

            else:
                top_in = np.concatenate(
                    top_in, np.zeros((nhours, nmerid, nzonal, 1)), axis=3
                )
                lat_in = np.concatenate(lat_in, data, axis=3)

    variables = [times, ljust_specs_in, lon, lat, top_in, lat_in]

    save_nc(
        lbc_file,
        variables,
        varnames,
        vardims,
        dimnames,
        dimlens,
        units,
        dtypes,
        mode=mode,
        attributes=attributes,
        format="NETCDF3_CLASSIC",
    )
