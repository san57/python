import copy
import os
import shutil

from netCDF4 import Dataset


def readnc(nc_file, variables):
    """Extract variables from a NetCDF file

    Args:
        file (str): path to the NetCDF file
        variables (list[str]): list of variables to extract

    Returns:
        list([np.array]): list of variable values

    """

    outvars = []

    with Dataset(nc_file, "r") as f:
        for var in variables:
            if var not in f.variables:
                raise ValueError(
                    "{} is not present in {}. Please check your "
                    "NetCDF file and variable names".format(var, nc_file)
                )
            outvars.append(f.variables[var][:])

    # If only one element, return it instead of a list-like output
    if len(outvars) == 1:
        outvars = outvars[0]

    return outvars


def save_nc(
    fout,
    variables,
    varnames,
    vardims,
    dimnames,
    dimlens,
    units,
    dtypes,
    format="NETCDF4",
    mode="w",
    attributes=None,
):
    """Saves variables to a NetCDF file

    Args:
        fout (str): file to create or to append
        variables (list[np.array]): list of variables to append
        varnames (list[str]): list of variable names
        vardims (list[tuple[str]]): list of tuples with
            the variable dimension names
        dimnames (list[str]): dimensions names
        units (list[str]): variables units
        dtypes (list[str]): variables types
        format: NetCDF format
        mode: IO mode; one of 'r', 'w' and 'a'; default is 'w' (write)

    """

    if os.path.isfile(fout) and mode == "w":
        shutil.rmtree(fout, ignore_errors=True)

    with Dataset(fout, mode, format=format) as fo:
        for kk, dim in enumerate(dimnames):
            if dim not in fo.dimensions:
                fo.createDimension(dim, dimlens[kk])

        for kk, var in enumerate(varnames):
            # Create variable if does not exist
            if var not in fo.variables:
                var1 = fo.createVariable(
                    var, datatype=dtypes[kk], dimensions=vardims[kk]
                )
            else:
                var1 = fo.variables[var]

            # Filling the variable in
            try:
                var1[:] = copy.copy(variables[kk])
            except:
                print(__file__)
                import code
                code.interact(local=dict(locals(), **globals()))
            setattr(var1, "units", units[kk])

            if attributes is not None:
                for key in attributes[kk]:
                    setattr(var1, key, attributes[kk][key])
