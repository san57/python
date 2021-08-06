import xarray as xr
from netCDF4 import Dataset


def read(
    self,
    name,
    tracdir,
    tracfile,
    varnames,
    dates,
    interpol_flx=False,
    comp_type=None,
    model=None,
    tracer=None,
    **kwargs
):
    """Get fluxes from pre-computed fluxes and load them into a pyCIF
    variables

    Args:
        self: the fluxes Plugin
        name: the name of the component
        tracdir, tracfile: flux directory and file format
        dates: list of dates to extract
        interpol_flx (bool): if True, interpolates fluxes at time t from
        values of surrounding available files

    """

    ic_file = min(dates).strftime("{}/{}".format(tracdir, tracfile))
    with Dataset(ic_file, "r") as f:
        if hasattr(tracer, "restart_id"):
            spec_id = "q{:02d}".format(tracer.restart_id)
        else:
            spec_id = "q{:02d}".format(
                getattr(model.chemistry.acspecies, name).restart_id
            )

        data = f.variables[spec_id][:]

    xmod = xr.DataArray(
        data, coords={"time": [min(dates)]}, dims=("time", "lev", "lat", "lon")
    )

    return xmod
