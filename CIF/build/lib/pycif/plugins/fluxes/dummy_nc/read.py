from pycif.utils.check import info


def read(
    self,
    name,
    tracdir,
    tracfile,
    varnames,
    dates,
    interpol_flx=False,
    **kwargs
):
    """Get fluxes from pre-computed fluxes and load them into a pycif
    variables

    Args:
        self: the model Plugin
        name: the name of the component
        tracdir, tracfile: flux directory and file format
        dates: list of dates to extract
        interpol_flx (bool): if True, interpolates fluxes at time t from
        values of surrounding available files

    """

    info("READING NetCDF dummy fluxes")

    raise Exception
