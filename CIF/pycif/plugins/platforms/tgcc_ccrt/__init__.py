import subprocess

# It is necessary to have some measurements and some info about the meteo
# to initialize the observation vector
requirements = {"model": {"any": True, "empty": True}}


def ini_data(plugin, **kwargs):
    """Initializes the platform

    Args:

    """

    # Loading NetCDF module according to the model that will be run
    if (
        plugin.model.plugin.name == "CHIMERE"
        and plugin.model.plugin.version == "std"
    ):
        module_load("load netcdf-fortran")
        module_load("load pnetcdf")
        module_load("load grib")

        # Compile variables
        plugin.model.NETCDFLIB = "${NETCDFFORTRAN_LIBDIR}"
        plugin.model.NETCDFINC = "${NETCDFFORTRAN_INCDIR}"
        plugin.model.GRIBLIB = "${GRIB_LIBDIR}"
        plugin.model.GRIBINC = "${GRIB_INCDIR}"
        plugin.model.LDFLAGS = (
            "-L.  -linitio -ltools -lmodel -ltools -L${NETCDFLIB} -lnetcdff"
        )

    if (
        plugin.model.plugin.name == "LMDZ"
        and plugin.model.plugin.version == "std"
    ):
        # Unload NetCDF
        module_load("unload netcdf")

        # Load NetCDF
        module_load("load netcdf/4")


def module_load(args):
    """Apply the shell command module load/unload to the shell"""

    process = subprocess.Popen(
        "module {}".format(args), shell=True, stdout=subprocess.PIPE
    )
    stdout = process.communicate()
