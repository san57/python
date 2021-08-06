import os
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
        # Unload NetCDF
        module_load("unload netcdf")

        # Load NetCDF3
        module_load("load netcdf/3")

        # Compile variables
        plugin.model.NETCDFLIB = "/usr/local/install/netcdf-3.6.2/lib"
        plugin.model.NETCDFINC = "/usr/local/install/netcdf-3.6.2/include"
        plugin.model.GRIBLIB = "/usr/local/install/grib_api-1.14/lib"
        plugin.model.GRIBINC = "/usr/local/install/grib_api-1.14/include"

    if (
        plugin.model.plugin.name == "LMDZ"
        and plugin.model.plugin.version == "std"
    ):
        # Unload NetCDF
        module_load("unload netcdf")

        # Load NetCDF3
        module_load("load netcdf/4")


def module_load(args):
    """Apply the shell command module load/unload to the shell"""

    process = subprocess.Popen(
        "/usr/bin/modulecmd tcsh {}".format(args),
        shell=True,
        stdout=subprocess.PIPE,
    )
    stdout = process.communicate()[0].decode("utf-8")

    # If nothing to do, return
    if len(stdout) == 0:
        return

    # Otherwise, loop over all sub-commands 'setenv', 'unsetenv'
    stdout = stdout.split(";")
    for cmd in stdout:
        args = cmd.split()

        if args == []:
            continue

        if args[0] == "setenv":
            value = args[2]
            if args[1] == "PATH":
                value = "{}:{}".format(os.getenv(args[1]), value)

            os.putenv(args[1], value)

        elif args[0] == "unsetenv":
            os.unsetenv(args[1])
