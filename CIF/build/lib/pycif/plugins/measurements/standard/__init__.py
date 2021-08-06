import shutil

import pandas as pd
from builtins import str

from pycif.utils import path
from pycif.utils.check import info
from pycif.utils.classes.obsparsers import ObsParser
from pycif.utils.datastores import dump
from pycif.utils.datastores.empty import init_empty


def parse_tracers(self, datei, datef, file_monitor="", workdir="", **kwargs):
    """Parses all observation files related to the tracers specified as
    inputs

    Args:
        self (Measurement) : dictionary of tracers as defined in the Yaml
        file
        datei (datetime.datetime): initial date for the inversion window
        datef (datetime.datetime): end date for the inversion window
        file_monitor (str): file with pre-compile observations if exists
        workdir (str): working directory
        logfile (str): path to the log file for verbose instances
        **kwargs (dictionary) : any additional argument that might be useful
                                for extracting observations.
                                Default contains config_dict

    Returns:
        dictionary : dictionary with all observations

    Notes:
        The data that are kept in the datastore are also saved in a
        monit_standard.txt file for debugging
    """

    # Dump type: default is nc
    self.dump_type = getattr(self, "dump_type", "nc")

    # If file_monitor is defined, tries reading it
    if hasattr(self, "file_monitor"):
        file_monitor = self.file_monitor

        try:
            info("Extracting measurements from {}".format(file_monitor))
            return dump.read_datastore(
                file_monitor, dump_type=self.dump_type, **kwargs
            )

        except IOError as e:
            info(e.message)
            info("Could not find a monitor file, reading observations")

        except Exception as e:
            info(e)
            info("Could not read the specified monitor file: {}", file_monitor)
            raise e

    # Otherwise, create the monitor from observations
    if hasattr(self, "workdir"):
        workdir = self.workdir

    file_monitor = workdir + "/obs/monit_standard.nc"

    # If the measurement definition is empty in the Yaml,
    # return an empty datastore
    if not hasattr(self, "species"):
        return init_empty()

    # Loops through tracers if monitor not found
    path.init_dir(workdir + "/obs/")
    shutil.rmtree(file_monitor, ignore_errors=True)

    datastore = {}

    # Merging species datastores into one data
    for spec in self.species.attributes:
        specattr = getattr(self.species, spec)

        if hasattr(specattr, "format") and hasattr(specattr, "provider"):
            info(
                "Extracting measurements for {}"
                " with a single provider {}".format(spec, specattr.provider)
            )

            parser = ObsParser.get_parser(specattr)
            datastore[spec] = parser.parse_multiple_files(spec=spec)
            datastore[spec] = datastore[spec].loc[str(datei): str(datef)]
            continue

        else:
            info(
                "Extracting measurements for {} from multiple providers".format(
                    spec
                )
            )

        # Looping over providers
        dataspec = {}
        for provider in specattr.attributes:
            # Get the observation parser
            providattr = getattr(specattr, provider)
            parser = ObsParser.get_parser(providattr)

            # Reading all files from provider
            dataspec[provider] = parser.parse_multiple_files(spec=spec)

            # Cropping to the inversion window
            dataspec[provider] = dataspec[provider].loc[
                str(datei): str(datef)
            ]

        datastore[spec] = pd.concat(list(dataspec.values()))

    # Grouping species into a single datastore
    datastore = pd.concat(list(datastore.values()))

    # Dumping
    dump.dump_datastore(
        datastore, file_monitor, workdir, dump_type=self.dump_type, **kwargs
    )

    return datastore
