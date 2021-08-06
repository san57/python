import datetime
import os
import pytest
import yaml
import shutil


@pytest.fixture()
def dummy_config_fwd(tmpdir):
    """
    Fixture that returns a temporary folder and a configuration file for the
    dummy_forward model.
    """
    current_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    root_dir = os.path.abspath(os.path.join(current_dir, "../../.."))
    data_dir = os.path.abspath(os.path.join(root_dir, "data/"))

    tmpdir_str = tmpdir.strpath
    # tmpdir_str = "/home/chimereges/aberchet/pytest/fwd/"
    # os.system("mkdir -p {}".format(tmpdir_str))

    config = {
        "verbose": 1,
        "logfile": "pycif.logtest",
        "workdir": tmpdir_str,
        "datei": datetime.date(2010, 1, 1),
        "datef": datetime.date(2010, 1, 5),
        "mode": {
            "plugin": {"name": "forward", "version": "std"},
            "perturb_obsvect": True,
            "obserror": 0.01,
        },
        "model": {
            "plugin": {"name": "dummy", "version": "std"},
            "file_pg": os.path.join(
                root_dir, "model_sources/dummy_gauss/Pasquill-Gifford.txt"
            ),
            "chemistry": {"acspecies": {"CH4": None}},
        },
        "measurements": {
            "plugin": {"name": "random", "version": "std"},
            "file_monitor": os.path.join(tmpdir_str, "monitor_reference.nc"),
            "dump_type": "nc",
            "species": {
                "CH4": {
                    "frequency": "1H",
                    "nstations": 50,
                    "duration": "1H",
                    "random_subperiod_shift": True,
                    "zmax": 100,
                }
            },
        },
        "obsvect": {
            "plugin": {"name": "standard", "version": "std"},
            "file_obsvect": os.path.join(tmpdir_str, "monitor_reference.nc"),
            "dump_type": "nc",
            "transform": {
                "timeavg": {
                    "plugin": {"name": "timeavg", "version": "std", "type": "transform"}
                }
            },
        },
        "statevect": {
            "plugin": {"name": "standard", "version": "std"},
            "components": {
                "fluxes": {
                    "parameters": {
                        "CH4": {
                            "plugin": {
                                "name": "dummy",
                                "version": "txt",
                                "type": "fluxes",
                            },
                            "hresol": "hpixels",
                            "type": "physical",
                            "errtype": "max",
                            "err": 1,
                            "period": "1D",
                            "dir": os.path.join(tmpdir_str, "statevect/"),
                            "file": "flx_real.txt",
                            "hcorrelations": {
                                "landsea": False,
                                "dump_hcorr": True,
                                "dircorrel": os.path.join(tmpdir_str, "statevect/"),
                                "sigma": 3000,
                            },
                            "tcorrelations": {"sigma_t": 5},
                            "flx_text": "CIF",
                        }
                    }
                },
                "meteo": {
                    "plugin": {"name": "dummy", "version": "csv", "type": "meteo"},
                    "dir": os.path.join(data_dir, "dummy_gauss/"),
                    "file": "meteo2.csv",
                    "resolution": "1H",
                },
            },
        },
        "domain": {
            "plugin": {"name": "dummy", "version": "std"},
            "xmin": 0,
            "xmax": 25000,
            "nlon": 30,
            "ymin": 0,
            "ymax": 20000,
            "nlat": 15,
        },
    }

    with open(os.path.join(tmpdir_str, "dummy_config.yml"), "w") as outfile:
        yaml.dump(config, outfile)

    # Doing the test
    from pycif.utils.classes.setup import Setup
    result = Setup.run_simu({"def_file": outfile.name})

    yield result
    
    # Copying monitor for later inversion
    shutil.copy(
        "{}/obsvect/monitor.nc".format(tmpdir_str),
        "{}/../".format(tmpdir_str)
    )
    
    # Flushing pycif before next test
    Setup.flushall()
