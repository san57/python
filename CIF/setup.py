#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import find_packages, setup
import subprocess
import sys
import os
import stat


def GDAL_with_version():
    try:
        tmp = subprocess.check_output(["gdal-config", "--version"])
        return "gdal==" + str(tmp.decode("utf8").replace("\n", ""))
    except NOT_FOUND_EXCEPTION:
        print("The gdal-config binary was not found")
        raise
    except:
        raise


if sys.version_info[0] < 3:
    NOT_FOUND_EXCEPTION = OSError
    NUMPY = "numpy<1.17"
    CFTIME = "cftime==1.1.1"
    setup_kwargs = {"setup_requires": [NUMPY]}
else:
    NOT_FOUND_EXCEPTION = FileNotFoundError
    NUMPY = "numpy"
    CFTIME = "cftime"
    setup_kwargs = {}

requirements = [
    GDAL_with_version(),
    NUMPY,
    CFTIME,
    'scipy',
    'matplotlib',
    'pandas',
    'netCDF4',
    'pytz',
    'python-dateutil',
    'xarray',
    'pyyaml',
    'pyproj',
    'psutil',
    'Pillow',
    'cfgrib<=0.9.8.1',
    'future'
]

setup(
    name='pyCIF',
    version='0.0.2',
    description=(
        'A Python-based interface to the Community Inversion Framework'
    ),
    author='Antoine Berchet',
    author_email='antoine.berchet@lsce.ipsl.fr',
    url='http://community-inversion.eu/',
    license='CeCILL-B',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[requirements],
    extras_require={
        'dev': [
            'black',
            'Sphinx',
            'sphinx_rtd_theme',
            'sphinxcontrib-napoleon'
        ],
        'test': [
            'coverage',
            'flake8',
            'flake8-html',
            'pytest',
            'pytest-benchmark',
            'pytest-cov',
            'pytest-html',
            'tox'
        ]
    },
    python_requires='>=2.7',
    entry_points={"console_scripts": ['pycif=pycif.__main__:main']},
    zip_safe=False,
    **setup_kwargs
)

# Initializing hooks
if not os.path.isfile(".git/hooks/post-checkout"):
    os.symlink(
        "{}/hooks/post-checkout".format(os.getcwd()),
        "{}/.git/hooks/post-checkout".format(os.getcwd())
    )
    os.chmod(".git/hooks/post-checkout", stat.S_IRWXU)

if not os.path.isfile(".git/hooks/post-merge"):
    os.symlink(
        "{}/hooks/post-merge".format(os.getcwd()),
        "{}/.git/hooks/post-merge".format(os.getcwd())
    )
    os.chmod(".git/hooks/post-merge", stat.S_IRWXU)

# Compiling a first time FLEXPART sources
process = subprocess.Popen(
    "f2py -c mod_flexpart.f90 -m mod_flexpart",
    shell=True,
    stdout=subprocess.PIPE,
    cwd=os.path.join(os.getcwd(), "pycif/plugins/models/flexpart/utils"),
    stderr=subprocess.PIPE,
)
stdout = process.communicate()