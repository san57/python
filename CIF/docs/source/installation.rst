##############
How to install
##############

.. role:: bash(code)
   :language: bash

.. contents:: Contents
    :local:


*********************
Requirements
*********************
Before installing pycif, the following libraries and packages must be installed
on your server:

    - python >=3.6:
        pycif is designed for python 3.6 or above; compatibility with python 2.7 has been lost and is not maintained anymore;

    - git:
        can be installed following instructions from `git-scm.com <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`__
    - geographical and binary tools:
        pycif is based on the following geographical and binary libraries:
            - `GEOS <https://trac.osgeo.org/geos/>`__
            - `GDAL <https://gdal.org/>`__
            - `PROJ <https://proj.org/>`__
            - `hdf5 <https://www.hdfgroup.org/solutions/hdf5/>`__
            - `netCDF4 <https://www.unidata.ucar.edu/software/netcdf/>`__
    - basic libraries and tools (most likely already installed on most systems):
        - `GCC <https://gcc.gnu.org/wiki/InstallingGCC>`__
        - `freetype <https://www.freetype.org/>`__
        - `libpng <http://www.libpng.org/pub/png/libpng.html>`__
        - `OpenBLAS <https://www.openblas.net/>`__
        - `psutils <https://psutil.readthedocs.io/en/latest/>`__
        - `py-yaml <https://pyyaml.org/wiki/PyYAMLDocumentation>`__ or python-yaml (depending on the distribution)
        - `linux-headers <https://wiki.gentoo.org/wiki/Linux-headers>`__
        - ethtool

The following python packages must be installed before getting pycif:

    - `numpy <https://www.numpy.org/>`__
    - `scipy <https://www.scipy.org/>`__
    - `pandas <http://pandas.pydata.org/>`__
    - `xarray <http://xarray.pydata.org/en/stable/>`__
    - `matplotlib <https://matplotlib.org/>`__
    - `pyproj <https://pyproj4.github.io/pyproj/stable/>`__
    - `setuptools <https://setuptools.readthedocs.io/>`__
    - `cython <http://cython.org/>`__
    - `shapely <https://shapely.readthedocs.io/en/stable/manual.html>`__
    - `psutil <https://pypi.org/project/psutil/>`__
    - `netCDF4 <http://unidata.github.io/netcdf4-python/netCDF4/index.html>`__
    - `cftime <https://pypi.org/project/cftime/>`__
    - `PyYAML <https://pyyaml.org/wiki/PyYAMLDocumentation>`__
    - `dateutil <https://dateutil.readthedocs.io/en/stable/>`__
    - `tz <https://pypi.org/project/tz/>`__
    - `Pillow <https://pillow.readthedocs.io/en/stable/installation.html>`__
    - `future <https://docs.python.org/3/library/__future__.html>`__
    - `cfgrib <https://github.com/ecmwf/cfgrib>`__

****************************
Using conda for requirements
****************************

`conda <https://docs.conda.io/en/latest/>`__ enables an easy handling of package installations and virtual environments.
If you are using conda, please find below some instructions to install required packages:

.. code-block:: bash

    # Packages installed in Conda base environment
    conda activate base
    conda install --yes numpy scipy matplotlib basemap h5py netcdf4 pyhdf gdal

    # Make a pyCIF environment
    conda deactivate
    conda create --clone base --name cif
    conda activate cif

    # Install remaining packages for pyCIF
    conda install --yes openblas pyyaml pandas xarray cython shapely psutil pillow future cfgrib

    # Check if the packages are installed
    conda list "(gdal|openblas|pyyaml|pandas|xarray|cython|shapely|psutil|pillow|future|cfgrib)"

You can also find :doc:`here <pycif_conda_specs>` a conda virtual environment file to automatically create a virtual environment.
Copy the content of the link in a text file and use the command:

.. code-block:: bash

    conda create --name <env> --file <the file>

Don't forget to activate your environment before installing pyCIF below:

.. code-block:: bash

    conda activate cif

****************
Getting the code
****************

The main server hosting the CIF is: `git.nilu.no/VERIFY/CIF <https://git.nilu.no/VERIFY/CIF>`__.
Anyone can clone sources and start using them.
At the top right of the git-page there’s a blue button with the word “Clone”.
If you click on that, you’ll see two options, one is “Clone with SSH”, the
other “Clone with HTTPS”.
Both options give you an address which you can use in a git command:

.. code-block:: bash

    git clone GIT/ADDRESS /where/you/want/the/sources

Please be aware that the SSH address works only for users registered at the NILU server.
Contributors to the CIF can ask for a login by sending a mail at `help@community-inversion.eu <help@community-inversion.eu>`__

For other users, please use the HTTPS address.

To get updates from the master version of the CIF, use the following colomn:

.. code-block:: bash

    git pull origin master

For further git commands: `try.github.io <https://try.github.io/>`__

****************
Installing pycif
****************

pyCIF is a python module that can be installed to your python environment.

You can install the CIF as a third-party package.
It will add it to the list of available packages exactly like when using :bash:`conda` or :bash:`pip`:

.. code-block:: bash

    cd /where/you/have/the/sources
    python setup.py install

Please be aware that the command will try to modify some system files, hence requiring super-user authorization.
If you have admin rights on your server, you may use the :bash:`sudo` prefix.

Python also offers the possibility of installing thrid-party packages on the user environment only, which do not requires admin rights:

.. code-block:: bash

    cd /where/you/have/the/sources
    python setup.py install --user

Despite all our efforts to have a distributed version of pycif as stable and robust as possible, thus
if you plan to do any modification
(to fix bugs, implement your own plugins, understand the code better by printing intermediate states, etc.),
use the following command:

.. code-block:: bash

    cd /where/you/have/the/sources
    python setup.py develop --user

Thus, any modification you do on the code will be directly accessible to python scripts using pycif.

*************************
Contributing to pycif
*************************

pycif is a research participative project. Any bug-fixes, plugin developments or general contribution to the present documentation are very welcome.
Guidelines for some typical developments are given :doc:`here <devtutos/developers-tutorials>`.

In any case, contributions should follow as closely as possible `PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ standards.

When doing developments, it is recommended to create a new branch to the git structure.
To do so, please follow the steps:

.. code-block:: bash

    cd /where/you/have/the/sources
    git checkout -b [name_of_your_new_branch]

This create a new local branch from the master.
You can then commit your modifications to this branch.

Once you are satisfied by your new version, you can push it to the main server:

.. code-block:: bash

    git push origin [name_of_your_new_branch]

Please note that this last step requires a login to the NILU server.
If you don't have one yet, you can send a mail to `help@community-inversion.eu <help@community-inversion.eu>`__

To propagate your changes to the master branch, you need to submit a merge request to the git server.
Please visit `this page <https://git.nilu.no/VERIFY/CIF/merge_requests>`__ to merge your branch.
