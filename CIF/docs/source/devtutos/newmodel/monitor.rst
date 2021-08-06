##################################################
Generate a monitor file compatible with the CIF
##################################################

.. role:: bash(code)
   :language: bash


The CIF compares observations to simulations to carry out data assimilation.
Observations are stored in a so-called :bash:`monitor.nc` file.
Please have a look :doc:`here<../../documentation/monitor>` to know more about :bash:`monitor.nc` files.

We prepare in the following an example of :bash:`monitor.nc`.

    1. Generate an arbitrary :bash:`monitor.nc` with random data inside:
        this can be done by using the Yaml example file available :doc:`here<Yaml_monitor>`;
        you need to change paths in this yaml and can change the extent of the domain you use;
        once modified, run pyCIF with this Yaml:
        :bash:`python -m pycif your_monitor_yaml.yml`

    2. Check the content of the :bash:`monitor.nc` you just created:
        here is a Python template to quick look the content of the :bash:`monitor.nc`:

        .. code:: python

            import matplotlib.pyplot as plt
            from pycif.utils.datastores import dump

            # Reading the monitor file into a Pandas dataframe
            datastore = dump.read_datastore(your_monitor_file)

            # Plot a map with the location of your stations
            sct = plt.scatter(datastore['lon'], datastore['lat'], c=datastore['alt'])

            plt.colorbar(sct)

            plt.show()

    3. Transform the observations of your test case into a pyCIF compatible file:
        if you do not want to start from scratch and rather use observations from a test case,
        you can use the random template to create a real-case :bash:`monitor.nc`;
        at this step, you only need a file gathering information on the raw observations;
        this includes the following columns:
            * date
            * duration
            * station
            * network
            * parameter
            * lon, lat, alt
            * obs
            * obserror

        To produce such a file, you need to fill manually a pandas Datastore compatible with pyCIF from
        your existing observations.
        Details are given :doc:`here<../../documentation/monitor>` about pyCIF datastores.

You now have a file storing all the raw observations you used for your test case.
