from types import MethodType

import numpy as np

from pycif.utils.check import info
from pycif.utils.check.errclass import PluginError
from pycif.utils.geometry.areas import calc_areas
from .baseclass import Plugin


class Domain(Plugin):
    @classmethod
    def register_plugin(cls, name, version, module, **kwargs):
        """Register a module for a plugin and version with possibly options

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            module (types.ModuleType): module defining the interface
                between pyCIF and the plugin
            plugin_type (str): type of plugin
            **kwargs (dictionary): default options for module

        """

        super(Domain, cls).register_plugin(
            name, version, module, plugin_type="domain"
        )

    def read_grid(self, *args, **kwargs):
        """Read a grid from an existing file

        Args:
            self (Domain): plugin defining the domain. Should include
            filegrid to be able to read the grid from a file

        Return:
            Grid domain with meshgrids for center lon/lat and corner lon/lat
        """
        raise PluginError("The function read_grid was not defined")

    def create_domain(self, *args, **kwargs):
        """Creates a grid if needed

        Args:
            domain (dictionary): dictionary defining the domain.
        """
        raise PluginError("The function create_domain was not defined")

    def calc_areas(self, *args, **kwargs):
        """Computes the area of each grid cell in your domain."""
        return calc_areas(self, **kwargs)

    def initiate_template(self, plg_type=None):
        """Initializes a Plugin given the template

        Args:
            plg_type (str): the plugin type

        """

        # Initializes the Plugin from the parent method
        module = super(Domain, self).initiate(plg_type="domain")

        # Replacing auxiliary functions
        if hasattr(module, "read_grid"):
            self.read_grid = MethodType(module.read_grid, self)

        if hasattr(module, "create_domain"):
            self.create_domain = MethodType(module.create_domain, self)

        if hasattr(module, "calc_areas"):
            self.calc_areas = MethodType(module.calc_areas, self)

    def ini_data(self, **kwargs):
        """Initializes the domain depending on the model used for the inversion.
        Defines a domain grid from a grid file or a set of parameters if the
        domain was not already defined. Domains are model dependant, so the
        outputs can be different depending
        on the model.

        Args:
            plugin (DomainPlugin): domain definition

        Returns:
            Updates on the fly the domain
        """

        # Read or create a domain
        try:
            # Read domain
            self.read_grid(**kwargs)

        except (IOError, AttributeError):
            # Generate a domain
            info("Couldn't read the domain. Generating it.")
            self.create_domain(**kwargs)

        # Compute areas that can be needed for emissions or diagnostics
        if not hasattr(self, "areas") and getattr(
            self, "compute_areas", False
        ):
            info("Computing areas")
            self.calc_areas(**kwargs)

        # Get side coordinates
        # self.get_sides()

    def get_sides(self):
        """
        Gets the sides of the domain

        :return:
        """

        # Concatenating together the longitudes and latitudes of the sides
        # Orders: West, East, South, North
        self.zlonc_side = np.concatenate(
            [
                self.zlonc[0, :],
                self.zlonc[-1, :],
                self.zlonc[:, 0],
                self.zlonc[:, -1],
            ]
        )[:, np.newaxis]
        self.zlatc_side = np.concatenate(
            [
                self.zlatc[0, :],
                self.zlatc[-1, :],
                self.zlatc[:, 0],
                self.zlatc[:, -1],
            ]
        )[:, np.newaxis]

        self.zlon_side = np.concatenate(
            [
                0.5 * (self.zlonc[0, :-1] + self.zlonc[0, 1:]),
                0.5 * (self.zlonc[-1, :-1] + self.zlonc[-1, 1:]),
                0.5 * (self.zlonc[:-1, 0] + self.zlonc[1:, 0]),
                0.5 * (self.zlonc[:-1, -1] + self.zlonc[1:, -1]),
            ]
        )[:, np.newaxis]

        self.zlat_side = np.concatenate(
            [
                0.5 * (self.zlatc[0, :-1] + self.zlatc[0, 1:]),
                0.5 * (self.zlatc[-1, :-1] + self.zlatc[-1, 1:]),
                0.5 * (self.zlatc[:-1, 0] + self.zlatc[1:, 0]),
                0.5 * (self.zlatc[:-1, -1] + self.zlatc[1:, -1]),
            ]
        )[:, np.newaxis]

        self.nlon_side = self.zlon_side.size
        self.nlat_side = 1
