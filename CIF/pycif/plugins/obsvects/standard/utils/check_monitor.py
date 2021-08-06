from builtins import str

from pycif.utils.check import info


def check_monitor(self):
    """Check the consistency between the observation datastore and the model
    configuration set up"""

    datastore = self.datastore

    # For old monitor with no nc_attributes, do nothing
    if not hasattr(datastore, "nc_attributes"):
        info(
            "Cannot check the datastore from a previous version."
            "Please be careful with the use of it"
        )
        return True, True, True, False

    # Otherwise, check what part of the monitor is to be re-computed
    nc_attributes = datastore.nc_attributes

    ok_hcoord = (
        nc_attributes.get("domain nlat", None) == str(self.model.domain.nlat)
    ) and (
        nc_attributes.get("domain nlon", None) == str(self.model.domain.nlon)
    )

    ok_vcoord = True

    ok_tstep = (
        nc_attributes.get("datei", None)
        == self.datei.strftime("%d-%m-%Y %H:%M:%S")
    ) and (
        nc_attributes.get("datef", None)
        == self.datef.strftime("%d-%m-%Y %H:%M:%S")
    )

    allcorrec = ok_hcoord and ok_tstep and ok_vcoord

    return allcorrec, ok_hcoord, ok_vcoord, not ok_tstep
