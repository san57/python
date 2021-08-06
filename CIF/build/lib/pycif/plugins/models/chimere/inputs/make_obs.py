import numpy as np
from builtins import str


def make_obs(self, datastore, runsubdir, mode):

    # If empty datastore, do nothing
    if datastore.size == 0:
        return

    # Otherwise, crop the datastore to active species
    self.dataobs = datastore

    mask = (
        self.dataobs["parameter"]
        .str.upper()
        .isin(self.chemistry.acspecies.attributes)
    )
    data2write = self.dataobs.loc[mask]

    # Include only part of the datastore
    val = "obs_incr" if mode == "adj" else "obs"

    # Write in a txt file
    obs_file = "{}/obs.txt".format(runsubdir)
    with open(obs_file, "w") as f:
        # total number of model's values that are necessary
        nbobs = len(data2write)
        nbdatatot = data2write["dtstep"].sum()

        # write header
        f.write(str(nbobs) + " " + str(nbdatatot) + "\n")

        # write data
        for d in data2write.iterrows():
            ddata = d[1]
            spec = ddata["parameter"].upper()

            # convert level to fortran
            ddata["level"] += 1

            # TODO: deal with altitude
            if np.isnan(ddata["level"]):
                ddata["level"] = 1

            # Otherwise, loop over required time steps
            for dt in np.arange(ddata["dtstep"]):
                nbtstepglo = ddata["tstep_glo"] + dt
                nh = self.nhour[nbtstepglo]
                tstep = self.subtstep[nbtstepglo]
                toprint = "{} {} {:.0f} {:.0f} {:.0f} {}".format(
                    nh - 1, tstep, ddata["i"], ddata["j"], ddata["level"], spec
                )

                # Append increment value in adjoint mode
                if mode == "adj":
                    toprint += " {:23.19e}".format(ddata[val])

                f.write(toprint + "\n")
