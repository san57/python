import numpy as np
from builtins import str


def make_obs(self, datastore, runsubdir, mode):
    # If empty datastore, do nothing
    if datastore.size == 0:
        return

    self.dataobs = datastore

    col2extract = "obs_incr" if mode == "adj" else "obs"

    # Include only part of the datastore
    data = datastore[
        ["tstep", "dtstep", "i", "j", col2extract, "level"]
    ].values

    # Converts level, tstep, i, j to fortran levels
    data[:, [0, 2, 3, 5]] += 1

    # Assumes that stations with no level are in first level
    # TODO: make it general
    data[np.isnan(data[:, -1]), -1] = 1
    data[data[:, -1] <= 0, -1] = 1
    data[:, -1] /= 100.0

    # Attribute ID to each species
    for k, s in enumerate(self.chemistry.acspecies.attributes):
        mask = (datastore["parameter"].str.lower() == s.lower()).values
        data[mask, -1] += k + 1

    # Write in a binary
    obs_file = "{}/obs.bin".format(runsubdir)
    data.T.tofile(obs_file)

    with open("{}/obs.txt".format(runsubdir), "w") as f:
        f.write(str(len(data)))
