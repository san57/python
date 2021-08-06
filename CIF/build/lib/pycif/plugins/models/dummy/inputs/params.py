import calendar
import datetime

from builtins import str
from builtins import zip


def make_rundef(filedef, datei, physic, runsubdir, lmdz_yearref=1979):
    # Reads the original run.def
    with open(filedef, "r") as f:
        data = f.read()

    # Days since 1st of January of the same year
    daydebjan = 1 + (datei - datetime.datetime(datei.year, 1, 1)).days

    # Number of days in the month
    nday = calendar.monthrange(datei.year, datei.month)[1]

    physic = "T" if physic else "F"

    # Replacing parameters
    pattern_init = ["$DD", "$YYYY", "$NBD", "${nday[${MO}-1]}", "$PHYSIC"]
    pattern_final = [daydebjan, lmdz_yearref, nday, nday, physic]

    for pi, pf in zip(pattern_init, pattern_final):
        data = data.replace(str(pi), str(pf))

    # Fills run.def
    with open("{}/run.def".format(runsubdir), "w") as f:
        f.write(data)


def make_totinput(self, runsubdir, datei, mode, footprint="F"):
    """Makes a totinput file, storing the main simulation parameters:
    - number of effective tracers
    - using SACS (T) or not (F)
    - spliting of dynamic timestep
    - spliting of physical timestep
    - read start file (T) or not (F)
    - forward (T) or backward (F)
    - output diagnostics
    - output wfunc or not (T or F)
    - if footprint (to change number of days 28 to 30..)
    - ndayloc, number of days in the month
    - convOH = T if vmr or F if molec/cm3

    """

    with open("{}/totinput".format(runsubdir), "w") as f:
        nbtr = self.nspecies
        sacs = "T" if hasattr(self, "chemistry") else "F"
        start = (
            "T"
            if mode in ["fwd", "tl"]
            or (mode == "adj" and hasattr(self, "chain"))
            else "F"
        )
        start = "T"
        forward = "T" if mode in ["fwd", "tl"] else "F"
        diag = 0
        wfunc = "F"
        convOH = (
            "T"
            if getattr(
                getattr(
                    getattr(
                        getattr(self, "chemistry", None), "prescrconcs", None
                    ),
                    "OH",
                    None,
                ),
                "convOH",
                True,
            )
            else "F"
        )

        # TODO: periodflux is hard-coded at the moment
        periodflux = 1

        # Number of days in the month
        nday = calendar.monthrange(datei.year, datei.month)[1]
        month = datei.month
        year = datei.year

        for var in [
            nbtr,
            sacs,
            self.domain.dsplit,
            self.domain.psplit,
            start,
            forward,
            diag,
            wfunc,
            footprint,
            nday,
            convOH,
            periodflux,
            nday,
            month,
            year,
        ]:
            f.write(str(var) + "\n")


def make_species(self, runsubdir):
    with open("{}/ACTIVE_SPECIES".format(runsubdir), "w") as f:
        # Two possible ways of defining species in the Yaml:
        # If one species, one can simply give it to species; then,
        # the argument self_ids must be specified at the model level
        # Else, see the regular definition in the wiki
        if type(self.species) == str:
            f.write("{} {}\n".format(self.species.upper(), self.restart_ids))
            self.nspecies = 1
        else:
            self.nspecies = len(self.species.attributes)
            for spec in self.species.attributes:
                i = getattr(self.species, spec).restart_id
                f.write("{} {}\n".format(spec.upper(), i))

    if hasattr(self, "chemistry"):
        if hasattr(self.chemistry, "prescrconcs"):
            with open("{}/PRESCRIBED_SPECIES".format(runsubdir), "w") as f:
                for spec in self.chemistry.prescrconcs.attributes:
                    f.write(spec.upper() + "\n")

        if hasattr(self.chemistry, "deposition"):
            with open("{}/DEPOSITED_SPECIES".format(runsubdir), "w") as f:
                for spec in self.chemistry.deposition.attributes:
                    f.write(spec.upper() + "\n")

        if hasattr(self.chemistry, "prodloss3d"):
            with open("{}/PRODLOSS_SPECIES".format(runsubdir), "w") as f:
                for spec in self.chemistry.prodloss3d.attributes:
                    f.write(spec.upper() + "\n")
