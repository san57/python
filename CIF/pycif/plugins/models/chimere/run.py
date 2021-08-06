import datetime
import os
import shutil
import subprocess

from builtins import str

from pycif.utils.check import info


def run(self, runsubdir, mode, workdir, nbproc=1, do_simu=True, **kwargs):
    """Run the CHIMERE model in forward mode

    Args:
        self: the model Plugin
        runsubdir (str): working directory for the current run
        mode (str): forward or backward
        workdir (str): pyCIF working directory
        do_simu (bool): re-run or not existing simulation

    """

    if not do_simu:
        if mode in ["fwd", "tl"]:
            # Keeps the running directory in memory for later adjoint
            # simulations
            self.adj_refdir = "{}/../".format(runsubdir)
        return

    # Cleaning the directory
    shutil.rmtree("{}/mod.txt".format(runsubdir), ignore_errors=True)

    # Number of processors
    nbproc = str(self.nzdoms * self.nmdoms + 1)

    # Running the model
    info("Running sub-simulation in {}".format(runsubdir))
    with open("{}/chimere.log".format(runsubdir), "w") as log:
        process = subprocess.Popen(
            [
                "time",
                "-p",
                "mpirun",
                "-np",
                nbproc,
                "-mca",
                "btl",
                "self,sm,tcp",
                "chimere.e",
            ],
            cwd=runsubdir,
            stdout=log,
            stderr=subprocess.PIPE,
        )
        _, err = process.communicate()

        if err != "" and not os.path.isfile("{}/all_good".format(runsubdir)):
            info("Exception in CHIMERE")
            info(str(err, "utf-8"))
            raise Exception(
                "CHIMERE did not run properly in {}".format(runsubdir)
            )

    # Process here initial conditions for the next simulation
    date = datetime.datetime.strptime(
        os.path.basename(runsubdir), "%Y-%m-%d_%H-%M"
    )
    nho = self.nho
    if mode in ["fwd", "tl"]:
        os.system(
            date.strftime(
                "mv -f {}/end.nc {}/../chain/end.%Y%m%d%H.{}.nc".format(
                    runsubdir, runsubdir, nho
                )
            )
        )
        # Keeps the running directory in memory for later adjoint simulations
        self.adj_refdir = "{}/../".format(runsubdir)

    if mode == "adj":
        os.system(
            date.strftime(
                "mv -f {}/aend.nc {}/../chain/aend.%Y%m%d%H.{}.nc".format(
                    runsubdir, runsubdir, nho
                )
            )
        )
