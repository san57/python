import os
import subprocess
from shutil import copytree, ignore_patterns, rmtree, copy
from pycif.utils.check import info

def compile(self):
    comp_dir = "{}/model".format(self.workdir)

    # Copying the executables
    if not getattr(self, "auto-recompile", False):
        source = "{}/src/fwdchimere.e".format(self.direxec)
        copy(source, comp_dir)

        source = "{}/src_tl/tlchimere.e".format(self.direxec)
        copy(source, comp_dir)

        source = "{}/src_ad/achimere.e".format(self.direxec)
        copy(source, comp_dir)
        return

    # Otherwise, re-compile
    # Copying sources locally; overwrite folder if exists
    if os.path.isdir("{}/sources".format(comp_dir)):
        rmtree("{}/sources".format(comp_dir))

    copytree(
        self.direxec,
        "{}/sources".format(comp_dir),
        ignore=ignore_patterns("*.a"),
    )

    # Modifying LDFLAGS if needed
    if hasattr(self, "LDFLAGS"):
        for fld in ["src", "src_tl", "src_ad"]:
            old_lines = open(
                "{}/sources/{}/Makefile_chimere".format(comp_dir, fld), "r"
            ).readlines()
            old_lines = [ln.rstrip() for ln in old_lines]

            for k, ln in enumerate(old_lines):
                if ln[:7] == "LDFLAGS":
                    old_lines[k] = "LDFLAGS =\t{}".format(self.LDFLAGS)

            with open(
                "{}/sources/{}/Makefile_chimere".format(comp_dir, fld), "w"
            ) as f:
                f.write("\n".join(old_lines))

    # Modifying Makefile header if needed
    old_lines = open(
        "{}/sources/Makefile.hdr.sed".format(comp_dir), "r"
    ).readlines()
    old_lines = [ln.rstrip() for ln in old_lines]

    for k, ln in enumerate(old_lines):
        if ln[:9] == "NETCDFLIB" and hasattr(self, "NETCDFLIB"):
            old_lines[k] = "NETCDFLIB\t=\t{}".format(self.NETCDFLIB)
        if ln[:9] == "NETCDFINC" and hasattr(self, "NETCDFINC"):
            old_lines[k] = "NETCDFINC\t=\t{}".format(self.NETCDFINC)
        if ln[:7] == "GRIBLIB" and hasattr(self, "GRIBLIB"):
            old_lines[k] = "GRIBLIB\t=\t{}".format(self.GRIBLIB)
        if ln[:7] == "GRIBINC" and hasattr(self, "GRIBINC"):
            old_lines[k] = "GRIBINC\t=\t{}".format(self.GRIBINC)

    with open("{}/sources/Makefile.hdr.sed".format(comp_dir), "w") as f:
        f.write("\n".join(old_lines))

    # Now compiling
    comp_mode = "N" if getattr(self, "compile-mode") == "PROD" else "Y"
    comp_clean = "Y" if getattr(self, "compile-clean") else "N"

    for mode, pref, suff in zip(
        ["A", "L", "D"], ["a", "tl", "fwd"], ["_ad", "_tl", ""]
    ):
        process = subprocess.Popen(
            "./compile-chimere {} {} {}".format(mode, comp_mode, comp_clean),
            shell=True,
            stdout=subprocess.PIPE,
            cwd="{}/sources/".format(comp_dir),
            stderr=subprocess.PIPE,
        )
        stdout = process.communicate()
        
        source = "{}/sources/src{}/{}chimere.e".format(comp_dir, suff, pref)
        copy(source, comp_dir)
