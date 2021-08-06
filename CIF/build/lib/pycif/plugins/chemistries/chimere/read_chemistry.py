import pandas as pd

from pycif.utils.check import info


def read_chemicalscheme(chemistry, **kwargs):
    """Reads a chemical scheme from existing files

    Args:
        chemistry (pycif.utils.classes.chemistries.Chemistry): the chemical
            scheme

    Return:
        Grid dictionary with characteristics of the chemical scheme

    Notes:

    """

    info("Reading Chemistry")

    workdir = chemistry.workdir
    dirchem_ref = "{}/chemical_scheme/{}/".format(workdir, chemistry.schemeid)

    # ACTIVE SPECIES
    file_chem = "{}/ACTIVE_SPECIES.{}".format(dirchem_ref, chemistry.schemeid)
    acspecies = pd.read_csv(
        file_chem, header=None, sep=" ", usecols=[0, 1], names=["ID", "name"]
    )
    chemistry.acspecies = chemistry.from_dict(
        {s: None for s in acspecies["name"]}
    )
    chemistry.nacspecies = len(acspecies)

    # ANTHROPIC
    file_chem = "{}/ANTHROPIC.{}".format(dirchem_ref, chemistry.schemeid)
    emis_species = pd.read_csv(
        file_chem, header=None, sep=" ", usecols=[0, 1], names=["ID", "name"]
    )
    chemistry.emis_species = chemistry.from_dict(
        {s: None for s in emis_species["name"]}
    )
    chemistry.nemisspec = len(emis_species)

    # BIOGENIC
    file_chem = "{}/BIOGENIC.{}".format(dirchem_ref, chemistry.schemeid)
    bio_species = pd.read_csv(
        file_chem, header=None, sep=" ", usecols=[0, 1], names=["ID", "name"]
    )
    chemistry.bio_species = chemistry.from_dict(
        {s: None for s in bio_species["name"]}
    )
    chemistry.nemisspec_interp = len(bio_species)

    # DEPO_SPEC
    file_chem = "{}/DEPO_SPEC.{}".format(dirchem_ref, chemistry.schemeid)
    dep_species = pd.read_csv(
        file_chem, header=None, sep=" ", usecols=[0], names=["name"]
    )
    chemistry.dep_species = chemistry.from_dict(
        {s: None for s in dep_species["name"]}
    )
    chemistry.ndepspecies = len(dep_species)

    # CHEMISTRY
    with open(dirchem_ref + "CHEMISTRY." + chemistry.schemeid, "r") as fsp:
        ln = fsp.readlines()
        chemistry.nreacs = len(ln)

    with open(dirchem_ref + "FAMILIES." + chemistry.schemeid, "r") as fsp:
        ln = fsp.readlines()
        chemistry.nfamilies = len(ln)

    # TODO: generalize number of prescribed species
    chemistry.nprspecies = 4
