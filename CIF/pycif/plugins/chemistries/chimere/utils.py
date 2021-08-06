import os
import re
from os.path import exists, getsize

import numpy as np
import pandas as pd
from builtins import range
from builtins import str


def create_mandchem(chemistry, mandatory_files):
    """
    Create the mandatory chemistry files using the chemistry config yml file
    """
    for mfile in mandatory_files:
        os.system("touch {}".format(mfile))

    # Chemical reactions
    with open(mandatory_files[0], "w") as f:
        if chemistry.nreacs > 0:
            for attr in chemistry.reactions.attributes:
                reac = getattr(chemistry.reactions, attr)
                f.write(reac + "\n")

    if hasattr(chemistry, "prescrconcs"):
        with open(mandatory_files[1], "w") as f:
            for attr in chemistry.prescrconcs.attributes:
                f.write(attr + "\n")

    if hasattr(chemistry, "prodloss3d"):
        with open(mandatory_files[2], "w") as f:
            for attr in chemistry.prodloss3d.attributes:
                f.write(attr + "\n")

    if hasattr(chemistry, "deposition"):
        with open(mandatory_files[3], "w") as f:
            for attr in chemistry.deposition.attributes:
                f.write(attr + "\n")


def create_optchem(chemistry, filer, fileps):
    """
    Create the chemistry files ALL_SPECIES, CHEMISTRY, REACTION_RATES,
                               STOICHIOMETRY and ACTIVE_SPECIES
    At the moment, FAMILIES is set empty

    ALL_SPECIES is created

    Args:
        chemistry (pycif.utils.classes.chemistries.Chemistry): chemical scheme
        filer (str) : path to the file with reactions
        fileps (str): path to the file with prescribed species

    """

    workdir = chemistry.workdir
    dirchem_ref = chemistry.dirchem_ref
    mecachim = chemistry.schemeid

    files = "{}/STOICHIOMETRY.{}".format(dirchem_ref, mecachim)
    filec = "{}/CHEMISTRY.{}".format(dirchem_ref, mecachim)
    filerr = "{}/REACTION_RATES.{}".format(dirchem_ref, mecachim)
    filej = "{}/PHOTO_RATES.{}".format(dirchem_ref, mecachim)
    filef = "{}/FAMILIES.{}".format(dirchem_ref, mecachim)
    fileals = "{}/ALL_SPECIES.{}".format(dirchem_ref, mecachim)
    fileas = "{}/ACTIVE_SPECIES.{}".format(dirchem_ref, mecachim)
    fileanth = "{}/ANTHROPIC.{}".format(dirchem_ref, mecachim)
    filebio = "{}/BIOGENIC.{}".format(dirchem_ref, mecachim)

    os.system("rm -f {} {} {} {}".format(files, filec, filerr, filef))

    # Read reactions
    df_reac = pd.read_csv(
        filer,
        index_col=False,
        header=None,
        comment="#",
        engine="python",
        delim_whitespace=True,
    )

    # Read prescribed species
    prescribed_species = np.array([])
    if exists(fileps) and getsize(fileps) > 0:
        df_prescr = pd.read_csv(
            fileps,
            index_col=False,
            header=None,
            comment="#",
            engine="python",
            delim_whitespace=True,
        )
        prescribed_species = df_prescr[0]

    nonactive = np.append(prescribed_species, ["O2", "X", None])

    # Create stoechiometry and chemistry
    parts = df_reac[0].str.split("->", n=2, expand=True)
    left_hand_sides = parts[0]
    right_hand_sides = parts[1]
    losses = left_hand_sides.str.split("+", expand=True)
    nlosses = left_hand_sides.str.split("+").apply(len)
    losses.insert(0, -1, nlosses)
    prods = right_hand_sides.str.split("+", expand=True)
    nprods = pd.Series(np.sum(~pd.isnull(prods).values, axis=1))
    prods.insert(0, -1, nprods)

    prods_spec = np.copy(prods.values)
    stoichiometry = np.zeros((1, 6))
    for (i, j), value in np.ndenumerate(prods.values[:, 1:]):
        try:
            prods_list = prods.values[i, j + 1].split("*")
            if len(prods_list) > 1:
                arr = np.array(
                    [
                        prods_list[1],
                        prods_list[0],
                        prods_list[0],
                        prods_list[0],
                        prods_list[0],
                        i + 1,
                    ]
                )
                stoichiometry = np.append(
                    stoichiometry, arr[np.newaxis, ...], axis=0
                )

            if prods_list[-1] not in nonactive:
                prods_spec[i, j + 1] = prods_list[-1]
            else:
                prods_spec[i, j + 1] = None
                prods_spec[i, 0] -= 1
        except AttributeError:
            pass

    file_chemistry = np.append(losses.values, prods_spec, axis=1)
    file_chemistry = pd.DataFrame(file_chemistry)
    stoichiometry = pd.DataFrame(stoichiometry[1:])

    # Create reactions
    reactions = df_reac[1].values
    reactions_rates, nphoto_rates = read_react(reactions, filej)
    reactions_rates = pd.DataFrame(reactions_rates)

    # Active species
    chemistry.nspecies = len(chemistry.acspecies.attributes)

    # Create active_species (output_species) and all_species
    output_species = np.array(chemistry.acspecies.attributes)
    all_species = np.append(output_species, prescribed_species, axis=0)
    type_spec = np.array([["type"]])
    acinfos = np.array([0, 0])
    type_spec = np.broadcast_to(type_spec, (all_species.shape[0], 1))
    acinfos = np.broadcast_to(acinfos, (output_species.shape[0], 2))
    all_species = np.append(all_species[:, np.newaxis], type_spec, axis=1)
    output_species = np.append(output_species[:, np.newaxis], acinfos, axis=1)
    for i in range(output_species.shape[0]):
        outspec = getattr(chemistry.acspecies, output_species[i, 0])
        output_species[i, 1] = str(getattr(outspec, "restart_id"))
        output_species[i, 2] = str(getattr(outspec, "mass"))

    for i in range(all_species.shape[0]):
        if i < chemistry.nspecies:
            all_species[i, 1] = "ac"
        else:
            all_species[i, 1] = "pr"

    all_species = pd.DataFrame(all_species)
    output_species = pd.DataFrame(output_species)

    # Create files to create
    stoichiometry.to_csv(files, sep=" ", header=False, index=False)
    file_chemistry.to_csv(filec, sep=" ", header=False, index=False)
    reactions_rates.to_csv(filerr, sep=" ", header=False, index=False)
    all_species.to_csv(fileals, sep=" ", header=False, index=False)
    output_species.to_csv(fileas, sep=" ", header=False, index=False)
    output_species.to_csv(fileanth, sep=" ", header=False, index=False)
    output_species.to_csv(filebio, sep=" ", header=False, index=False)

    return all_species.shape[0], nphoto_rates


def read_react(reactions, filej):
    """
    Parse the types of reactions in the REACTIONS file and extract the data
    """
    nreacts = reactions.shape[0]
    reactions_rates = np.zeros((nreacts, 10))
    reactions_rates.fill(None)
    photo_rates = np.zeros((1, 2))
    nphoto_rates = 0

    for i, value in np.ndenumerate(reactions):
        # Type 1 : Constant rate
        if re.search(r"k=", value):
            val_list = re.split(r"=", value)
            reactions_rates[i, 0] = 1
            reactions_rates[i, 1] = val_list[1]

        # Type 2 : Simplified Arrhenius
        if re.search(r"k\(T\)=Aexp\(-B\/T\)", value):
            val_list = re.split(r"=|,", value)
            reactions_rates[i, 0] = 2
            reactions_rates[i, 1] = val_list[3]
            reactions_rates[i, 2] = val_list[5]

        # Type 3 : Complete Arrhenius
        if re.search(r"k\(T\)=Aexp\(-B\/T\)\(300\/T\)\*\*N", value):
            val_list = re.split(r"=|,", value)
            reactions_rates[i, 0] = 3
            reactions_rates[i, 1] = val_list[3]
            reactions_rates[i, 2] = val_list[5]
            reactions_rates[i, 3] = val_list[7]

        # Type 4 : Relative pressure
        if re.search(r"k\(P\)=A\(B\+C\*P/Pref\)", value):
            val_list = re.split(r"=|,", value)
            reactions_rates[i, 0] = 4
            reactions_rates[i, 1] = val_list[3]
            reactions_rates[i, 2] = val_list[5]
            reactions_rates[i, 3] = val_list[7]

        # Type 5 : Simple Photolysis
        if re.search(r"J=", value):
            val_list = re.split(r"=|,", value)
            reactions_rates[i, 0] = 5
            reactions_rates[i, 1] = val_list[1]
            photo_rates = np.append(
                photo_rates, [[str(i[0] + 1), "j" + val_list[1]]], axis=0
            )
            nphoto_rates += 1

    photo_rates = pd.DataFrame(photo_rates[1:])
    photo_rates.to_csv(filej, sep=" ", header=False, index=False)

    return reactions_rates, nphoto_rates
