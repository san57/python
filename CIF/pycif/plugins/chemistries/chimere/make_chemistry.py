import os

from pycif.utils.path import init_dir
from .utils import create_mandchem, create_optchem


def create_chemicalscheme(self):
    """
    Read the mandatory chemistry files, create the optional ones and set the
    species attributes
    consistently with the data retrieved
    """

    # Common variables
    workdir = self.workdir
    mecachim = self.schemeid
    dirchem_ref = self.dirchem_ref

    # Initializes number of species and reactions
    self.nacspecies = (
        len(self.acspecies.attributes) if hasattr(self, "acspecies") else 0
    )
    self.nemisspec = (
        len(self.emis_species.attributes)
        if hasattr(self, "emis_species")
        else 0
    )
    self.nemisspec_interp = (
        len(self.emis_species_interp.attributes)
        if hasattr(self, "emis_species_interp")
        else 0
    )
    self.nprspecies = (
        len(self.prescrconcs.attributes) if hasattr(self, "prescrconcs") else 0
    )
    self.nprodspecies = (
        len(self.prodloss3d.attributes) if hasattr(self, "prodloss3d") else 0
    )
    self.ndepspecies = (
        len(self.deposition.attributes) if hasattr(self, "deposition") else 0
    )
    self.nreacs = (
        len(self.reactions.attributes) if hasattr(self, "reactions") else 0
    )
    self.nfamilies = (
        len(self.families.attributes) if hasattr(self, "families") else 0
    )

    # Cleaning the target directory
    os.system("rm -rf  {}".format(dirchem_ref))
    init_dir(dirchem_ref)

    # List of files
    finf = "{}/chemical_scheme.nml".format(dirchem_ref)
    filer = "{}/REACTIONS.{}".format(dirchem_ref, mecachim)
    fileps = "{}/PRESCRIBED_SPECIES.{}".format(dirchem_ref, mecachim)
    filepl = "{}/PRODLOSS_SPECIES.{}".format(dirchem_ref, mecachim)
    filedp = "{}/DEPO_SPECIES.{}".format(dirchem_ref, mecachim)
    filea = "{}/ANTHROPIC.{}".format(dirchem_ref, mecachim)
    fileb = "{}/BIOGENIC.{}".format(dirchem_ref, mecachim)

    mandatory_files = [filer, fileps, filepl, filedp]
    create_mandchem(self, mandatory_files)

    files = "{}/STOICHIOMETRY.{}".format(dirchem_ref, mecachim)
    filec = "{}/CHEMISTRY.{}".format(dirchem_ref, mecachim)
    filerr = "{}/REACTION_RATES.{}".format(dirchem_ref, mecachim)
    filej = "{}/PHOTO_RATES.{}".format(dirchem_ref, mecachim)
    filef = "{}/FAMILIES.{}".format(dirchem_ref, mecachim)
    fileals = "{}/ALL_SPECIES.{}".format(dirchem_ref, mecachim)
    fileas = "{}/ACTIVE_SPECIES.{}".format(dirchem_ref, mecachim)

    nallqmax, nphoto_rates = create_optchem(self, filer, fileps)

    # Create chemical_scheme.nml namelist
    os.system('echo "&args" > {}'.format(finf))
    os.system("echo \"fnacspec = '{}'\" >> {}".format(fileas, finf))
    os.system("echo \"fnallspec = '{}'\" >> {}".format(fileals, finf))
    os.system("echo \"fnprescr = '{}'\" >> {}".format(fileps, finf))
    os.system("echo \"fnprodl = '{}'\" >> {}".format(filepl, finf))
    os.system("echo \"fndep = '{}'\" >> {}".format(filedp, finf))
    os.system("echo \"fnchem = '{}'\" >> {}".format(filec, finf))
    os.system("echo \"fnstoi = '{}'\" >> {}".format(files, finf))
    os.system("echo \"fnrates = '{}'\" >> {}".format(filerr, finf))
    os.system("echo \"fnjrates = '{}'\" >> {}".format(filej, finf))
    os.system("echo \"fnfamilies = '{}'\" >> {}".format(filef, finf))
    os.system('echo "iqmax = {}" >> {}'.format(self.nspecies, finf))
    os.system('echo "iallqmax = {}" >> {}'.format(nallqmax, finf))
    os.system('echo "iprescrmax = {}" >> {}'.format(self.nprspecies, finf))
    os.system('echo "iprodmax = {}" >> {}'.format(self.nprodspecies, finf))
    os.system('echo "idepmax = {}" >> {}'.format(self.ndepspecies, finf))
    os.system('echo "nreac = {}" >> {}'.format(self.nreacs, finf))
    os.system('echo "ijratesmax = {}" >> {}'.format(nphoto_rates, finf))
    os.system('echo "/" >> {}'.format(finf))
