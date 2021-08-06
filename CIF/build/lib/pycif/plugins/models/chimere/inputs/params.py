def make_nml(self, runsubdir, sdc, mode):
    # Default values
    nvegtype = getattr(self, "nvegtype", 16)

    # Fills in chimere.nml
    domain = self.domain
    chemistry = self.chemistry
    dirchem = chemistry.dirchem_ref
    schemeid = chemistry.schemeid
    nho = self.nho

    # NML string to write
    with open("{}/chimere.nml".format(runsubdir), "w") as f:
        f.write("&args\n")

        f.write("version   = '{}'\n".format(mode))
        f.write("domain = {}\n".format(domain.domid))
        f.write("nphour_ref = {}\n".format(self.nphour_ref))
        f.write("ichemstep = {}\n".format(self.ichemstep))
        f.write("nzonal_domain = {}\n".format(domain.nlon))
        f.write("nmerid_domain = {}\n".format(domain.nlat))
        f.write("nvert_raw = {}\n".format(domain.nlev))
        f.write("nlevemis = {}\n".format(self.nlevemis))
        f.write("nspec = {}\n".format(chemistry.nacspecies))
        f.write("nemisa= {}\n".format(chemistry.nemisspec))
        f.write("nemisb= {}\n".format(chemistry.nemisspec_interp))
        f.write("ndep= {}\n".format(chemistry.ndepspecies))
        f.write("nreac= {}\n".format(chemistry.nreacs))
        f.write("nfam = {}\n".format(chemistry.nfamilies))
        f.write("nspresc = {}\n".format(chemistry.nprspecies))
        f.write("nreactamax = {}\n".format(chemistry.nreactamax))
        f.write("ntemps = {}\n".format(chemistry.ntemps))
        f.write("ntabmax = {}\n".format(chemistry.ntabmax))
        f.write("nlevphotmax = {}\n".format(chemistry.nlevphotmax))
        f.write("ntabuzenmax = {}\n".format(chemistry.ntabuzenmax))
        f.write("nphotmax = {}\n".format(chemistry.nphotmax))
        f.write("ivsurf = {}\n".format(domain.emissublayer + 1))
        f.write("ideepconv = {}\n".format(self.ideepconv))
        f.write("idatestart = {}\n".format(sdc))
        f.write("nhourrun = {}\n".format(nho))
        f.write("nzdoms = {}\n".format(self.nzdoms))
        f.write("nmdoms = {}\n".format(self.nmdoms))

        f.write("fninit = '{}/INI_CONCS.0.nc'\n".format(runsubdir))

        f.write(
            "fnoutspec = '{}/OUTPUT_SPECIES.{}'\n".format(dirchem, schemeid)
        )
        f.write("fnspec = '{}/ACTIVE_SPECIES.{}'\n".format(dirchem, schemeid))
        f.write("fnchem = '{}/CHEMISTRY.{}'\n".format(dirchem, schemeid))
        f.write("fnstoi = '{}/STOICHIOMETRY.{}'\n".format(dirchem, schemeid))
        f.write("fnrates = '{}/REACTION_RATES.{}'\n".format(dirchem, schemeid))
        f.write("fnfamilies = '{}/FAMILIES.{}'\n".format(dirchem, schemeid))
        f.write(
            "fnphot = '{}/PHOTO_PARAMETERS.{}'\n".format(dirchem, schemeid)
        )
        f.write("fnanthro = '{}/ANTHROPIC.{}'\n".format(dirchem, schemeid))
        f.write("fnbiogen = '{}/BIOGENIC.{}'\n".format(dirchem, schemeid))

        f.write("fnemisa = '{}/AEMISSIONS.nc'\n".format(runsubdir))
        f.write("fnemisb = '{}/BEMISSIONS.nc'\n".format(runsubdir))
        f.write("fnbounconc = '{}/BOUN_CONCS.nc'\n".format(runsubdir))
        f.write("fnmeteo = '{}/METEO.nc'\n".format(runsubdir))

        f.write("fndepoespe = '{}/DEPO_SPEC.{}'\n".format(dirchem, schemeid))
        f.write("fndepopars = '{}/DEPO_PARS.{}'\n".format(dirchem, schemeid))
        f.write("fnwetd = '{}/WETD_SPEC.{}'\n".format(dirchem, schemeid))

        f.write(
            "fnlanduse = '{}/LANDUSE/LANDUSE_{}'\n".format(
                domain.repgrid, domain.domid
            )
        )

        f.write("fnout = '{}/out.nc'\n".format(runsubdir))
        f.write("fnothers = '{}/par.nc'\n".format(runsubdir))
        if mode != "adj":
            f.write("fnconcs = '{}/end.nc'\n".format(runsubdir))
        else:
            f.write("fnconcs = '{}/aend.nc'\n".format(runsubdir))

        f.write("fndepos = '{}/dep.nc'\n".format(runsubdir))

        f.write("fniCOOcorn = '{}'\n".format(domain.corners))
        f.write("usechemistry = {}\n".format(self.usechemistry))
        f.write("usedepos = {}\n".format(self.usedepos))
        f.write("useemissions = {}\n".format(self.useemissions))
        f.write("usetransmix = {}\n".format(self.usetransmix))
        f.write("usewetdepos = {}\n".format(self.usewetdepos))
        f.write("useabsclipconc = {}\n".format(self.useabsclipconc))
        f.write("dryairout = {}\n".format(self.dryairout))
        f.write("nvegtype = {}\n".format(self.nvegtype))
        f.write("nlduse = {}\n".format(self.nlduse))
        f.write("nparammax = {}\n".format(self.nparammax))
        f.write("clipconc = {}\n".format(self.clipconc))
        f.write("ntyperate = {}\n".format(self.ntyperate))
        f.write("ihoursu = {}\n".format(self.ihoursu))
        f.write("nivout = {}\n".format(self.nivout))
        f.write("nsaveconcs = {}\n".format(self.nsaveconcs))
        f.write("nsavedepos = {}\n".format(self.nsavedepos))

        if not self.optemisb:
            f.write("optemisb = 0\n")
        else:
            f.write("optemisb = 1\n")

        f.write("nitgs = {}\n".format(self.nitgs))
        f.write("nitgssu = {}\n".format(self.nitgssu))

        if mode != "tl":
            f.write("fnemisaincr = 'dummy'\n")
            f.write("fnemisbincr = 'dummy'\n")
            f.write("fnbounconcincr = 'dummy'\n")
            f.write("fninitincr = 'dummy'\n")

        if mode == "tl":
            f.write(
                "fnemisaincr = '{}/AEMISSIONS.increment.nc'\n".format(
                    runsubdir
                )
            )
            f.write(
                "fnemisbincr = '{}/BEMISSIONS.increment.nc'\n".format(
                    runsubdir
                )
            )
            f.write(
                "fnbounconcincr = '{}/BOUN_CONCS.increment.nc'\n".format(
                    runsubdir
                )
            )
            f.write(
                "fninitincr = '{}/INI_CONCS.0.increment.nc'\n".format(
                    runsubdir
                )
            )

        f.write("afnout = '{}/aout.'\n".format(runsubdir))

        ainit = 1 if (mode == "adj" and hasattr(self, "chain")) else 0
        f.write("iopainit = {}\n".format(ainit))

        netcdfoutput = 1 if self.dumpncoutput else 0
        f.write("NetCDF_output = {}\n".format(netcdfoutput))

        netcdfpar = 1 if self.dumpncpar else 0
        f.write("NetCDF_parout = {}\n".format(netcdfpar))

        netcdf_outtype = 1 if self.dumpnctype == "double" else 0
        f.write("NetCDF_outtype = {}\n".format(netcdf_outtype))

        f.write("fwd = {}\n".format(1 if mode == "fwd" else 0))
        f.write("tl = {}\n".format(1 if mode == "tl" else 0))
        f.write("ad = {}\n".format(1 if mode == "adj" else 0))
        # pulse hour Elise Potier
        f.write("hpulse = {}\n".format(self.hpulse))

        #        f.write("idatestart = "+str(datei)+'\n')
        #        f.write("nhourrun = "+str(nho)+'\n')
        #        f.write("nzo  ="+str(domain.nzonal))
        #        f.write("nme  ="+str(domain.nmerid))
        #        nsho = int(nzo*nme/2+nzo)
        #        f.write("nshow="+str(nsho))
        f.write("/\n")
