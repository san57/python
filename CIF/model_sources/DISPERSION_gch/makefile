#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Definitions de Macros pour Make
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Repertoires :

LOCAL_DIR= ${GCM}
FILTRE   = filtrez
DIM = 3
DIRMAIN=dyn$(DIM)d
RM=rm

# Les differentes librairies pour l"edition des liens:

dyn3d      = $(LIBO)/libdyn3d.a $(LIBO)/lib$(FILTRE).a
dyn2d      = $(LIBO)/libdyn2d.a
dyn1d      = $(LIBO)/libdyn1d.a
L_DYN      = -ldyn$(DIM)d
L_FILTRE   = -l$(FILTRE)
L_PHY = -lphylmd 
L_BIBIO    = -lbibio

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Option de compilation FORTRAN
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPTIMTRU90 = -cpp -O3 -DCPP_PARA -assume byterecl -r8 -DNC_DOUBLE -assume buffered_io #-check all -traceback
INCLUDE = -I$(LIBF)/grid -I$(LIBF)/bibio -I.  -I$(LIBF)/dyn3d  -I/usr/local/include -I$(LIBF)/phylmd -I${LOCAL_DIR}/libo -module $(LIBO) -I $(LIBO)
F90 = mpif90
COMPTRU90 = $(F90) $(OPTIMTRU90) $(INCLUDE) -c
LINK    = $(F90)  -cpp -O -DCPP_PARA -DUSE_VT -align all -assume byterecl -r8 -DNC_DOUBLE -assume buffered_io  -extend_source -I${LOCAL_DIR}/libo -module $(LIBO) -I $(LIBO)
AR      = ar
#OPTION_LINK = -L/usr/local/lib -lnetcdff -lnetcdf -I/usr/local/include
OPTION_LINK = -lnetcdff -lnetcdf

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creation des differents executables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Executables:
# ------------

PROG = dispersion

main : dyn bibio phys 
	cd $(LIBO) ; ls lib*.a ; cd $(GCM) ;\
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/$(DIRMAIN)/$(PROG).F90 -o $(PROG).o ; \
	$(LINK) $(PROG).o -L$(LIBO) $(L_DYN) $(L_FILTRE) $(L_PHY) $(L_BIBIO) $(OPTION_LINK) -o $(LOCAL_DIR)/$(PROG).e ; $(RM) $(PROG).o 

dyn : $(LIBO)/libdyn$(DIM)d.a $(FILTRE)$(DIM)d

phys : $(LIBO)/libphylmd.a

bibio : $(LIBO)/libbibio.a

$(FILTRE)3d : $(LIBO)/lib$(FILTRE).a

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Contenu des differentes bibliotheques
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory bibio
#=======================================================================

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(vampir.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(write_field.o)

.PRECIOUS	: $(LIBO)/libbibio.a


# Compilation des membres de la bibliotheque libbibio.a

$(LIBO)/libbibio.a(vampir.o) : $(LIBF)/bibio/vampir.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/bibio/vampir.F90 ; \
	$(AR) r $(LIBO)/libbibio.a vampir.o ; $(RM) vampir.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(write_field.o) : $(LIBF)/bibio/write_field.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/bibio/write_field.F90 ; \
	$(AR) r $(LIBO)/libbibio.a write_field.o ; $(RM) write_field.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory dyn3d
#=======================================================================

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(coefpoly.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(cray.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(disvert.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(diverg.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dynstate0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(dynrestart.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(exner_hyb.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxhyp.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxy.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxyhyper.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fxysinus.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(fyhyp.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_dyn_fi_ad_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_dyn_fi.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_dyn_fi_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_ecrit_fi.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_fi_dyn_ad_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gr_fi_dyn_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(heavyside.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(iniconst.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(inigeom.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(initial0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(intex.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(intex_tl.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(intex_ad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(readchem.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(readfluxnc_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(massbar.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(pression.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(read_dstoke.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(read_dstoke_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(read_fstoke0.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(read_fstoke0_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(timeloop.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(timeloop_ad.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(tropo.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vlsplt_ad_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vlsplt_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(vlsplt_tl_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(bands.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(gather_tabobs.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(mod_const_para.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(mod_hallo.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(parallel.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(times.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(write_field_p.o)

$(LIBO)/libdyn3d.a :  $(LIBO)/libdyn3d.a(write_func_int.o)

.PRECIOUS	: $(LIBO)/libdyn3d.a


# Compilation des membres de la bibliotheque libdyn3d.a

$(LIBO)/libdyn3d.a(coefpoly.o) : $(LIBF)/dyn3d/coefpoly.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/coefpoly.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a coefpoly.o ; $(RM) coefpoly.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(cray.o) : $(LIBF)/dyn3d/cray.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/cray.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a cray.o ; $(RM) cray.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(disvert.o) : $(LIBF)/dyn3d/disvert.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/disvert.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a disvert.o ; $(RM) disvert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(diverg.o) : $(LIBF)/dyn3d/diverg.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/diverg.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a diverg.o ; $(RM) diverg.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dynstate0.o) : $(LIBF)/dyn3d/dynstate0.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/logic.h \
$(LIBO)/libphylmd.a(modchem.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/dynstate0.F90; \
	$(AR) r $(LIBO)/libdyn3d.a dynstate0.o ; $(RM) dynstate0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(dynrestart.o) : $(LIBF)/dyn3d/dynrestart.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/dynrestart.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a dynrestart.o ; $(RM) dynrestart.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(exner_hyb.o) : $(LIBF)/dyn3d/exner_hyb.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/serre.h \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/exner_hyb.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a exner_hyb.o ; $(RM) exner_hyb.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxhyp.o) : $(LIBF)/dyn3d/fxhyp.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/fxhyp.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a fxhyp.o ; $(RM) fxhyp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxy.o) : $(LIBF)/dyn3d/fxy.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxy_new.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/fxy.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a fxy.o ; $(RM) fxy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxyhyper.o) : $(LIBF)/dyn3d/fxyhyper.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/fxyhyper.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a fxyhyper.o ; $(RM) fxyhyper.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fxysinus.o) : $(LIBF)/dyn3d/fxysinus.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/grid/fxy_sin.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/fxysinus.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a fxysinus.o ; $(RM) fxysinus.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(fyhyp.o) : $(LIBF)/dyn3d/fyhyp.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/fyhyp.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a fyhyp.o ; $(RM) fyhyp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_dyn_fi_ad_p.o) : $(LIBF)/dyn3d/gr_dyn_fi_ad_p.F90 \
$(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o) \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libphylmd.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/gr_dyn_fi_ad_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_dyn_fi_ad_p.o ; $(RM) gr_dyn_fi_ad_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_dyn_fi.o) : $(LIBF)/dyn3d/gr_dyn_fi.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/gr_dyn_fi.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_dyn_fi.o ; $(RM) gr_dyn_fi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_dyn_fi_p.o) : $(LIBF)/dyn3d/gr_dyn_fi_p.F90 \
$(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o) \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/gr_dyn_fi_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_dyn_fi_p.o ; $(RM) gr_dyn_fi_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_ecrit_fi.o) : $(LIBF)/dyn3d/gr_ecrit_fi.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/gr_ecrit_fi.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_ecrit_fi.o ; $(RM) gr_ecrit_fi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_fi_dyn_ad_p.o) : $(LIBF)/dyn3d/gr_fi_dyn_ad_p.F90 \
$(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o) \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/gr_fi_dyn_ad_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_fi_dyn_ad_p.o ; $(RM) gr_fi_dyn_ad_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gr_fi_dyn_p.o) : $(LIBF)/dyn3d/gr_fi_dyn_p.F90 \
$(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o) \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libphylmd.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/gr_fi_dyn_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a gr_fi_dyn_p.o ; $(RM) gr_fi_dyn_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(heavyside.o) : $(LIBF)/dyn3d/heavyside.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/heavyside.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a heavyside.o ; $(RM) heavyside.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(iniconst.o) : $(LIBF)/dyn3d/iniconst.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/iniconst.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a iniconst.o ; $(RM) iniconst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(inigeom.o) : $(LIBF)/dyn3d/inigeom.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comdissnew.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/inigeom.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a inigeom.o ; $(RM) inigeom.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(initial0.o) : $(LIBF)/dyn3d/initial0.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/initial0.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a initial0.o ; $(RM) initial0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(readchem.o) : $(LIBF)/dyn3d/readchem.F90 \
$(LIBF)/dyn3d/ajout.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libphylmd.a(modchem.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/readchem.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a readchem.o ; $(RM) readchem.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(intex.o) : $(LIBF)/dyn3d/intex.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/intex.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a intex.o ; $(RM) intex.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(intex_tl.o) : $(LIBF)/dyn3d/intex_tl.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/intex_tl.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a intex_tl.o ; $(RM) intex_tl.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(intex_ad.o) : $(LIBF)/dyn3d/intex_ad.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/intex_ad.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a intex_ad.o ; $(RM) intex_ad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(readfluxnc_p.o) : $(LIBF)/dyn3d/readfluxnc_p.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/tracstoke.h \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libbibio.a(vampir.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/readfluxnc_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a readfluxnc_p.o ; $(RM) readfluxnc_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(massbar.o) : $(LIBF)/dyn3d/massbar.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/massbar.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a massbar.o ; $(RM) massbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(pression.o) : $(LIBF)/dyn3d/pression.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/pression.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a pression.o ; $(RM) pression.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(read_dstoke.o) : $(LIBF)/dyn3d/read_dstoke.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/read_dstoke.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a read_dstoke.o ; $(RM) read_dstoke.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(read_dstoke_p.o) : $(LIBF)/dyn3d/read_dstoke_p.F90 \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/read_dstoke_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a read_dstoke_p.o ; $(RM) read_dstoke_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(read_fstoke0.o) : $(LIBF)/dyn3d/read_fstoke0.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/read_fstoke0.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a read_fstoke0.o ; $(RM) read_fstoke0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(read_fstoke0_p.o) : $(LIBF)/dyn3d/read_fstoke0_p.F90 \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/read_fstoke0_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a read_fstoke0_p.o ; $(RM) read_fstoke0_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(timeloop.o) : $(LIBF)/dyn3d/timeloop.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/tracstoke.h \
$(LIBF)/dyn3d/ajout.h \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libdyn3d.a(mod_hallo.o) \
$(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libdyn3d.a(times.o) \
$(LIBO)/libdyn3d.a(write_field_p.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/timeloop.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a timeloop.o ; $(RM) timeloop.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(timeloop_ad.o) : $(LIBF)/dyn3d/timeloop_ad.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/tracstoke.h \
$(LIBF)/dyn3d/ajout.h \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libdyn3d.a(mod_hallo.o) \
$(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libdyn3d.a(times.o) \
$(LIBO)/libdyn3d.a(write_field_p.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/timeloop_ad.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a timeloop_ad.o ; $(RM) timeloop_ad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(tropo.o) : $(LIBF)/dyn3d/tropo.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/tropo.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a tropo.o ; $(RM) tropo.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vlsplt_ad_p.o) : $(LIBF)/dyn3d/vlsplt_ad_p.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libdyn3d.a(mod_hallo.o) \
$(LIBO)/libdyn3d.a(write_field_p.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/vlsplt_ad_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a vlsplt_ad_p.o ; $(RM) vlsplt_ad_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vlsplt_p.o) : $(LIBF)/dyn3d/vlsplt_p.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libdyn3d.a(mod_hallo.o) \
$(LIBO)/libdyn3d.a(write_field_p.o) \
$(LIBO)/libdyn3d.a(times.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/vlsplt_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a vlsplt_p.o ; $(RM) vlsplt_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(vlsplt_tl_p.o) : $(LIBF)/dyn3d/vlsplt_tl_p.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libdyn3d.a(mod_hallo.o) \
$(LIBO)/libdyn3d.a(write_field_p.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/vlsplt_tl_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a vlsplt_tl_p.o ; $(RM) vlsplt_tl_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(bands.o) : $(LIBF)/dyn3d/bands.F90 \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libdyn3d.a(times.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/bands.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a bands.o ; $(RM) bands.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(gather_tabobs.o) : $(LIBF)/dyn3d/gather_tabobs.F90 \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/gather_tabobs.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a gather_tabobs.o ; $(RM) gather_tabobs.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(mod_const_para.o) : $(LIBF)/dyn3d/mod_const_para.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/mod_const_para.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a mod_const_para.o ; $(RM) mod_const_para.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(mod_hallo.o) : $(LIBF)/dyn3d/mod_hallo.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/mod_hallo.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a mod_hallo.o ; $(RM) mod_hallo.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(mod_interface_dyn_phys.o) : $(LIBF)/dyn3d/mod_interface_dyn_phys.F90 \
$(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/mod_interface_dyn_phys.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a mod_interface_dyn_phys.o ; $(RM) mod_interface_dyn_phys.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(parallel.o) : $(LIBF)/dyn3d/parallel.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libdyn3d.a(mod_const_para.o) \
$(LIBO)/libbibio.a(vampir.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/parallel.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a parallel.o ; $(RM) parallel.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(times.o) : $(LIBF)/dyn3d/times.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/times.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a times.o ; $(RM) times.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(write_field_p.o) : $(LIBF)/dyn3d/write_field_p.F90 \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libbibio.a(write_field.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/write_field_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a write_field_p.o ; $(RM) write_field_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d.a(write_func_int.o) : $(LIBF)/dyn3d/write_func_int.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libdyn3d.a(parallel.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/dyn3d/write_func_int.F90 ; \
	$(AR) r $(LIBO)/libdyn3d.a write_func_int.o ; $(RM) write_func_int.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory filtrez
#=======================================================================

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(acc.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(eigen_sort.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(filtreg.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(inifgn.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(inifilr.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(jacobi.o)

.PRECIOUS	: $(LIBO)/libfiltrez.a


# Compilation des membres de la bibliotheque libfiltrez.a

$(LIBO)/libfiltrez.a(acc.o) : $(LIBF)/filtrez/acc.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/filtrez/acc.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a acc.o ; $(RM) acc.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(eigen_sort.o) : $(LIBF)/filtrez/eigen_sort.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/filtrez/eigen_sort.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a eigen_sort.o ; $(RM) eigen_sort.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(filtreg.o) : $(LIBF)/filtrez/filtreg.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/filtrez/parafilt.h \
$(LIBF)/filtrez/coefils.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/filtrez/filtreg.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a filtreg.o ; $(RM) filtreg.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(inifgn.o) : $(LIBF)/filtrez/inifgn.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/filtrez/coefils.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/filtrez/inifgn.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a inifgn.o ; $(RM) inifgn.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(inifilr.o) : $(LIBF)/filtrez/inifilr.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/filtrez/parafilt.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/filtrez/coefils.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/serre.h
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/filtrez/inifilr.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a inifilr.o ; $(RM) inifilr.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(jacobi.o) : $(LIBF)/filtrez/jacobi.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/filtrez/jacobi.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a jacobi.o ; $(RM) jacobi.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phylmd
#=======================================================================

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(cltrac.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(cvltr.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(thermcell_dq.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(nflxtr.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(phytrac.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(read_pstoke0.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(readsdepvel.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(readsflux.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(chemistry.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(rates.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(init_chem.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(modchem.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(suphec.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(comgeomphy.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(dimphy.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(init_phys_lmdz.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(mod_grid_phy_lmdz.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_data.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_transfert.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(mod_phys_lmdz_omp_data.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(mod_phys_lmdz_omp_transfert.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(mod_phys_lmdz_para.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(mod_phys_lmdz_transfert_para.o)

$(LIBO)/libphylmd.a :  $(LIBO)/libphylmd.a(write_field_phy.o)

.PRECIOUS	: $(LIBO)/libphylmd.a


# Compilation des membres de la bibliotheque libphylmd.a

$(LIBO)/libphylmd.a(cltrac.o) : $(LIBF)/phylmd/cltrac.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBO)/libphylmd.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/cltrac.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a cltrac.o ; $(RM) cltrac.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(nflxtr.o) : $(LIBF)/phylmd/nflxtr.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libphylmd.a(write_field_phy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/nflxtr.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a nflxtr.o ; $(RM) nflxtr.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(cvltr.o) : $(LIBF)/phylmd/cvltr.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libphylmd.a(write_field_phy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/cvltr.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a cvltr.o ; $(RM) cvltr.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(thermcell_dq.o) : $(LIBF)/phylmd/thermcell_dq.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libphylmd.a(write_field_phy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/thermcell_dq.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a thermcell_dq.o ; $(RM) thermcell_dq.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(phytrac.o) : $(LIBF)/phylmd/phytrac.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBF)/phylmd/indicesol.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libphylmd.a(write_field_phy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/phytrac.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a phytrac.o ; $(RM) phytrac.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(read_pstoke0.o) : $(LIBF)/phylmd/read_pstoke0.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/phylmd/indicesol.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/read_pstoke0.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a read_pstoke0.o ; $(RM) read_pstoke0.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(readsdepvel.o) : $(LIBF)/phylmd/readsdepvel.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBO)/libphylmd.a(modchem.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/readsdepvel.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a readsdepvel.o ; $(RM) readsdepvel.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(readsflux.o) : $(LIBF)/phylmd/readsflux.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phylmd/dimphy.h \
$(LIBO)/libphylmd.a(modchem.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/readsflux.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a readsflux.o ; $(RM) readsflux.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(chemistry.o) : $(LIBF)/phylmd/chemistry.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libphylmd.a(modchem.o) \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libdyn3d.a(parallel.o) \
$(LIBO)/libphylmd.a(write_field_phy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/chemistry.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a chemistry.o ; $(RM) chemistry.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(rates.o) : $(LIBF)/phylmd/rates.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBO)/libphylmd.a(modchem.o) \
$(LIBO)/libphylmd.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/rates.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a rates.o ; $(RM) rates.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(modchem.o) : $(LIBF)/phylmd/modchem.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/modchem.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a modchem.o ; $(RM) modchem.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(init_chem.o) : $(LIBF)/phylmd/init_chem.F90 \
$(LIBO)/libphylmd.a(modchem.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/init_chem.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a init_chem.o ; $(RM) init_chem.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(suphec.o) : $(LIBF)/phylmd/suphec.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/suphec.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a suphec.o ; $(RM) suphec.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(comgeomphy.o) : $(LIBF)/phylmd/comgeomphy.F90 \
$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/comgeomphy.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a comgeomphy.o ; $(RM) comgeomphy.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(dimphy.o) : $(LIBF)/phylmd/dimphy.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/dimphy.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a dimphy.o ; $(RM) dimphy.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(init_phys_lmdz.o) : $(LIBF)/phylmd/init_phys_lmdz.F90 \
$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphylmd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/init_phys_lmdz.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a init_phys_lmdz.o ; $(RM) init_phys_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(mod_grid_phy_lmdz.o) : $(LIBF)/phylmd/mod_grid_phy_lmdz.F90
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/mod_grid_phy_lmdz.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a mod_grid_phy_lmdz.o ; $(RM) mod_grid_phy_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_data.o) : $(LIBF)/phylmd/mod_phys_lmdz_mpi_data.F90 \
$(LIBO)/libdyn3d.a(mod_const_para.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/mod_phys_lmdz_mpi_data.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a mod_phys_lmdz_mpi_data.o ; $(RM) mod_phys_lmdz_mpi_data.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_transfert.o) : $(LIBF)/phylmd/mod_phys_lmdz_mpi_transfert.F90 \
$(LIBO)/libphylmd.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_data.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/mod_phys_lmdz_mpi_transfert.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a mod_phys_lmdz_mpi_transfert.o ; $(RM) mod_phys_lmdz_mpi_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(mod_phys_lmdz_omp_data.o) : $(LIBF)/phylmd/mod_phys_lmdz_omp_data.F90 \
$(LIBO)/libphylmd.a(dimphy.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/mod_phys_lmdz_omp_data.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a mod_phys_lmdz_omp_data.o ; $(RM) mod_phys_lmdz_omp_data.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(mod_phys_lmdz_omp_transfert.o) : $(LIBF)/phylmd/mod_phys_lmdz_omp_transfert.F90 \
$(LIBO)/libphylmd.a(mod_phys_lmdz_omp_data.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/mod_phys_lmdz_omp_transfert.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a mod_phys_lmdz_omp_transfert.o ; $(RM) mod_phys_lmdz_omp_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o) : $(LIBF)/phylmd/mod_phys_lmdz_para.F90 \
$(LIBO)/libphylmd.a(mod_phys_lmdz_transfert_para.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphylmd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/mod_phys_lmdz_para.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a mod_phys_lmdz_para.o ; $(RM) mod_phys_lmdz_para.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(mod_phys_lmdz_transfert_para.o) : $(LIBF)/phylmd/mod_phys_lmdz_transfert_para.F90 \
$(LIBO)/libphylmd.a(mod_phys_lmdz_mpi_transfert.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_omp_transfert.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/mod_phys_lmdz_transfert_para.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a mod_phys_lmdz_transfert_para.o ; $(RM) mod_phys_lmdz_transfert_para.o ; \
	cd $(GCM)

$(LIBO)/libphylmd.a(write_field_phy.o) : $(LIBF)/phylmd/write_field_phy.F90 \
$(LIBO)/libphylmd.a(dimphy.o) \
$(LIBO)/libphylmd.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphylmd.a(mod_grid_phy_lmdz.o)
	cd $(LOCAL_DIR); \
	$(COMPTRU90) $(LIBF)/phylmd/write_field_phy.F90 ; \
	$(AR) r $(LIBO)/libphylmd.a write_field_phy.o ; $(RM) write_field_phy.o ; \
	cd $(GCM)


clean : 
	rm -rf $(LIBO)/*
	rm -f $(LOCAL_DIR)/$(PROG).e
