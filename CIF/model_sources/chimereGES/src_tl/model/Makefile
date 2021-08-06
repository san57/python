#-*-makefile-*-
# model

include		../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


OBJS_CORE	= chemistry.o chemistry_tl.o\
		deposition.o deposition_tl.o depvel.o emissions.o emissions_tl.o searchinterv.o \
		htransport.o iniworker_tl.o locvalues.o  \
		mixing.o photorates.o physics2_tl.o physloc_tl.o \
		ppm.o prep_outprint_tl.o prodloss.o prodloss_tl.o rates.o \
		twostep_mod_tl.o transmix.o transmix_tl.o vanleer.o vanleer_tl.o  vertav.o \
		vtransport.o wdeposition.o wdeposition_tl.o zenith.o worker_tl.o worker_message_subs.o \
		worker_mpi_finalize.o


OBJS		= worker_common.o worker_message_subs.o

all:		libmodel.a


libmodel.a:	$(OBJS_CORE)
			$(RM) libmodel.a
		$(AR) libmodel.a $(OBJS_CORE) $(OBJS)

# specific rules
worker_common.o:	worker_common.f90
			$(MPIFC) -o worker_common.o $(FFLAGS) -I.. -c  worker_common.f90


worker_message_subs.o:		worker_message_subs.f90 worker_common.o
				$(MPIFC) -o worker_message_subs.o $(FFLAGS) -I..  -c worker_message_subs.f90

twostep_mod_tl.o:		worker_common.o worker_message_subs.o

worker_tl.o:		twostep_mod_tl.o






# common rules for other object targets

%.o:			%.f90 $(OBJS)
			$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -c $<

clean:
	$(RM) $(RM_LIST)


