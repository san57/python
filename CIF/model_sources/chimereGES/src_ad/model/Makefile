#-*-makefile-*-
# model

include		../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

OBJS_CORE	= chemistry.o achemistry.o deposition.o depvel.o worker.o emissions.o \
		htransport.o ainiworker.o locvalues.o  \
		mixing.o photorates.o physics2.o physloc.o ppm.o prep_outprint.o prodloss.o rates.o \
		twostep_mod.o transmix.o vanleer.o  vertav.o vtransport.o \
		wdeposition.o zenith.o aworker.o atwostep_mod.o alocvalues.o aprodloss.o \
		awdeposition.o atransmix.o aemissions.o adeposition.o avanleer.o \
		worker.o worker_message_subs.o worker_mpi_finalize.o \
		awrite_concs_hour.o awrite_concs_phys.o aphysics2.o aphysloc.o


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

twostep_mod.o:		worker_common.o worker_message_subs.o

atwostep_mod.o:		worker_message_subs.o

worker.o:		twostep_mod.o

aworker.o:		atwostep_mod.o




# common rules for other object targets

%.o:			%.f90 $(OBJS)
			$(MPIFC) -I.. -I. -o $@ $(FFLAGS) -c $<

clean:
	$(RM) $(RM_LIST)


