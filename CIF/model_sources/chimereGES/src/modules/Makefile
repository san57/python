#-*-makefile-*-

# Cuisine et dependances ...

include		../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 

all:	chimere_common.o wholedomain_common.o chimere_consts.o master_message_subs.o message_defs.o

chimere_common.o: chimere_consts.o message_defs.o

master_message_subs.o: \
	chimere_common.o wholedomain_common.o message_defs.o master_message_subs.f90
	$(MPIFC) -o master_message_subs.o $(FFLAGS)  -c master_message_subs.f90

%.o:	%.f90
	$(MPIFC) -o $@ $(FFLAGS)  -c $<

clean:
	$(RM) $(RM_LIST)


