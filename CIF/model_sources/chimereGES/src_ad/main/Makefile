#-*-makefile-*-

# Cuisine et dependances ...

include		../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


all:	chimere.o aintegrun.o dirintegrun.o master_locvalues.o amaster_alocvalues.o


%.o:	%.f90
	$(MPIFC) -I.. -o $@ $(FFLAGS)  -c $<

clean:
	$(RM) $(RM_LIST)


