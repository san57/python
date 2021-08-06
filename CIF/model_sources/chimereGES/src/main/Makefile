#-*-makefile-*-

# Cuisine et dependances ...

include		../../Makefile.hdr

RM_LIST =	*.a *.o *.e *.mod core 


all:	chimere.o integrun.o master_locvalues.o


%.o:	%.f90
	$(MPIFC) -I.. -o $@ $(FFLAGS)  -c $<

clean:
	$(RM) $(RM_LIST)


