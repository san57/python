#-*-makefile-*-
# tools

include		../../Makefile.hdr

LDFLAGS	=	 -L${NETCDFLIB} -lnetcdf -lnetcdff 

RM_LIST =	*.a *.o *.e core


LIBOBJS	=	calendar.o io.o subs.o setgas.o layave.o intex.o


all:		libtools.a

libtools.a:	$(LIBOBJS)
	       	$(RM) libtools.a
		$(AR) libtools.a  $(LIBOBJS)

%.o :		%.f90
		$(REALFC) -o $@ $(FFLAGS) -I.. -c $<


clean:
		$(RM) $(RM_LIST)


