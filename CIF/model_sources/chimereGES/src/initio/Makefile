#-*-makefile-*-
# initio

include		../../Makefile.hdr

LDFLAGS	=	 -L${NETCDFLIB} -lnetcdf -lnetcdff 

RM_LIST =	*.a *.o *.e *.mod core

OBJS	=	checkcfl.o \
		endchimere.o findspec.o \
		iniboun.o inicdf.o inichem.o inichimere.o \
		iniconc.o inidepo.o inidepout.o inidoms.o iniemis.o iniend.o \
		inigeom.o iniout.o iniparam.o iniphys.o \
		iniread.o outprint.o readhour.o renewhour.o \
		write_depo.o

all:		libinitio.a

libinitio.a:	$(OBJS)
	       	$(RM) libinitio.a
		$(AR) libinitio.a  $(OBJS)

%.o :		%.f90
		$(REALFC) -I.. -o $@ $(FFLAGS) -I.. -c $<

clean:
		$(RM) $(RM_LIST)


