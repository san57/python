#-*-makefile-*-
# initio

include		../../Makefile.hdr

LDFLAGS	=	 -L${NETCDFLIB} -lnetcdf -lnetcdff 

RM_LIST =	*.a *.o *.e *.mod core

OBJS	=	checkcfl.o \
		endchimere_tl.o findspec.o \
		iniboun_tl.o inicdf.o inicdf_tl.o inichem.o inichimere_tl.o \
		iniconc_tl.o inidepo.o inidepout.o inidoms.o iniemis_tl.o iniend_tl.o \
		inigeom.o iniout.o iniparam.o iniphys.o \
		iniread.o outprint_tl.o readhour_tl.o renewhour_tl.o \
		write_depo.o

all:		libinitio.a

libinitio.a:	$(OBJS)
	       	$(RM) libinitio.a
		$(AR) libinitio.a  $(OBJS)

%.o :		%.f90
		$(REALFC) -I.. -o $@ $(FFLAGS) -I.. -c $<

clean:
		$(RM) $(RM_LIST)


