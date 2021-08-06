#-*-makefile-*-
# initio

include		../../Makefile.hdr

LDFLAGS	=	 -L${NETCDFLIB} -lnetcdf -lnetcdff 

RM_LIST =	*.a *.o *.e *.mod core

OBJS	=	checkcfl.o \
		endchimere.o findspec.o \
		iniboun.o inicdf.o inichem.o inichimere.o \
		iniconc.o inidepo.o inidoms.o iniemis.o inigeom.o iniphys.o \
		iniparam.o iniend.o inidepout.o \
		iniread.o outprint.o readhour.o renewhour.o iniout.o \
		aendchimere.o ainichimere.o ainiaemis.o ainiaboun.o \
		ainiaconc.o ainiemis.o ainiconc.o \
		ainiend.o ainiout.o ainiphys.o areadhour.o arenewhour.o \
		awrite_adj.o awrite_aconcs.o

all:		libinitio.a

libinitio.a:	$(OBJS)
	       	$(RM) libinitio.a
		$(AR) libinitio.a  $(OBJS)

%.o :		%.f90
		$(REALFC) -I.. -o $@ $(FFLAGS) -I.. -c $<

clean:
		$(RM) $(RM_LIST)


