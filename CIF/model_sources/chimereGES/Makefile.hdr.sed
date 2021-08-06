#-*-makefile-*-

### This header file is automatically included in the secondary Makefiles.
### Please tune it to your own installation

### Specify where the headers and libraries of your netCDF package reside.
# To avoid trouble, netCDF should have been compiled with the 
# same compiler you use to build CHIMERE
# In most Linux distributions, netCDF has been compiled using g77.
# This may not be compatible with the f90 compilers required for CHIMERE.
#
NETCDFLIB	=	/usr/local/install/netcdf-3.6.2/lib
NETCDFINC	=	/usr/local/install/netcdf-3.6.2/include

### If you want to build the ECMWF meteo interface
# with grib_api (for grib 1 and 2)
GRIBLIB = 	/usr/local/install/grib_api-1.14/lib
GRIBINC =	/usr/local/install/grib_api-1.14/include
F2C	=	g2c


MONITORING = _MONITOR_
ifeq	($(MONITORING),YES)
### Where is your compiler located
### You can get it by issuing the command "which ifort"
REALFC	=	/usr/local/install/scorep-1.1/bin/scorep ifort
### Where is your mpif77 wrapper located
### You can get it by issuing the command "which mpif77"
MF77	=	/usr/local/install/scorep-1.1/bin/scorep /usr/local/bin/mpif77
endif
ifeq	($(MONITORING),NO)
REALFC	=	ifort
MF77	=	mpif77
endif

### Choose your execution mode { PROD | DEVEL }
### PROD is fast, DEVEL allows for more checking and error tracking
MODE	=	_DEBUG_

### If you work with a high resolution grid and many levels, the size of
#   your data segment may be higher than 2 GB. In this case, CHIMERE shall
#   be compiled with special options. If you choose BIGARRAY = YES, then
#   these special options will be selected, at the expense of a slower
#   run-time execution. Please note that the lam and netCDF libraries
#   shall also be built for large addressing space. Refer to the LAM and NETCDF
#   HOWTOs in this directory.
BIGARRAY = NO

### If you use the Fedora Core 4 GNU/Linux distribution, you may
#   experience problems with ifort and Interprocedural Optimisation. 
#   If this is the case, you should disable it.
#   Otherwise just comment out the following line to get the maximum
#   performance of CHIMERE.
#FC4_BUG = -no-ipo
FC4_BUG	= -no-ipo

#########################################################################
### In principle, you should not have to modify too many things below ...

MPIFC	=	export LAMMPIF77=$(REALFC); $(MF77)
FC	=	$(MPIFC)
MPIFLAG	=	MPI

##### IFORT #####
COMPILO	=	FINE
F77=$(FC)
ifeq	($(MODE),DEVEL)
# For debug/development
F77FLAGS1 = -I${NETCDFINC} -I${GRIBINC} -fpe0 -fpp -ftrapuv -g -p -traceback -DIFORT -D$(MPIFLAG) $(FC4_BUG) -fp-model strict -check all
endif
ifeq	($(MODE),PROD)
# for production
F77FLAGS1 = -I${NETCDFINC} -I${GRIBINC} -fpe0 -fpp -O2 -ip -mp1 -prec_div -DIFORT -D$(MPIFLAG) $(FC4_BUG)
endif
ifeq	($(BIGARRAY),YES)
# For data segment > 2GB
F77FLAGS = $(F77FLAGS1) -mcmodel=medium -i-dynamic
else
F77FLAGS = $(F77FLAGS1)
endif
FFLAGS_BIG = $(F77FLAGS) -free
FFLAGS = $(F77FLAGS) -free


# Misc. commands
RM	=	/bin/rm -f
AR	=	/usr/bin/ar r
CPP	=	/usr/bin/cpp
LN	=	/bin/ln -sf
CD	=	cd

.SUFFIXES:
