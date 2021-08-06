SUBROUTINE readsflux(filename)

USE variables
USE SPECIES_NAME
USE CONSTANTS

IMPLICIT NONE

INCLUDE 'netcdf.inc'
INCLUDE "dimensions.h"
INCLUDE "paramet.h"
INCLUDE "dimphy.h"

!-------------------------------------------------------------
!       ... Dummy arguments
!-------------------------------------------------------------
CHARACTER(LEN=*), INTENT(in) :: filename
!-------------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------------
INTEGER :: iret
INTEGER :: error
INTEGER :: ncid
INTEGER :: varid
INTEGER :: time   ! ndays ou nmonths

!  LOGICAL :: notentered = .true.
PRINT *, 'Read data from file ', filename

! ... Open the file
iret = nf_open(filename, 0, ncid)
PRINT*,'iret=',iret
CALL check_err(iret)

!---------------------------------------------------------
! Lecture des flux des traceurs
!---------------------------------------------------------

! Get the vector lengths

iret = nf_inq_dimid(ncid, 'time', varid)
CALL check_err(iret)
iret = nf_inq_dimlen(ncid, varid, time_len)
CALL check_err(iret)
time = time_len

iret = nf_inq_dimid(ncid, 'vector', varid)
CALL check_err(iret)
iret = nf_inq_dimlen(ncid, varid, vect)
CALL check_err(iret)

!  Allocations
!TODO : deal with it
!ALLOCATE(flx_CH2O(vect,time_len), STAT=error)
!IF (error /= 0) STOP 'Space requested not possible'
!if (id_CH2O <= iqmax.and.id_CH2O > 0) then
!iret = nf_inq_varid(ncid, 'flx_ch2o', varid)
!CALL check_err(iret)
!iret = nf_get_var_double(ncid, varid, flx_CH2O)
!CALL check_err(iret)
!endif

! ... Close the file
iret = nf_close(ncid)
CALL check_err(iret)

END SUBROUTINE readsflux
