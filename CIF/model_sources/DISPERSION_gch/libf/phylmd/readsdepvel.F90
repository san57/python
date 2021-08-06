SUBROUTINE readsdepvel(nom,filename,depvel)

USE variables
USE CONSTANTS

IMPLICIT NONE

INCLUDE "netcdf.inc"
INCLUDE "dimensions.h"
INCLUDE "paramet.h"
INCLUDE "dimphy.h"

!---------------------------------------------------------
! Reading deposition velocity
!---------------------------------------------------------

CHARACTER(LEN=*), INTENT(in) :: filename
CHARACTER(len=*), INTENT(in) :: nom
      
!-------------------------------------------------------------
!  Local variables
!-------------------------------------------------------------

INTEGER :: iret
REAL :: rcode
INTEGER :: ncid
REAL, DIMENSION(klon,llm) :: depvel
REAL, DIMENSION(iim,jjm+1,llm) :: tempo
INTEGER, DIMENSION(4) :: start, count	
INTEGER :: iddep

PRINT *, 'Read data from file ', filename
! ... Open the file
ncid = NCOPN(filename,NCNOWRIT,rcode)
PRINT*,'Variable to read=','Dep_'//nom
iddep = NCVID(ncid,'Dep_'//nom,rcode)
start(1)=1
start(2)=1
start(3)=1
start(4)=1
count(1)=iim
count(2)=jjm+1
count(3)=12 ! Assuming deposition file with 12 months
count(4)=1
iret = nf_get_vara_double(ncid, iddep, start, count, tempo)
CALL check_err(iret)
CALL gr_ecrit_fi(llm,klon,iim,jjm+1,tempo,depvel)
! ... Close the file
iret = nf_close(ncid)
CALL check_err(iret)
        
RETURN
END SUBROUTINE readsdepvel
