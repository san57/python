! ==========================================================================
! LMDZT - Atmospheric tracer transport, forward, tangent-linear and adjoint
!             for use within the PYVAR inversion system
!
! Copyright Laboratoire des Sciences du Climat et de l'Environnement (LSCE)
!           Unite mixte CEA-CNRS-UVSQ
! Initial version from the LMDZ.3.3 model, developed by IPSL
!
! Code manager:
! Frederic Chevallier, LSCE, frederic.chevallier@lsce.ipsl.fr
!
! This software is governed by the CeCILL  license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and,  more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
! ==================================================================

SUBROUTINE read_dstoke(irec,zdtvr,ziadvtrac,ziadvtrac2)
  
  
  IMPLICIT NONE

 include "netcdf.inc"
 include "dimensions.h"
 include "paramet.h"
 include "comgeom.h"
 include "comvert.h"
  
  
  INTEGER,save :: ncidfd
  INTEGER,save :: variddt,varididvt,varididvp
  REAL :: dtv(1,1,1),adv2(1,1,1),adv(1,1,1)
  REAL :: zdtvr,ziadvtrac,ziadvtrac2
  INTEGER :: status,irec
  REAL :: rcode
  INTEGER :: epais(2), debut(2)
  
  IF (irec == 0) THEN
      
      ncidfd=NCOPN('defstoke.nc',NCNOWRIT,rcode)
      
      variddt=NCVID(ncidfd,'dtvr',rcode)
      PRINT*,'ncidfd,variddt',ncidfd,variddt
      
      varididvt=NCVID(ncidfd,'istdyn',rcode)
      PRINT*,'ncidfd,varididvt',ncidfd,varididvt
      
      varididvp=NCVID(ncidfd,'istphy',rcode)
      PRINT*,'ncidfd,varididvp',ncidfd,varididvp
      
      ! lecture de zdtvr et ziadvtrac
      
      epais(1) = 1
      epais(2) = 1
      debut(1) = 1
      debut(2) = 1
      
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidfd,variddt,debut,epais,dtv)
#else
      status=NF_GET_VARA_REAL(ncidfd,variddt,debut,epais,dtv)
#endif
      zdtvr=dtv(1,1,1)
      
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidfd,varididvt,debut,epais,adv)
#else
      status=NF_GET_VARA_REAL(ncidfd,varididvt,debut,epais,adv)
#endif
      ziadvtrac= adv(1,1,1)
      
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidfd,varididvp,debut,epais,adv2)
#else
      status=NF_GET_VARA_REAL(ncidfd,varididvp,debut,epais,adv2)
#endif
      ziadvtrac2= adv2(1,1,1)
      
      WRITE(*,*) 'ds read_dstoke zdtvr = ', zdtvr
      WRITE(*,*) 'ds read_dstoke ziadvtrac = ', ziadvtrac
      WRITE(*,*) 'ds read_dstoke ziadvtrac2 = ', ziadvtrac2
      
      !           status=NF_CLOSE(ncidfd)
      
  ELSE
      STOP'Pas bon irec ne 0'
  ENDIF
  
  RETURN
  
END SUBROUTINE read_dstoke

