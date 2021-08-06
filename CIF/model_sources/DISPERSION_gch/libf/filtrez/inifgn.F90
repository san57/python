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

SUBROUTINE inifgn(dv)
  !  
  !    ...  H.Upadyaya , O.Sharma  ... 
  !
  IMPLICIT NONE
  !
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "comgeom.h"
  INCLUDE "serre.h"
  
  !
  REAL :: vec(iim,iim),vec1(iim,iim)
  REAL :: dlonu(iim),dlonv(iim)
  REAL :: du(iim),dv(iim),d(iim)
  REAL :: pi
  INTEGER :: i,j,k,imm1,nrot
  !
  INCLUDE "coefils.h"
  !
  EXTERNAL SSUM, acc,eigen,jacobi
  REAL :: SSUM
  !
  
  imm1  = iim -1
  pi = 2.* ASIN(1.)
  !
  DO i=1,iim
    dlonu(i)=  xprimu( i )
    dlonv(i)=  xprimv( i )
  ENDDO
  
  DO i=1,iim
    sddv(i)   = SQRT(dlonv(i))
    sddu(i)   = SQRT(dlonu(i))
    unsddu(i) = 1./sddu(i)
    unsddv(i) = 1./sddv(i)
  ENDDO
  !
  DO j=1,iim
    DO i=1,iim
      vec(i,j)     = 0.
      vec1(i,j)    = 0.
      eignfnv(i,j) = 0.
      eignfnu(i,j) = 0.
    ENDDO
  ENDDO
  !
  !
  eignfnv(1,1)    = -1.
  eignfnv(iim,1)  =  1.
  DO i=1,imm1
    eignfnv(i+1,i+1)= -1.
    eignfnv(i,i+1)  =  1.
  ENDDO
  DO j=1,iim
    DO i=1,iim
      eignfnv(i,j) = eignfnv(i,j)/(sddu(i)*sddv(j))
    ENDDO
  ENDDO
  DO j=1,iim
    DO i=1,iim
      eignfnu(i,j) = -eignfnv(j,i)
    ENDDO
  ENDDO
  !
#ifdef CRAY
  CALL MXM(eignfnu,iim,eignfnv,iim,vec ,iim)
  CALL MXM(eignfnv,iim,eignfnu,iim,vec1,iim)
#else
  DO j = 1, iim
    DO i = 1, iim
      vec (i,j) = 0.0
      vec1(i,j) = 0.0
      DO k = 1, iim
        vec (i,j) = vec(i,j)  + eignfnu(i,k) * eignfnv(k,j)
        vec1(i,j) = vec1(i,j) + eignfnv(i,k) * eignfnu(k,j)
      ENDDO
    ENDDO
  ENDDO
#endif
  
  !
  CALL jacobi(vec,iim,iim,dv,eignfnv,nrot)
  CALL acc(eignfnv,d,iim)
  CALL eigen_sort(dv,eignfnv,iim,iim)
  !
  CALL jacobi(vec1,iim,iim,du,eignfnu,nrot)
  CALL acc(eignfnu,d,iim)
  CALL eigen_sort(du,eignfnu,iim,iim)
  
  !c   ancienne version avec appels IMSL
  !
  !     CALL MXM(eignfnu,iim,eignfnv,iim,vec,iim)
  !     CALL MXM(eignfnv,iim,eignfnu,iim,vec1,iim)
  !     CALL EVCSF(iim,vec,iim,dv,eignfnv,iim)
  !     CALL acc(eignfnv,d,iim)
  !     CALL eigen(eignfnv,dv)
  !
  !     CALL EVCSF(iim,vec1,iim,du,eignfnu,iim)
  !     CALL acc(eignfnu,d,iim)
  !     CALL eigen(eignfnu,du)
  
  RETURN
END SUBROUTINE inifgn

