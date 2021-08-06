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

SUBROUTINE fxysinus (rlatu,yprimu,rlatv,yprimv,rlatu1,yprimu1, &
   rlatu2,yprimu2, &
   rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)


  IMPLICIT NONE
  !
  !     Calcul  des longitudes et des latitudes  pour une fonction f(x,y)
  !            avec y = Asin( j )  .
  !
  !     Auteur  :  P. Le Van
  !
  !
 include "dimensions.h"
 include "paramet.h"
 include "comconst.h"
  
  INTEGER :: i,j
  
  REAL :: rlatu(jjp1), yprimu(jjp1),rlatv(jjm), yprimv(jjm), &
     rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
  REAL :: rlonu(iip1),xprimu(iip1),rlonv(iip1),xprimv(iip1), &
     rlonm025(iip1),xprimm025(iip1), rlonp025(iip1),xprimp025(iip1)
  
 include "fxy_sin.h"
  
  
  !    ......  calcul  des  latitudes  et de y'   .....
  !
  DO j = 1, jjm + 1 
    rlatu(j) = fy    ( FLOAT( j )        )
    yprimu(j) = fyprim( FLOAT( j )        )
  ENDDO
  
  
  DO j = 1, jjm
    
    rlatv(j)  = fy    ( FLOAT( j ) + 0.5  )
    rlatu1(j) = fy    ( FLOAT( j ) + 0.25 ) 
    rlatu2(j) = fy    ( FLOAT( j ) + 0.75 ) 
    
    yprimv(j)  = fyprim( FLOAT( j ) + 0.5  ) 
    yprimu1(j) = fyprim( FLOAT( j ) + 0.25 )
    yprimu2(j) = fyprim( FLOAT( j ) + 0.75 )
    
  ENDDO

  !
  !     .....  calcul   des  longitudes et de  x'   .....
  !
  DO i = 1, iim + 1
    rlonv(i)     = fx    (   FLOAT( i )          )
    rlonu(i)     = fx    (   FLOAT( i ) + 0.5    )
    rlonm025(i)     = fx    (   FLOAT( i ) - 0.25  )
    rlonp025(i)     = fx    (   FLOAT( i ) + 0.25  )
    
    xprimv  (i)    = fxprim (  FLOAT( i )          )
    xprimu  (i)    = fxprim (  FLOAT( i ) + 0.5    )
    xprimm025(i)    = fxprim (  FLOAT( i ) - 0.25   )
    xprimp025(i)    = fxprim (  FLOAT( i ) + 0.25   )
  ENDDO

  !
  RETURN
END SUBROUTINE fxysinus

