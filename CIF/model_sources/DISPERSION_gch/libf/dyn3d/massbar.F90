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

SUBROUTINE massbar(  masse, massebx, masseby )
  !
  ! **********************************************************************
  !
  !  Calcule les moyennes en x et  y de la masse d'air dans chaque maille.
  ! **********************************************************************
  !    Auteurs : P. Le Van , Fr. Hourdin  .
  !   ..........
  !
  !  ..  masse                 est  un argum. d'entree  pour le s-pg ...
  !  ..  massebx,masseby      sont des argum. de sortie pour le s-pg ...
  !     
  !
  USE parallel
  IMPLICIT NONE
  !
 include "dimensions.h"
 include "paramet.h"
 include "comconst.h"
 include "comgeom.h"
  !
  REAL :: masse( ip1jmp1,llm ), massebx( ip1jmp1,llm )  , &
     masseby(   ip1jm,llm )
  INTEGER :: ij,l,ijb,ije
  !
  !
  !   Methode pour calculer massebx et masseby .
  !   ----------------------------------------
  !
  !    A chaque point scalaire P (i,j) est affecte 4 coefficients d'aires
  !       alpha1(i,j)  calcule  au point ( i+1/4,j-1/4 )
  !       alpha2(i,j)  calcule  au point ( i+1/4,j+1/4 )
  !       alpha3(i,j)  calcule  au point ( i-1/4,j+1/4 )
  !       alpha4(i,j)  calcule  au point ( i-1/4,j-1/4 )
  !
  !    Avec  alpha1(i,j) = area(i+1/4,j-1/4)/ area(i,j)
  !
  !    N.B .  Pour plus de details, voir s-pg  ...  iniconst ...
  !
  !
  !
  !   alpha4 .         . alpha1    . alpha4
  !    (i,j)             (i,j)       (i+1,j)
  !
  !             P .        U .          . P
  !           (i,j)       (i,j)         (i+1,j)
  !
  !   alpha3 .         . alpha2    .alpha3 
  !    (i,j)              (i,j)     (i+1,j)
  !
  !             V .        Z .          . V
  !           (i,j)
  !
  !   alpha4 .         . alpha1    .alpha4
  !   (i,j+1)            (i,j+1)   (i+1,j+1) 
  !
  !             P .        U .          . P
  !          (i,j+1)                    (i+1,j+1)
  !
  !
  !
  !                       On  a :
  !
  !    massebx(i,j) = masse(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))   +
  !                   masse(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
  !     localise  au point  ... U (i,j) ...
  !
  !    masseby(i,j) = masse(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )  +
  !                   masse(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
  !     localise  au point  ... V (i,j) ...
  !
  !
  !=======================================================================
      
      
      
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)  
  DO l = 1 , llm
    !
    ijb=ij_begin
    ije=ij_end
    
    DO  ij = ijb, ije - 1
      massebx(ij,l) =  masse( ij, l) * alpha1p2( ij  )     +  &
         masse(ij+1, l) * alpha3p4(ij+1 )
    ENDDO
    
    !    .... correction pour massebx( iip1,j) .....
    !    ...    massebx(iip1,j)= massebx(1,j) ...
    !
    !DIR$ IVDEP
    
    
    
    DO  ij = ijb+iim, ije+iim, iip1
      massebx( ij,l ) = massebx( ij - iim,l )
    ENDDO
    
    
    
    ijb=ij_begin
    ije=ij_end
    IF (south_pole) ije=ij_end-iip1
    
    DO  ij = ijb,ije
      masseby( ij,l ) = masse(  ij   , l ) * alpha2p3(   ij    )  + &
         masse(ij+iip1, l ) * alpha1p4( ij+iip1 )
    ENDDO
    
  ENDDO
!$OMP END DO NOWAIT
!
  RETURN
END SUBROUTINE massbar
