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

SUBROUTINE pression( ngrid, ap, bp, ps, p )
  !

  !      Auteurs : P. Le Van , Fr.Hourdin  .

  !  ************************************************************************
  !     Calcule la pression p(l) aux differents niveaux l = 1 ( niveau du
  !     sol) a l = llm +1 ,ces niveaux correspondant aux interfaces des (llm) 
  !     couches , avec  p(ij,llm +1) = 0.  et p(ij,1) = ps(ij)  .      
  !  ************************************************************************
  !
  USE parallel
  IMPLICIT NONE
  !
 include "dimensions.h"
 include "paramet.h"
  !
  INTEGER :: ngrid
  INTEGER :: l,ij
  
  REAL :: ap( llmp1 ), bp( llmp1 ), ps( ngrid ), p( ngrid,llmp1 ) 
  
  DO    l    = 1, llmp1
    DO  ij   = ij_begin, ij_end
      p(ij,l) = ap(l) + bp(l) * ps(ij)
    ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE pression
