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

SUBROUTINE INTEX_TL( &
   klevi , &  ! in
   klevf , &  ! in
   presi , &  ! in
   presf , &  ! in
   veci  , & ! in
   vecf  , & ! out
   veci_tl  , &  ! in
   vecf_tl  )   ! out
    
  IMPLICIT NONE
  
  INTEGER, INTENT(in) :: klevi      ! number of levels of the initial grid
  INTEGER, INTENT(in) :: klevf      ! number of levels of the final grid
  REAL(kind=8), INTENT(in), DIMENSION(klevi)  :: presi ! initial grid
  REAL(kind=8), INTENT(in), DIMENSION(klevf)  :: presf ! final grid
  REAL(kind=8), INTENT(in), DIMENSION(klevi)  :: veci  ! initial vec array
  REAL(kind=8), INTENT(out),DIMENSION(klevf)  :: vecf  ! final vec array
  REAL(kind=8), INTENT(in), DIMENSION(klevi)  :: veci_tl  ! initial vec array
  REAL(kind=8), INTENT(out), DIMENSION(klevf) :: vecf_tl  ! final vec array
  INTEGER :: jki, jkf
  REAL (kind=8)   :: slope, slope_tl, t1, t1_tl, t2, t2_tl, p1, p2, lp1, lp2
  REAL(kind=8), DIMENSION(klevi)  :: lpresi
  REAL(kind=8), DIMENSION(klevf)  :: lpresf
  !
  !- End of header --------------------------------------------------------
  
  vecf(:) = -1000.
  lpresi(:) = LOG( presi(:) )
  lpresf(:) = LOG( presf(:) )
  
  DO jkf = 1,klevf
    DO jki = 1,klevi-1
      p1 = presi(jki)
      p2 = presi(jki+1)
      lp1 = lpresi(jki)
      lp2 = lpresi(jki+1)
      IF (presf(jkf) >= p1 .AND. presf(jkf) < p2) THEN
          t1 = veci(jki)
          t2 = veci(jki+1)
	  t1_tl = veci_tl(jki)
          t2_tl = veci_tl(jki+1)
	  slope_tl = (t1_tl-t2_tl)/(lp1-lp2)
	  slope = (t1-t2)/(lp1-lp2)
          
          vecf_tl(jkf) = t1_tl + slope_tl*(lpresf(jkf)-lp1)
	  vecf(jkf) = t1 + slope*(lpresf(jkf)-lp1)
          !
      ELSE IF (jki == 1 .AND. presf(jkf) < p1) THEN
          vecf_tl(jkf) = veci_tl(jki)
	  vecf(jkf) = veci(jki)
      ELSE IF (jki == (klevi-1) .AND. vecf(jkf) == -1000. ) THEN
          vecf_tl(jkf) = veci_tl(klevi)
	  vecf(jkf) = veci(klevi)
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE INTEX_TL
