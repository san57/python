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

Subroutine INTEX_AD( &
 klevi , & ! in
 klevf , & ! in
 presi , & ! in
 presf , & ! in
 veci  , & ! in
 vecf  , & ! out
 veci_ad,& ! inout
 vecf_ad)  ! inout

Implicit None

Integer, Intent(in) :: klevi      ! number of levels of the initial grid
Integer, Intent(in) :: klevf      ! number of levels of the final grid
Real, Intent(in), Dimension(klevi)  :: presi ! initial grid
Real, Intent(in), Dimension(klevf)  :: presf ! final grid
Real, Intent(in), Dimension(klevi)  :: veci  ! initial vec array
Real, Intent(out), Dimension(klevf) :: vecf  ! final vec array
Real, Intent(inout), Dimension(klevi)  :: veci_ad  ! initial vec array
Real, Intent(inout), Dimension(klevf) :: vecf_ad  ! final vec array
Integer :: jki, jkf
Real    :: slope, slope_ad, t1, t1_ad, t2, t2_ad, p1, p2, lp1, lp2
Real, Dimension(klevi)  :: lpresi
Real, Dimension(klevf)  :: lpresf
!
!- End of header --------------------------------------------------------

vecf(:) = -1000.
lpresi(:) = Log( presi(:) )
lpresf(:) = Log( presf(:) )

Do jkf = 1,klevf
   Do jki = 1,klevi-1
      p1 = presi(jki)
      p2 = presi(jki+1)
      lp1 = lpresi(jki)
      lp2 = lpresi(jki+1)
      If (presf(jkf) >= p1 .And. presf(jkf) < p2) Then
         t1 = veci(jki)
         t2 = veci(jki+1)
         slope = (t1-t2)/(lp1-lp2)
         vecf(jkf) = t1 + slope*(lpresf(jkf)-lp1)
         !
      Else If (jki == 1 .And. presf(jkf) < p1) Then
         vecf(jkf) = veci(jki)
      End If
   End Do
End Do
Do jkf = klevf,1,-1
   Do jki = klevi-1,1,-1
      p1 = presi(jki)
      p2 = presi(jki+1)
      lp1 = lpresi(jki)
      lp2 = lpresi(jki+1)
      If (presf(jkf) >= p1 .And. presf(jkf) < p2) Then
        t1_ad = 0.
        t2_ad = 0.
        slope_ad = 0.
        t1_ad = t1_ad + vecf_ad(jkf)
        slope_ad = slope_ad + (lpresf(jkf)-lp1) * vecf_ad(jkf)
        vecf_ad(jkf) = 0.
        t1_ad = t1_ad + slope_ad / (lp1-lp2)
        t2_ad = t2_ad - slope_ad / (lp1-lp2)
        veci_ad(jki+1) = veci_ad(jki+1) + t2_ad
        t2_ad = 0.
        veci_ad(jki) = veci_ad(jki) + t1_ad
        t1_ad = 0.
      Else If (jki == 1 .And. presf(jkf) < p1) Then
        veci_ad(jki) = veci_ad(jki) + vecf_ad(jkf)
        vecf_ad(jkf) = 0.
      Else If (jki == (klevi-1) .And. vecf(jkf) == -1000. ) Then
        veci_ad(klevi) = veci_ad(klevi) + vecf_ad(jkf)
        vecf_ad(jkf) = 0.
        vecf(jkf) = veci(klevi)
      End If
   End Do
End Do
Return
End Subroutine INTEX_AD
