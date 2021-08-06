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

SUBROUTINE coefpoly ( Xf1, Xf2, Xprim1, Xprim2, xtild1,xtild2, a0,a1,a2,a3 )
  IMPLICIT NONE
  !
  !   ...  Auteur :   P. Le Van  ...
  !
  !
  !    Calcul des coefficients a0, a1, a2, a3 du polynome de degre 3 qui
  !      satisfait aux 4 equations  suivantes :
  !
  !    a0 + a1*xtild1 + a2*xtild1*xtild1 + a3*xtild1*xtild1*xtild1 = Xf1
  !    a0 + a1*xtild2 + a2*xtild2*xtild2 + a3*xtild2*xtild2*xtild2 = Xf2
  !               a1  +     2.*a2*xtild1 +     3.*a3*xtild1*xtild1 = Xprim1
  !               a1  +     2.*a2*xtild2 +     3.*a3*xtild2*xtild2 = Xprim2
  !
  !  On en revient a resoudre un systeme de 4 equat.a 4 inconnues a0,a1,a2,a3

  REAL*8 :: Xf1, Xf2,Xprim1,Xprim2, xtild1,xtild2
  REAL*8 :: a1,a2,a3,a0, xtil1car, xtil2car,derr,x1x2car

  xtil1car = xtild1 * xtild1
  xtil2car = xtild2 * xtild2 
  
  derr= 2. *(Xf2-Xf1)/( xtild1-xtild2)
  
  x1x2car = ( xtild1-xtild2)*(xtild1-xtild2)
  
  a3 = (derr + Xprim1+Xprim2 )/x1x2car
  a2     = ( Xprim1 - Xprim2 + 3.* a3 * ( xtil2car-xtil1car ) ) / &
     (  2.* ( xtild1 - xtild2 )  )

  a1     = Xprim1 -3.* a3 * xtil1car     -2.* a2 * xtild1
  a0     =  Xf1 - a3 * xtild1* xtil1car -a2 * xtil1car - a1 *xtild1
  
  RETURN
END SUBROUTINE coefpoly
