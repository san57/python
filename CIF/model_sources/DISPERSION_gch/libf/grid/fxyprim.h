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

!-----------------------------------------------------------------------
! INCLUDE 'fxyprim.h'
!
!    ................................................................
!    ................  Fonctions in line  ...........................
!    ................................................................
!
      REAL :: fy, fx, fxprim, fyprim
      REAL :: ri, rj
!
!
      fy    ( rj ) =    pi/FLOAT(jjm) * ( 0.5 * FLOAT(jjm) +  1. - rj  )
      fyprim( rj ) =    pi/FLOAT(jjm)

!     fy(rj)=ASIN(1.+2.*((1.-rj)/FLOAT(jjm)))
!     fyprim(rj)=1./SQRT((rj-1.)*(jjm+1.-rj))

      fx    ( ri ) = 2.*pi/FLOAT(iim) * ( ri - 0.5*  FLOAT(iim) - 1. )
!     fx    ( ri ) = 2.*pi/FLOAT(iim) * ( ri - 0.5* ( FLOAT(iim) + 1.) )
      fxprim( ri ) = 2.*pi/FLOAT(iim)
!
!
!    La valeur de pi est passee par le common/const/ou /const2/ .
!    Sinon, il faut la calculer avant d'appeler ces fonctions .
!
!   ----------------------------------------------------------------
!     Fonctions a changer eventuellement, selon x(x) et y(y) choisis .
!   -----------------------------------------------------------------
!
!    .....  ici, on a l'application particuliere suivante   ........
!
!                **************************************
!                **     x = 2. * pi/iim *  X         **
!                **     y =      pi/jjm *  Y         **
!                **************************************
!
!   ..................................................................
!   ..................................................................
!
!
!
!-----------------------------------------------------------------------
