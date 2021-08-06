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

!     ------------------------------------------------------------------
!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!     ------------------------------------------------------------------
      REAL PTARG, PDELARG, P5ARG, PQSARG, PCOARG
      REAL FOEEW, FOEDE, qsats, qsatl, dqsats, dqsatl
      LOGICAL thermcep
      PARAMETER (thermcep=.TRUE.)
!
      FOEEW ( PTARG,PDELARG ) = EXP ( &
      (R3LES*(1.-PDELARG)+R3IES*PDELARG) * (PTARG-RTT) &
      /(PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG)) )
!
      FOEDE ( PTARG,PDELARG,P5ARG,PQSARG,PCOARG ) = PQSARG*PCOARG*P5ARG &
      / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG))**2
!
      qsats(ptarg) = 100.0 * 0.622 * 10.0 &
      ** (2.07023 - 0.00320991 * ptarg &
      - 2484.896 / ptarg + 3.56654 * LOG10(ptarg))
      qsatl(ptarg) = 100.0 * 0.622 * 10.0 &
      ** (23.8319 - 2948.964 / ptarg &
      - 5.028 * LOG10(ptarg) &
      - 29810.16 * EXP( - 0.0699382 * ptarg) &
      + 25.21935 * EXP( - 2999.924 / ptarg))
!
      dqsats(ptarg,pqsarg) = RLVTT/RCPD*pqsarg * (3.56654/ptarg &
      +2484.896*LOG(10.)/ptarg**2 &
      -0.00320991*LOG(10.))
      dqsatl(ptarg,pqsarg) = RLVTT/RCPD*pqsarg*LOG(10.)* &
      (2948.964/ptarg**2-5.028/LOG(10.)/ptarg &
      +25.21935*2999.924/ptarg**2*EXP(-2999.924/ptarg) &
      +29810.16*0.0699382*EXP(-0.0699382*ptarg))
