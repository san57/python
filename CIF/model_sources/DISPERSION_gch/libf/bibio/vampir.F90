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

MODULE Vampir
  
  INTEGER,PARAMETER :: VTcaldyn=1
  INTEGER,PARAMETER :: VTintegre=2
  INTEGER,PARAMETER :: VTadvection=3
  INTEGER,PARAMETER :: VTdissipation=4
  INTEGER,PARAMETER :: VThallo=5
  INTEGER,PARAMETER :: VTphysiq=6
  INTEGER,PARAMETER :: VTinca=7
  INTEGER,PARAMETER :: VTvlx=8
  INTEGER,PARAMETER :: VTvly=9
  INTEGER,PARAMETER :: VTvlz=10
  INTEGER,PARAMETER :: VTlecfluxnc=11
  INTEGER,PARAMETER :: VTread_source=12
  INTEGER,PARAMETER :: VTread_fstoke0=13
  INTEGER,PARAMETER :: VTread_pstoke0=14
  INTEGER,PARAMETER :: VTadvection_ad=15
CONTAINS
  
  SUBROUTINE InitVampir
    IMPLICIT NONE
#ifdef USE_VT
    !    include 'VT.inc'
    INTEGER :: ierr
    
    CALL VTSYMDEF(VTcaldyn,"caldyn","caldyn",ierr)
    CALL VTSYMDEF(VTintegre,"integre","integre",ierr)
    CALL VTSYMDEF(VTadvection,"advection","advection",ierr)
    CALL VTSYMDEF(VTdissipation,"dissipation","dissipation",ierr)
    CALL VTSYMDEF(VThallo,"hallo","hallo",ierr)
    CALL VTSYMDEF(VTphysiq,"physiq","physiq",ierr)
    CALL VTSYMDEF(VTinca,"inca","inca",ierr)
    CALL VTSYMDEF(VTvlx,"vlx","vlx",ierr)
    CALL VTSYMDEF(VTvly,"vly","vly",ierr)
    CALL VTSYMDEF(VTvlz,"vlz","vlz",ierr)
    CALL VTSYMDEF(VTlecfluxnc,"lecfluxnc","lecfluxnc",ierr)
    call VTSYMDEF(VTread_source,"read_source","read_source",ierr)
    CALL VTSYMDEF(VTread_fstoke0,"read_fstoke0","read_fstoke0",ierr)
    CALL VTSYMDEF(VTread_pstoke0,"read_pstoke0","read_pstoke0",ierr)
    CALL VTSYMDEF(VTadvection_ad,"VTadvection_ad","VTadvection_ad",ierr)
#endif  
  END SUBROUTINE InitVampir
  
  SUBROUTINE VTb(number) 
    IMPLICIT NONE
    INTEGER :: number
#ifdef USE_VT    
    !    include 'VT.inc'
    INTEGER :: ierr
    
    CALL VTBEGIN(number,ierr)
#endif
  END SUBROUTINE VTb
  
  SUBROUTINE VTe(number)
    IMPLICIT NONE
    INTEGER :: Number
#ifdef USE_VT    
    !    include 'VT.inc'
    INTEGER :: ierr
    
    CALL VTEND(number,ierr)
#endif    
    
  END SUBROUTINE VTe
  
END MODULE Vampir
  
