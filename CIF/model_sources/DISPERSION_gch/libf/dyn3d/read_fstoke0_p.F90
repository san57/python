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

SUBROUTINE read_fstoke0_p(irec,zrec,zim,zjm,zlm, aready,phis, masse,pbaru,pbarv,w,teta,phi)
  USE parallel
  USE mod_const_mpi
  USE vampir
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INCLUDE 'dimensions.h'
  INCLUDE 'paramet.h'
  
  INTEGER :: irec,zrec,zim,zjm,zlm
  REAL :: pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm)
  REAL :: teta(iip1,jjp1,llm),phis(iip1,jjp1),phi(iip1,jjp1,llm)
  REAL :: masse(iip1,jjp1,llm),w(iip1,jjp1,llm)
  REAL :: aready(iip1,jjp1)

  CALL VTb(VTread_fstoke0)
  CALL read_fstoke0(irec, zrec,zim,zjm,zlm, aready,phis,masse,pbaru,pbarv,w,teta,phi)
  CALL VTe(VTread_fstoke0)
  
END SUBROUTINE read_fstoke0_p
        
