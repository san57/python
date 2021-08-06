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

MODULE mod_phys_lmdz_para
  USE mod_phys_lmdz_transfert_para
  USE mod_phys_lmdz_mpi_data
  USE mod_phys_lmdz_omp_data
    
  INTEGER,SAVE :: klon_loc
  LOGICAL,SAVE :: is_sequential
  LOGICAL,SAVE :: is_parallel
  
!$OMP THREADPRIVATE(klon_loc)
  
CONTAINS

  SUBROUTINE Init_phys_lmdz_para(iim,jjp1,nb_proc,distrib)
  IMPLICIT NONE
    INTEGER,INTENT(in) :: iim
    INTEGER,INTENT(in) :: jjp1
    INTEGER,INTENT(in) :: nb_proc
    INTEGER,INTENT(in) :: distrib(0:nb_proc-1)

    CALL Init_phys_lmdz_mpi_data(iim,jjp1,nb_proc,distrib)
!$OMP PARALLEL
    CALL Init_phys_lmdz_omp_data(klon_mpi)
    klon_loc=klon_omp
    CALL Test_transfert
!$OMP END PARALLEL    
     IF (is_ok_mpi .OR. is_ok_omp) THEN
       is_sequential=.FALSE.
       is_parallel=.TRUE.
     ELSE
       is_sequential=.TRUE.
       is_parallel=.FALSE.
     ENDIF
     
  END SUBROUTINE Init_phys_lmdz_para

  SUBROUTINE Test_transfert
  USE mod_grid_phy_lmdz
  IMPLICIT NONE
  
    REAL :: Test_Field1d_glo(klon_glo,nbp_lev)
    REAL :: tmp1d_glo(klon_glo,nbp_lev)
    REAL :: Test_Field2d_glo(nbp_lon,nbp_lat,nbp_lev)
    REAL :: Test_Field1d_loc(klon_loc,nbp_lev)
    REAL :: Test_Field2d_loc(nbp_lon,jj_nb,nbp_lev)
    REAL :: CheckSum
    
    INTEGER :: i,l
  
    Test_Field1d_glo = 0.
    Test_Field2d_glo = 0.
    Test_Field1d_loc = 0.
    Test_Field2d_loc = 0.
  
    IF (is_mpi_root) THEN
!$OMP MASTER
      DO l=1,nbp_lev
        DO i=1,klon_glo
!          Test_Field1d_glo(i,l)=MOD(i,10)+10*(l-1)
           Test_Field1d_glo(i,l)=1
        ENDDO
      ENDDO
!$OMP END MASTER  
    ENDIF
  
    CALL Scatter(Test_Field1d_glo,Test_Field1d_loc)
    CALL Gather(Test_Field1d_loc,tmp1d_glo)
  
    IF (is_mpi_root) THEN
!$OMP MASTER  
      Checksum=sum(Test_Field1d_glo-tmp1d_glo)
      PRINT *, "------> Checksum =",Checksum," MUST BE 0"
!$OMP END MASTER
    ENDIF
    
    CALL grid1dTo2d_glo(Test_Field1d_glo,Test_Field2d_glo)
    CALL scatter2D(Test_Field2d_glo,Test_Field1d_loc)
    CALL gather2d(Test_Field1d_loc,Test_Field2d_glo)
    CALL grid2dTo1d_glo(Test_Field2d_glo,tmp1d_glo)

    IF (is_mpi_root) THEN
!$OMP MASTER  
      Checksum=sum(Test_Field1d_glo-tmp1d_glo)
      PRINT *, "------> Checksum =",Checksum," MUST BE 0"
!$OMP END MASTER
    ENDIF

    CALL bcast(Test_Field1d_glo)
    CALL reduce_sum(Test_Field1d_glo,tmp1d_glo)

    IF (is_mpi_root) THEN
!$OMP MASTER  
      Checksum=sum(Test_Field1d_glo*omp_size*mpi_size-tmp1d_glo)
      PRINT *, "------> Checksum =",Checksum," MUST BE 0"
!$OMP END MASTER
    ENDIF
    
     
   END SUBROUTINE Test_transfert
  
END MODULE mod_phys_lmdz_para
    
