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

MODULE mod_phys_lmdz_transfert_para

  USE mod_phys_lmdz_mpi_transfert
  USE mod_phys_lmdz_omp_transfert 



  INTERFACE bcast
    MODULE PROCEDURE bcast_c,                                     &
                     bcast_i,bcast_i1,bcast_i2,bcast_i3,bcast_i4, &
                     bcast_r,bcast_r1,bcast_r2,bcast_r3,bcast_r4, &
		     bcast_l,bcast_l1,bcast_l2,bcast_l3,bcast_l4
  END INTERFACE

  INTERFACE scatter
    MODULE PROCEDURE scatter_i,scatter_i1,scatter_i2,scatter_i3, &
                     scatter_r,scatter_r1,scatter_r2,scatter_r3, &
		     scatter_l,scatter_l1,scatter_l2,scatter_l3
  END INTERFACE

  
  INTERFACE gather
    MODULE PROCEDURE gather_i,gather_i1,gather_i2,gather_i3, &
                     gather_r,gather_r1,gather_r2,gather_r3, &
		     gather_l,gather_l1,gather_l2,gather_l3  
  END INTERFACE
  
  INTERFACE scatter2D
    MODULE PROCEDURE scatter2D_i,scatter2D_i1,scatter2D_i2,scatter2D_i3, &
                     scatter2D_r,scatter2D_r1,scatter2D_r2,scatter2D_r3, &
		     scatter2D_l,scatter2D_l1,scatter2D_l2,scatter2D_l3
  END INTERFACE

  INTERFACE gather2D
    MODULE PROCEDURE gather2D_i,gather2D_i1,gather2D_i2,gather2D_i3, &
                     gather2D_r,gather2D_r1,gather2D_r2,gather2D_r3, &
		     gather2D_l,gather2D_l1,gather2D_l2,gather2D_l3
  END INTERFACE 
  
  INTERFACE reduce_sum
    MODULE PROCEDURE reduce_sum_i,reduce_sum_i1,reduce_sum_i2,reduce_sum_i3,reduce_sum_i4, &
                     reduce_sum_r,reduce_sum_r1,reduce_sum_r2,reduce_sum_r3,reduce_sum_r4
  END INTERFACE 

   
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Broadcast --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Strings -- !!

  SUBROUTINE bcast_c(var)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_c

!! -- Les entiers -- !!
  
  SUBROUTINE bcast_i(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i

  SUBROUTINE bcast_i1(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i1


  SUBROUTINE bcast_i2(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i2


  SUBROUTINE bcast_i3(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i3


  SUBROUTINE bcast_i4(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_i4

 
!! -- Les reels -- !!
  
  SUBROUTINE bcast_r(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var

!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r

  SUBROUTINE bcast_r1(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r1


  SUBROUTINE bcast_r2(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r2


  SUBROUTINE bcast_r3(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r3


  SUBROUTINE bcast_r4(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_r4 


!! -- Les booleens -- !!
  
  SUBROUTINE bcast_l(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l

  SUBROUTINE bcast_l1(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l1


  SUBROUTINE bcast_l2(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l2


  SUBROUTINE bcast_l3(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l3


  SUBROUTINE bcast_l4(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:,:)
   
!$OMP MASTER
    CALL bcast_mpi(Var)
!$OMP END MASTER
    CALL bcast_omp(Var)
    
  END SUBROUTINE bcast_l4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter_i(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    INTEGER,DIMENSION(klon_mpi) :: Var_tmp
    
!$OMP MASTER
      CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER

      CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_i


  SUBROUTINE scatter_i1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(Varout,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      INTEGER,DIMENSION(klon_mpi,s1) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body

  END SUBROUTINE scatter_i1


  SUBROUTINE scatter_i2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(Varout,2),SIZE(Varout,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      INTEGER,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body
    
  END SUBROUTINE scatter_i2


  SUBROUTINE scatter_i3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      INTEGER,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body
    
  END SUBROUTINE scatter_i3


  SUBROUTINE scatter_r(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    REAL,DIMENSION(klon_mpi) :: Var_tmp
    
!$OMP MASTER
      CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER

      CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_r


  SUBROUTINE scatter_r1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(Varout,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      REAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body
    
  END SUBROUTINE scatter_r1


  SUBROUTINE scatter_r2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(Varout,2),SIZE(Varout,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      REAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body
    
  END SUBROUTINE scatter_r2


  SUBROUTINE scatter_r3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      REAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body
    
  END SUBROUTINE scatter_r3
  
  

  SUBROUTINE scatter_l(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    LOGICAL,DIMENSION(klon_mpi) :: Var_tmp
    
!$OMP MASTER
      CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER

      CALL scatter_omp(Var_tmp,Varout)
    
  END SUBROUTINE scatter_l


  SUBROUTINE scatter_l1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(Varout,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      LOGICAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body
    
  END SUBROUTINE scatter_l1


  SUBROUTINE scatter_l2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(Varout,2),SIZE(Varout,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      LOGICAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body

  END SUBROUTINE scatter_l2


  SUBROUTINE scatter_l3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(Varout,2),SIZE(Varout,3),SIZE(Varout,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      LOGICAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
!$OMP MASTER
        CALL scatter_mpi(VarIn,Var_tmp)
!$OMP END MASTER
        CALL scatter_omp(Var_tmp,Varout)
    END SUBROUTINE body
    
  END SUBROUTINE scatter_l3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Gather   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!! --> Les entiers

  SUBROUTINE gather_i(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    INTEGER, DIMENSION(klon_mpi) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
    CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER
  
  END SUBROUTINE gather_i


  SUBROUTINE gather_i1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      INTEGER,DIMENSION(klon_mpi,s1) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_i1


  SUBROUTINE gather_i2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      INTEGER,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_i2


  SUBROUTINE gather_i3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      INTEGER,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_i3


!!!!! --> Les reels

  SUBROUTINE gather_r(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    REAL, DIMENSION(klon_mpi) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
    CALL gather_mpi(Var_tmp,VarOut)
!$OMP END MASTER
  
  END SUBROUTINE gather_r


  SUBROUTINE gather_r1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      REAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_r1


  SUBROUTINE gather_r2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      REAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_r2


  SUBROUTINE gather_r3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      REAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_r3


!!!!! --> Les booleens

  SUBROUTINE gather_l(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    LOGICAL, DIMENSION(klon_mpi) :: Var_tmp
    
    CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
    CALL gather_mpi(Var_tmp,VarOut)
!$OMP END MASTER
  
  END SUBROUTINE gather_l


  SUBROUTINE gather_l1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      LOGICAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_l1


  SUBROUTINE gather_l2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      LOGICAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_l2


  SUBROUTINE gather_l3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      LOGICAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather_mpi(Var_tmp,Varout)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE gather_l3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter2D   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! --> Les entiers

  SUBROUTINE scatter2D_i(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    INTEGER,DIMENSION(klon_mpi) :: Var_tmp    

!$OMP MASTER    
    CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
    CALL scatter_omp(Var_tmp,VarOut)

  END SUBROUTINE scatter2D_i


  SUBROUTINE scatter2D_i1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      INTEGER,DIMENSION(klon_mpi,s1) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_i1
  

  SUBROUTINE scatter2D_i2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2),SIZE(VarOut,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      INTEGER,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_i2  


  SUBROUTINE scatter2D_i3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2),SIZE(VarOut,3),SIZE(VarOut,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      INTEGER,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_i3
  

!!!!! --> Les reels

  SUBROUTINE scatter2D_r(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    REAL,DIMENSION(klon_mpi) :: Var_tmp    

!$OMP MASTER    
    CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
    CALL scatter_omp(Var_tmp,VarOut)

  END SUBROUTINE scatter2D_r


  SUBROUTINE scatter2D_r1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      REAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_r1
  

  SUBROUTINE scatter2D_r2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2),SIZE(VarOut,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      REAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_r2  


  SUBROUTINE scatter2D_r3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2),SIZE(VarOut,3),SIZE(VarOut,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      REAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_r3
    
    
!!!!! --> Les booleens


  SUBROUTINE scatter2D_l(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    LOGICAL,DIMENSION(klon_mpi) :: Var_tmp    

!$OMP MASTER    
    CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
    CALL scatter_omp(Var_tmp,VarOut)

  END SUBROUTINE scatter2D_l


  SUBROUTINE scatter2D_l1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      LOGICAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_l1
  

  SUBROUTINE scatter2D_l2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2),SIZE(VarOut,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      LOGICAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_l2  


  SUBROUTINE scatter2D_l3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut

    CALL body(VarIn,VarOut,SIZE(VarOut,2),SIZE(VarOut,3),SIZE(VarOut,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      LOGICAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
!$OMP MASTER    
      CALL scatter2D_mpi(VarIn,Var_tmp)
!$OMP END MASTER
      CALL scatter_omp(Var_tmp,VarOut)

    END SUBROUTINE body

  END SUBROUTINE scatter2D_l3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Gather2D   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! --> Les entiers

  SUBROUTINE gather2D_i(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER,DIMENSION(klon_mpi) :: Var_tmp

    CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
    CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

  END SUBROUTINE gather2D_i
  

  SUBROUTINE gather2D_i1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      INTEGER,DIMENSION(klon_mpi,s1) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_i1

  
  SUBROUTINE gather2D_i2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      INTEGER,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_i2


  SUBROUTINE gather2D_i3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      INTEGER,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_i3


!!!!! --> Les reels

  SUBROUTINE gather2D_r(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL,DIMENSION(klon_mpi) :: Var_tmp

    CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
    CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

  END SUBROUTINE gather2D_r
  

  SUBROUTINE gather2D_r1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      REAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_r1

  
  SUBROUTINE gather2D_r2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      REAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_r2


  SUBROUTINE gather2D_r3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      REAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_r3
  

!!!!! --> Les booleens

  SUBROUTINE gather2D_l(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    LOGICAL,DIMENSION(klon_mpi) :: Var_tmp

    CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
    CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

  END SUBROUTINE gather2D_l
  

  SUBROUTINE gather2D_l1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      LOGICAL,DIMENSION(klon_mpi,s1) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_l1

  
  SUBROUTINE gather2D_l2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      LOGICAL,DIMENSION(klon_mpi,s1,s2) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_l2


  SUBROUTINE gather2D_l3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      LOGICAL,DIMENSION(klon_mpi,s1,s2,s3) :: Var_tmp
      
      CALL gather_omp(VarIn,Var_tmp)
!$OMP MASTER
      CALL gather2D_mpi(Var_tmp,VarOut)
!$OMP END MASTER    

    END SUBROUTINE body

  END SUBROUTINE gather2D_l3
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des reduce_sum   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Les entiers

  SUBROUTINE reduce_sum_i(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  :: VarIn
    INTEGER,INTENT(OUT) :: VarOut
    
    INTEGER             :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
    CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER
  
  END SUBROUTINE reduce_sum_i  


  SUBROUTINE reduce_sum_i1(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      INTEGER,DIMENSION(s1) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE reduce_sum_i1  


  SUBROUTINE reduce_sum_i2(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1),SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      INTEGER,DIMENSION(s1,s2) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE reduce_sum_i2  
  

  SUBROUTINE reduce_sum_i3(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      INTEGER,DIMENSION(s1,s2,s3) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE reduce_sum_i3  


  SUBROUTINE reduce_sum_i4(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3,s4)
      INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3,s4
      INTEGER,DIMENSION(s1,s2,s3,s4) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
   
  END SUBROUTINE reduce_sum_i4  


! Les reels

  SUBROUTINE reduce_sum_r(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN)  :: VarIn
    REAL,INTENT(OUT) :: VarOut
    
    REAL             :: Var_tmp
           
    CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
    CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER
  
  END SUBROUTINE reduce_sum_r  


  SUBROUTINE reduce_sum_r1(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1)
      REAL,INTENT(IN),DIMENSION(:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:) :: VarOut
      INTEGER,INTENT(IN) :: s1
      REAL,DIMENSION(s1) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE reduce_sum_r1  


  SUBROUTINE reduce_sum_r2(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1),SIZE(VarIn,2))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2)
      REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2
      REAL,DIMENSION(s1,s2) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE reduce_sum_r2  
  

  SUBROUTINE reduce_sum_r3(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
      REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3
      REAL,DIMENSION(s1,s2,s3) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
  
  END SUBROUTINE reduce_sum_r3  


  SUBROUTINE reduce_sum_r4(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,SIZE(VarIn,1),SIZE(VarIn,2),SIZE(VarIn,3),SIZE(VarIn,4))
    
  CONTAINS
    SUBROUTINE body(VarIn,VarOut,s1,s2,s3,s4)
      REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
      REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
      INTEGER,INTENT(IN) :: s1,s2,s3,s4
      REAL,DIMENSION(s1,s2,s3,s4) :: Var_tmp
      
      CALL reduce_sum_omp(VarIn,Var_tmp)
!$OMP MASTER      
      CALL reduce_sum_mpi(Var_tmp,VarOut)
!$OMP END MASTER

    END SUBROUTINE body
   
  END SUBROUTINE reduce_sum_r4  

   
END MODULE mod_phys_lmdz_transfert_para

