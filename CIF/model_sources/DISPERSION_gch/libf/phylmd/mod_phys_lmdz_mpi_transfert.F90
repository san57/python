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

MODULE mod_phys_lmdz_mpi_transfert


  INTERFACE bcast_mpi
    MODULE PROCEDURE bcast_mpi_c,                                                     &
                     bcast_mpi_i,bcast_mpi_i1,bcast_mpi_i2,bcast_mpi_i3,bcast_mpi_i4, &
                     bcast_mpi_r,bcast_mpi_r1,bcast_mpi_r2,bcast_mpi_r3,bcast_mpi_r4, &
		     bcast_mpi_l,bcast_mpi_l1,bcast_mpi_l2,bcast_mpi_l3,bcast_mpi_l4
  END INTERFACE

  INTERFACE scatter_mpi
    MODULE PROCEDURE scatter_mpi_i,scatter_mpi_i1,scatter_mpi_i2,scatter_mpi_i3, &
                     scatter_mpi_r,scatter_mpi_r1,scatter_mpi_r2,scatter_mpi_r3, &
		     scatter_mpi_l,scatter_mpi_l1,scatter_mpi_l2,scatter_mpi_l3
  END INTERFACE

  
  INTERFACE gather_mpi
    MODULE PROCEDURE gather_mpi_i,gather_mpi_i1,gather_mpi_i2,gather_mpi_i3, &
                     gather_mpi_r,gather_mpi_r1,gather_mpi_r2,gather_mpi_r3, &
		     gather_mpi_l,gather_mpi_l1,gather_mpi_l2,gather_mpi_l3  
  END INTERFACE
  
  INTERFACE scatter2D_mpi
    MODULE PROCEDURE scatter2D_mpi_i,scatter2D_mpi_i1,scatter2D_mpi_i2,scatter2D_mpi_i3, &
                     scatter2D_mpi_r,scatter2D_mpi_r1,scatter2D_mpi_r2,scatter2D_mpi_r3, &
		     scatter2D_mpi_l,scatter2D_mpi_l1,scatter2D_mpi_l2,scatter2D_mpi_l3
  END INTERFACE

  INTERFACE gather2D_mpi
    MODULE PROCEDURE gather2D_mpi_i,gather2D_mpi_i1,gather2D_mpi_i2,gather2D_mpi_i3, &
                     gather2D_mpi_r,gather2D_mpi_r1,gather2D_mpi_r2,gather2D_mpi_r3, &
		     gather2D_mpi_l,gather2D_mpi_l1,gather2D_mpi_l2,gather2D_mpi_l3
  END INTERFACE 
  
  INTERFACE reduce_sum_mpi
    MODULE PROCEDURE reduce_sum_mpi_i,reduce_sum_mpi_i1,reduce_sum_mpi_i2,reduce_sum_mpi_i3,reduce_sum_mpi_i4, &
                     reduce_sum_mpi_r,reduce_sum_mpi_r1,reduce_sum_mpi_r2,reduce_sum_mpi_r3,reduce_sum_mpi_r4
  END INTERFACE 

 INTERFACE grid1dTo2d_mpi
    MODULE PROCEDURE grid1dTo2d_mpi_i,grid1dTo2d_mpi_i1,grid1dTo2d_mpi_i2,grid1dTo2d_mpi_i3, &
                     grid1dTo2d_mpi_r,grid1dTo2d_mpi_r1,grid1dTo2d_mpi_r2,grid1dTo2d_mpi_r3, &
		     grid1dTo2d_mpi_l,grid1dTo2d_mpi_l1,grid1dTo2d_mpi_l2,grid1dTo2d_mpi_l3
 END INTERFACE 

 INTERFACE grid2dTo1d_mpi
    MODULE PROCEDURE grid2dTo1d_mpi_i,grid2dTo1d_mpi_i1,grid2dTo1d_mpi_i2,grid2dTo1d_mpi_i3, &
                     grid2dTo1d_mpi_r,grid2dTo1d_mpi_r1,grid2dTo1d_mpi_r2,grid2dTo1d_mpi_r3, &
		     grid2dTo1d_mpi_l,grid2dTo1d_mpi_l1,grid2dTo1d_mpi_l2,grid2dTo1d_mpi_l3
 END INTERFACE 
    
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Broadcast --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Les chaine de charactere -- !!

  SUBROUTINE bcast_mpi_c(var1)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var1
   
#ifndef CPP_PARA
    RETURN
#endif

    CALL bcast_mpi_cgen(Var1,len(Var1))
  END SUBROUTINE bcast_mpi_c

!! -- Les entiers -- !!
  
  SUBROUTINE bcast_mpi_i(var1)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var1
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_igen(Var1,1)
  END SUBROUTINE bcast_mpi_i

  SUBROUTINE bcast_mpi_i1(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_igen(Var,size(Var))
  END SUBROUTINE bcast_mpi_i1

  SUBROUTINE bcast_mpi_i2(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_igen(Var,size(Var))
  END SUBROUTINE bcast_mpi_i2

  SUBROUTINE bcast_mpi_i3(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_igen(Var,size(Var))
  END SUBROUTINE bcast_mpi_i3

  SUBROUTINE bcast_mpi_i4(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_igen(Var,size(Var))
  END SUBROUTINE bcast_mpi_i4


!! -- Les reels -- !!

  SUBROUTINE bcast_mpi_r(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_rgen(Var,1)
  END SUBROUTINE bcast_mpi_r

  SUBROUTINE bcast_mpi_r1(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_rgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_r1

  SUBROUTINE bcast_mpi_r2(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_rgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_r2

  SUBROUTINE bcast_mpi_r3(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_rgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_r3

  SUBROUTINE bcast_mpi_r4(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_rgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_r4
  
!! -- Les booleans -- !!

  SUBROUTINE bcast_mpi_l(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_lgen(Var,1)
  END SUBROUTINE bcast_mpi_l

  SUBROUTINE bcast_mpi_l1(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_lgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_l1

  SUBROUTINE bcast_mpi_l2(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_lgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_l2

  SUBROUTINE bcast_mpi_l3(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_lgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_l3

  SUBROUTINE bcast_mpi_l4(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:,:)
   
#ifndef CPP_PARA
    RETURN
#endif
    CALL bcast_mpi_lgen(Var,size(Var))
  END SUBROUTINE bcast_mpi_l4
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter_mpi_i(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    INTEGER :: dummy

#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif

     IF (is_mpi_root) THEN
      CALL scatter_mpi_igen(VarIn,Varout,1)
     ELSE
      CALL scatter_mpi_igen(dummy,Varout,1)
    ENDIF
    
  END SUBROUTINE scatter_mpi_i

  SUBROUTINE scatter_mpi_i1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER :: dummy

#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_igen(VarIn,Varout,Size(VarOut,2))
    ELSE
      CALL scatter_mpi_igen(dummy,Varout,Size(VarOut,2))
    ENDIF
    
  END SUBROUTINE scatter_mpi_i1
  
  SUBROUTINE scatter_mpi_i2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3))
    ELSE
      CALL scatter_mpi_igen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3))
    ENDIF
  END SUBROUTINE scatter_mpi_i2

  SUBROUTINE scatter_mpi_i3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ELSE
      CALL scatter_mpi_igen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ENDIF
  
  END SUBROUTINE scatter_mpi_i3


  SUBROUTINE scatter_mpi_r(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_rgen(VarIn,Varout,1)
    ELSE
      CALL scatter_mpi_rgen(dummy,Varout,1)
    ENDIF
  
  END SUBROUTINE scatter_mpi_r

  SUBROUTINE scatter_mpi_r1(VarIn, VarOut)
  USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
  IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_rgen(VarIn,Varout,Size(VarOut,2))
    ELSE
      CALL scatter_mpi_rgen(dummy,Varout,Size(VarOut,2))      
    ENDIF
  
  END SUBROUTINE scatter_mpi_r1
  
  SUBROUTINE scatter_mpi_r2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3))
    ELSE
      CALL scatter_mpi_rgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3))
    ENDIF
  
  END SUBROUTINE scatter_mpi_r2

  SUBROUTINE scatter_mpi_r3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ELSE
      CALL scatter_mpi_rgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ENDIF
  
  END SUBROUTINE scatter_mpi_r3


  SUBROUTINE scatter_mpi_l(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_lgen(VarIn,Varout,1)
    ELSE
      CALL scatter_mpi_lgen(dummy,Varout,1)
    ENDIF
    
  END SUBROUTINE scatter_mpi_l

  SUBROUTINE scatter_mpi_l1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_lgen(VarIn,Varout,Size(VarOut,2))
    ELSE
      CALL scatter_mpi_lgen(dummy,Varout,Size(VarOut,2))      
    ENDIF
  
  END SUBROUTINE scatter_mpi_l1
  
  SUBROUTINE scatter_mpi_l2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3))
    ELSE
      CALL scatter_mpi_lgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3))
    ENDIF
  
  END SUBROUTINE scatter_mpi_l2

  SUBROUTINE scatter_mpi_l3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif
    IF (is_mpi_root) THEN
      CALL scatter_mpi_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ELSE
      CALL scatter_mpi_lgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4))
    ENDIF
  
  END SUBROUTINE scatter_mpi_l3  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Gather   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!! --> Les entiers

  SUBROUTINE gather_mpi_i(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    INTEGER :: dummy

#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_igen(VarIn,VarOut,1)
    ELSE
      CALL gather_mpi_igen(VarIn,dummy,1)
    ENDIF
  
  END SUBROUTINE gather_mpi_i
  
  
  
  

!!!!!

  SUBROUTINE gather_mpi_i1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_igen(VarIn,VarOut,Size(VarIn,2))
    ELSE
      CALL gather_mpi_igen(VarIn,dummy,Size(VarIn,2))
    ENDIF
  
  END SUBROUTINE gather_mpi_i1

!!!!!
  
  SUBROUTINE gather_mpi_i2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_igen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3))
    ELSE
      CALL gather_mpi_igen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3))
    ENDIF
  
  END SUBROUTINE gather_mpi_i2

!!!!!

  SUBROUTINE gather_mpi_i3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_igen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ELSE
      CALL gather_mpi_igen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ENDIF
  
  END SUBROUTINE gather_mpi_i3

!!!!! --> Les reels

  SUBROUTINE gather_mpi_r(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_rgen(VarIn,VarOut,1)
    ELSE
      CALL gather_mpi_rgen(VarIn,dummy,1)
    ENDIF
  
  END SUBROUTINE gather_mpi_r

!!!!!

  SUBROUTINE gather_mpi_r1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_rgen(VarIn,VarOut,Size(VarIn,2))
    ELSE
      CALL gather_mpi_rgen(VarIn,dummy,Size(VarIn,2))
    ENDIF
  
  END SUBROUTINE gather_mpi_r1

!!!!!
  
  SUBROUTINE gather_mpi_r2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_rgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3))
    ELSE
      CALL gather_mpi_rgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3))      
    ENDIF
  
  END SUBROUTINE gather_mpi_r2

!!!!!

  SUBROUTINE gather_mpi_r3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_rgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ELSE
      CALL gather_mpi_rgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ENDIF
  
  END SUBROUTINE gather_mpi_r3

!!!!! --> Les booleen

  SUBROUTINE gather_mpi_l(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_lgen(VarIn,VarOut,1)
    ELSE
      CALL gather_mpi_lgen(VarIn,dummy,1)      
    ENDIF
  
  END SUBROUTINE gather_mpi_l

!!!!!

  SUBROUTINE gather_mpi_l1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_lgen(VarIn,VarOut,Size(VarIn,2))
    ELSE
      CALL gather_mpi_lgen(VarIn,dummy,Size(VarIn,2))
    ENDIF
  
  END SUBROUTINE gather_mpi_l1

!!!!!
  
  SUBROUTINE gather_mpi_l2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_lgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3))
    ELSE
      CALL gather_mpi_lgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3))
    ENDIF
  
  END SUBROUTINE gather_mpi_l2

!!!!!

  SUBROUTINE gather_mpi_l3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    LOGICAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL gather_mpi_lgen(VarIn,VarOut,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))
    ELSE
      CALL gather_mpi_lgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4))      
    ENDIF
  
  END SUBROUTINE gather_mpi_l3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter2D   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter2D_mpi_i(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    INTEGER,DIMENSION(klon_glo) :: Var_tmp    
    
    CALL grid2dTo1d_glo(VarIn,Var_tmp)
    CALL scatter_mpi(Var_tmp,VarOut)

  END SUBROUTINE scatter2D_mpi_i

  SUBROUTINE scatter2D_mpi_i1(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1)
        INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
        INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1
	
	INTEGER,DIMENSION(klon_glo,s1) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE scatter2D_mpi_i1

  SUBROUTINE scatter2D_mpi_i2(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2),size(VarOut,3))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2)
        INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
        INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2
	
        INTEGER,DIMENSION(klon_glo,s1,s2) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE scatter2D_mpi_i2
  
  SUBROUTINE scatter2D_mpi_i3(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2),size(VarOut,3),size(VarOut,4))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
        INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
        INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2,s3
	
        INTEGER,DIMENSION(klon_glo,s1,s2,s3) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body
  
 
  END SUBROUTINE scatter2D_mpi_i3



  SUBROUTINE scatter2D_mpi_r(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    REAL,DIMENSION(klon_glo) :: Var_tmp    
    
    CALL grid2dTo1d_glo(VarIn,Var_tmp)
    CALL scatter_mpi(Var_tmp,VarOut)

  END SUBROUTINE scatter2D_mpi_R


  SUBROUTINE scatter2D_mpi_r1(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut

    CALL body(VarIn,VarOut,size(VarOut,2))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1)
        REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
        REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1
	
        REAL,DIMENSION(klon_glo,s1) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE scatter2D_mpi_r1


  SUBROUTINE scatter2D_mpi_r2(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2),size(VarOut,3))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2)
        REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
        REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2
	
        REAL,DIMENSION(klon_glo,s1,s2) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE scatter2D_mpi_r2
  
  SUBROUTINE scatter2D_mpi_r3(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2),size(VarOut,3),size(VarOut,4))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
        REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
        REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2,s3
	
        REAL,DIMENSION(klon_glo,s1,s2,s3) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body
 
  END SUBROUTINE scatter2D_mpi_r3
  
  
  SUBROUTINE scatter2D_mpi_l(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    LOGICAL,DIMENSION(klon_glo) :: Var_tmp    
    
    CALL grid2dTo1d_glo(VarIn,Var_tmp)
    CALL scatter_mpi(Var_tmp,VarOut)

  END SUBROUTINE scatter2D_mpi_l


  SUBROUTINE scatter2D_mpi_l1(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1)
        LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
        LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1
	
        LOGICAL,DIMENSION(klon_glo,s1) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body
  
  END SUBROUTINE scatter2D_mpi_l1


  SUBROUTINE scatter2D_mpi_l2(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2),size(VarOut,3))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2)
        LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
        LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2
	
        LOGICAL,DIMENSION(klon_glo,s1,s2) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE scatter2D_mpi_l2
  
  SUBROUTINE scatter2D_mpi_l3(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,2),size(VarOut,3),size(VarOut,4))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
        LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
        LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2,s3
	
        LOGICAL,DIMENSION(klon_glo,s1,s2,s3) :: Var_tmp
        
        CALL grid2dTo1d_glo(VarIn,Var_tmp)
        CALL scatter_mpi(Var_tmp,VarOut)
      END SUBROUTINE body
 
  END SUBROUTINE scatter2D_mpi_l3
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Gather2D   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE gather2D_mpi_i(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER,DIMENSION(klon_glo) :: Var_tmp
    
    CALL gather_mpi(VarIn,Var_tmp)
    CALL grid1dTo2d_glo(Var_tmp,VarOut)

  END SUBROUTINE gather2D_mpi_i

  SUBROUTINE gather2D_mpi_i1(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1)
        INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
        INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1
	
        INTEGER,DIMENSION(klon_glo,s1) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_i1

  SUBROUTINE gather2D_mpi_i2(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3),SIZE(VarOut,4))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2)
        INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
        INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2
	
        INTEGER,DIMENSION(klon_glo,s1,s2) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_i2
  
  SUBROUTINE gather2D_mpi_i3(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3),SIZE(VarOut,4),SIZE(VarOut,5))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
        INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
        INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2,s3
	
        INTEGER,DIMENSION(klon_glo,s1,s2,s3) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_i3



  SUBROUTINE gather2D_mpi_r(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL,DIMENSION(klon_glo) :: Var_tmp
    
    CALL gather_mpi(VarIn,Var_tmp)
    CALL grid1dTo2d_glo(Var_tmp,VarOut)

  END SUBROUTINE gather2D_mpi_r

  SUBROUTINE gather2D_mpi_r1(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1)
        REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
        REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1
	
        REAL,DIMENSION(klon_glo,s1) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_r1

  SUBROUTINE gather2D_mpi_r2(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3),SIZE(VarOut,4))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2)
        REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
        REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2
	
        REAL,DIMENSION(klon_glo,s1,s2) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_r2
  
  SUBROUTINE gather2D_mpi_r3(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3),SIZE(VarOut,4),SIZE(VarOut,5))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
        REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
        REAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2,s3
	
        REAL,DIMENSION(klon_glo,s1,s2,s3) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_r3

  
  
  SUBROUTINE gather2D_mpi_l(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    LOGICAL,DIMENSION(klon_glo) :: Var_tmp
    
    CALL gather_mpi(VarIn,Var_tmp)
    CALL grid1dTo2d_glo(Var_tmp,VarOut)

  END SUBROUTINE gather2D_mpi_l

  SUBROUTINE gather2D_mpi_l1(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1)
        LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
        LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1
	
        LOGICAL,DIMENSION(klon_glo,s1) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_l1

  SUBROUTINE gather2D_mpi_l2(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3),SIZE(VarOut,4))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2)
        LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
        LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2
	
        LOGICAL,DIMENSION(klon_glo,s1,s2) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_l2
  
  SUBROUTINE gather2D_mpi_l3(VarIn, VarOut)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
    
    CALL body(VarIn,VarOut,size(VarOut,3),SIZE(VarOut,4),SIZE(VarOut,5))
  
    CONTAINS 
      SUBROUTINE body(VarIn,VarOut,s1,s2,s3)
        LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
        LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:) :: VarOut
        INTEGER,INTENT(IN) :: s1,s2,s3
	
        LOGICAL,DIMENSION(klon_glo,s1,s2,s3) :: Var_tmp
        
        CALL gather_mpi(VarIn,Var_tmp)
        CALL grid1dTo2d_glo(Var_tmp,VarOut)
      END SUBROUTINE body

  END SUBROUTINE gather2D_mpi_l3
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des reduce_sum   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE reduce_sum_mpi_i(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  :: VarIn
    INTEGER,INTENT(OUT) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut=VarIn
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_igen(VarIn,Varout,1)
    ELSE
      CALL reduce_sum_mpi_igen(VarIn,dummy,1)
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_i

  SUBROUTINE reduce_sum_mpi_i1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_i1

  SUBROUTINE reduce_sum_mpi_i2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_i2

  SUBROUTINE reduce_sum_mpi_i3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_i3

  SUBROUTINE reduce_sum_mpi_i4(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    INTEGER :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_igen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_igen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_i4                  
  
  
  SUBROUTINE reduce_sum_mpi_r(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN)  :: VarIn
    REAL,INTENT(OUT) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut=VarIn
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_rgen(VarIn,Varout,1)
    ELSE
      CALL reduce_sum_mpi_rgen(VarIn,dummy,1)
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_r

  SUBROUTINE reduce_sum_mpi_r1(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:)=VarIn(:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_r1

  SUBROUTINE reduce_sum_mpi_r2(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:)=VarIn(:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_r2

  SUBROUTINE reduce_sum_mpi_r3(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:)=VarIn(:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_r3

  SUBROUTINE reduce_sum_mpi_r4(VarIn, VarOut)
    USE mod_phys_lmdz_mpi_data, ONLY :  is_mpi_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    REAL :: dummy
    
#ifndef CPP_PARA
    VarOut(:,:,:,:)=VarIn(:,:,:,:)
    RETURN
#endif

    IF (is_mpi_root) THEN
      CALL reduce_sum_mpi_rgen(VarIn,Varout,SIZE(VarIn))
    ELSE
      CALL reduce_sum_mpi_rgen(VarIn,dummy,SIZE(VarIn))      
    ENDIF
  
  END SUBROUTINE reduce_sum_mpi_r4 
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE grid1dTo2d  !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE grid1dTo2d_mpi_i(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_igen(VarIn,VarOut,1)
  
  END SUBROUTINE grid1dTo2d_mpi_i
  

  SUBROUTINE grid1dTo2d_mpi_i1(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_igen(VarIn,VarOut,size(VarIn,2))
  
  END SUBROUTINE grid1dTo2d_mpi_i1

  SUBROUTINE grid1dTo2d_mpi_i2(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_igen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3))
  
  END SUBROUTINE grid1dTo2d_mpi_i2
  
  SUBROUTINE grid1dTo2d_mpi_i3(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_igen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid1dTo2d_mpi_i3


  SUBROUTINE grid1dTo2d_mpi_r(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_rgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid1dTo2d_mpi_r
  

  SUBROUTINE grid1dTo2d_mpi_r1(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_rgen(VarIn,VarOut,size(VarIn,2))
  
  END SUBROUTINE grid1dTo2d_mpi_r1

  SUBROUTINE grid1dTo2d_mpi_r2(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_rgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3))
  
  END SUBROUTINE grid1dTo2d_mpi_r2
  
  SUBROUTINE grid1dTo2d_mpi_r3(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:,:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_rgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid1dTo2d_mpi_r3
  
  
  
  SUBROUTINE grid1dTo2d_mpi_l(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_lgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid1dTo2d_mpi_l
  

  SUBROUTINE grid1dTo2d_mpi_l1(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_lgen(VarIn,VarOut,size(VarIn,2))
  
  END SUBROUTINE grid1dTo2d_mpi_l1

  SUBROUTINE grid1dTo2d_mpi_l2(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_lgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3))
  
  END SUBROUTINE grid1dTo2d_mpi_l2
  
  SUBROUTINE grid1dTo2d_mpi_l3(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_mpi_lgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid1dTo2d_mpi_l3


  SUBROUTINE grid2dTo1d_mpi_i(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:)  :: VarOut
    
    CALL grid2dTo1d_mpi_igen(VarIn,VarOut,1)
  
  END SUBROUTINE grid2dTo1d_mpi_i
  

  SUBROUTINE grid2dTo1d_mpi_i1(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_igen(VarIn,VarOut,size(VarIn,3))
  
  END SUBROUTINE grid2dTo1d_mpi_i1

  SUBROUTINE grid2dTo1d_mpi_i2(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_igen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid2dTo1d_mpi_i2
  
  SUBROUTINE grid2dTo1d_mpi_i3(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_igen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4)*size(VarIn,5))
  
  END SUBROUTINE grid2dTo1d_mpi_i3
 



  SUBROUTINE grid2dTo1d_mpi_r(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:)  :: VarOut
    
    CALL grid2dTo1d_mpi_rgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid2dTo1d_mpi_r
  

  SUBROUTINE grid2dTo1d_mpi_r1(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_rgen(VarIn,VarOut,size(VarIn,3))
  
  END SUBROUTINE grid2dTo1d_mpi_r1

  SUBROUTINE grid2dTo1d_mpi_r2(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_rgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid2dTo1d_mpi_r2
  
  SUBROUTINE grid2dTo1d_mpi_r3(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_rgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4)*size(VarIn,5))
  
  END SUBROUTINE grid2dTo1d_mpi_r3



  SUBROUTINE grid2dTo1d_mpi_l(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:)  :: VarOut
    
    CALL grid2dTo1d_mpi_lgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid2dTo1d_mpi_l
  

  SUBROUTINE grid2dTo1d_mpi_l1(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_lgen(VarIn,VarOut,size(VarIn,3))
  
  END SUBROUTINE grid2dTo1d_mpi_l1



  SUBROUTINE grid2dTo1d_mpi_l2(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_lgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid2dTo1d_mpi_l2

  
  SUBROUTINE grid2dTo1d_mpi_l3(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid2dTo1d_mpi_lgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4)*size(VarIn,5))
  
  END SUBROUTINE grid2dTo1d_mpi_l3


                            
END MODULE mod_phys_lmdz_mpi_transfert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEFINITION DES FONCTIONS DE TRANSFERT GENERIQUES !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE bcast_mpi_cgen(var,nb)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    IMPLICIT NONE
    
    CHARACTER(LEN=*),INTENT(INOUT) :: Var
    INTEGER,INTENT(IN) :: nb
    
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    INTEGER :: ierr

    IF (.not.is_ok_mpi) RETURN
    
#ifdef CPP_PARA
    CALL MPI_BCAST(Var,nb,MPI_CHARACTER,mpi_root_x,COMM_LMDZ_PHY,ierr)
#endif
        
  END SUBROUTINE bcast_mpi_cgen


      
  SUBROUTINE bcast_mpi_igen(var,nb)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    IMPLICIT NONE
    
    INTEGER,DIMENSION(nb),INTENT(INOUT) :: Var
    INTEGER,INTENT(IN) :: nb
    
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    INTEGER :: ierr

    IF (.not.is_ok_mpi) RETURN

#ifdef CPP_PARA
    CALL MPI_BCAST(Var,nb,MPI_INTEGER,mpi_root_x,COMM_LMDZ_PHY,ierr)
#endif
        
  END SUBROUTINE bcast_mpi_igen



  
  SUBROUTINE bcast_mpi_rgen(var,nb)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    IMPLICIT NONE
    
    REAL,DIMENSION(nb),INTENT(INOUT) :: Var
    INTEGER,INTENT(IN) :: nb
    
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    INTEGER :: ierr

    IF (.not.is_ok_mpi) RETURN

#ifdef CPP_PARA
    CALL MPI_BCAST(Var,nb,MPI_REAL_LMDZ,mpi_root_x,COMM_LMDZ_PHY,ierr)
#endif
    
  END SUBROUTINE bcast_mpi_rgen
  



  SUBROUTINE bcast_mpi_lgen(var,nb)
    USE mod_phys_lmdz_mpi_data ,  mpi_root_x=>mpi_root
    IMPLICIT NONE
    
    LOGICAL,DIMENSION(nb),INTENT(INOUT) :: Var
    INTEGER,INTENT(IN) :: nb
    
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    INTEGER :: ierr

    IF (.not.is_ok_mpi) RETURN

#ifdef CPP_PARA
    CALL MPI_BCAST(Var,nb,MPI_LOGICAL,mpi_root_x,COMM_LMDZ_PHY,ierr)
#endif

  END SUBROUTINE bcast_mpi_lgen

  

  SUBROUTINE scatter_mpi_igen(VarIn, VarOut, dimsize)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(klon_glo,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
  
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    INTEGER,DIMENSION(dimsize*klon_glo) :: VarTmp
    INTEGER :: nb,i,index,rank
    INTEGER :: ierr


    IF (.not.is_ok_mpi) THEN
      VarOut(:,:)=VarIn(:,:)
      RETURN
    ENDIF

    
    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        displs(rank)=Index-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp(Index:Index+nb-1)=VarIn(klon_mpi_para_begin(rank):klon_mpi_para_end(rank),i)
          Index=Index+nb
        ENDDO
      ENDDO
    ENDIF
      
#ifdef CPP_PARA 
    CALL MPI_SCATTERV(VarTmp,counts,displs,MPI_INTEGER,VarOut,klon_mpi*dimsize,   &
                      MPI_INTEGER,mpi_root_x, COMM_LMDZ_PHY,ierr)
#endif

  END SUBROUTINE scatter_mpi_igen

  SUBROUTINE scatter_mpi_rgen(VarIn, VarOut, dimsize)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(klon_glo,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
  
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    REAL,DIMENSION(dimsize*klon_glo) :: VarTmp
    INTEGER :: nb,i,index,rank
    INTEGER :: ierr

    IF (.not.is_ok_mpi) THEN
      VarOut(:,:)=VarIn(:,:)
      RETURN
    ENDIF
    
    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        displs(rank)=Index-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp(Index:Index+nb-1)=VarIn(klon_mpi_para_begin(rank):klon_mpi_para_end(rank),i)
          Index=Index+nb
        ENDDO
      ENDDO
    ENDIF
      
#ifdef CPP_PARA 
    CALL MPI_SCATTERV(VarTmp,counts,displs,MPI_REAL_LMDZ,VarOut,klon_mpi*dimsize,   &
                      MPI_REAL_LMDZ,mpi_root_x, COMM_LMDZ_PHY,ierr)

#endif

  END SUBROUTINE scatter_mpi_rgen

  
  SUBROUTINE scatter_mpi_lgen(VarIn, VarOut, dimsize)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(klon_glo,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
  
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    LOGICAL,DIMENSION(dimsize*klon_glo) :: VarTmp
    INTEGER :: nb,i,index,rank
    INTEGER :: ierr

    IF (.not.is_ok_mpi) THEN
      VarOut(:,:)=VarIn(:,:)
      RETURN
    ENDIF
    
    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        displs(rank)=Index-1
        counts(rank)=nb*dimsize
        DO i=1,dimsize
          VarTmp(Index:Index+nb-1)=VarIn(klon_mpi_para_begin(rank):klon_mpi_para_end(rank),i)
          Index=Index+nb
        ENDDO
      ENDDO
    ENDIF
      
#ifdef CPP_PARA
    CALL MPI_SCATTERV(VarTmp,counts,displs,MPI_LOGICAL,VarOut,klon_mpi*dimsize,   &
                      MPI_LOGICAL,mpi_root_x, COMM_LMDZ_PHY,ierr)
#endif

  END SUBROUTINE scatter_mpi_lgen  




  SUBROUTINE gather_mpi_igen(VarIn, VarOut, dimsize)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(klon_glo,dimsize) :: VarOut
  
    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    INTEGER,DIMENSION(dimsize*klon_glo) :: VarTmp
    INTEGER :: nb,i,index,rank
    INTEGER :: ierr

    IF (.not.is_ok_mpi) THEN
      VarOut(:,:)=VarIn(:,:)
      RETURN
    ENDIF

    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        displs(rank)=Index-1
        counts(rank)=nb*dimsize
	Index=Index+nb*dimsize
      ENDDO
     
    ENDIF
    
#ifdef CPP_PARA
    CALL MPI_GATHERV(VarIn,klon_mpi*dimsize,MPI_INTEGER,VarTmp,counts,displs,   &
                     MPI_INTEGER,mpi_root_x, COMM_LMDZ_PHY,ierr)
#endif

		          
    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        DO i=1,dimsize
          VarOut(klon_mpi_para_begin(rank):klon_mpi_para_end(rank),i)=VarTmp(Index:Index+nb-1)
	  Index=Index+nb
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE gather_mpi_igen  

  SUBROUTINE gather_mpi_rgen(VarIn, VarOut, dimsize)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(klon_glo,dimsize) :: VarOut
  
    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    REAL,DIMENSION(dimsize*klon_glo) :: VarTmp
    INTEGER :: nb,i,index,rank
    INTEGER :: ierr

    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        displs(rank)=Index-1
        counts(rank)=nb*dimsize
	Index=Index+nb*dimsize
      ENDDO
    ENDIF
    
    IF (.not.is_ok_mpi) THEN
      VarOut(:,:)=VarIn(:,:)
      RETURN
    ENDIF

#ifdef CPP_PARA
    CALL MPI_GATHERV(VarIn,klon_mpi*dimsize,MPI_REAL_LMDZ,VarTmp,counts,displs,   &
                      MPI_REAL_LMDZ,mpi_root_x, COMM_LMDZ_PHY,ierr)
#endif
		          
    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        DO i=1,dimsize
          VarOut(klon_mpi_para_begin(rank):klon_mpi_para_end(rank),i)=VarTmp(Index:Index+nb-1)
	  Index=Index+nb
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE gather_mpi_rgen  

  SUBROUTINE gather_mpi_lgen(VarIn, VarOut, dimsize)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
  
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(klon_glo,dimsize) :: VarOut
  
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif

    INTEGER,DIMENSION(0:mpi_size-1) :: displs
    INTEGER,DIMENSION(0:mpi_size-1) :: counts
    LOGICAL,DIMENSION(dimsize*klon_glo) :: VarTmp
    INTEGER :: nb,i,index,rank
    INTEGER :: ierr
    
    IF (.not.is_ok_mpi) THEN
      VarOut(:,:)=VarIn(:,:)
      RETURN
    ENDIF

    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        displs(rank)=Index-1
        counts(rank)=nb*dimsize
	Index=Index+nb*dimsize
      ENDDO
    ENDIF
    

#ifdef CPP_PARA
    CALL MPI_GATHERV(VarIn,klon_mpi*dimsize,MPI_LOGICAL,VarTmp,counts,displs,   &
                      MPI_LOGICAL,mpi_root_x, COMM_LMDZ_PHY,ierr)
#endif
		          
    IF (is_mpi_root) THEN
      Index=1
      DO rank=0,mpi_size-1
        nb=klon_mpi_para_nb(rank)
        DO i=1,dimsize
          VarOut(klon_mpi_para_begin(rank):klon_mpi_para_end(rank),i)=VarTmp(Index:Index+nb-1)
	  Index=Index+nb
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE gather_mpi_lgen
  


  SUBROUTINE reduce_sum_mpi_igen(VarIn,VarOut,nb)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
   
    INTEGER,DIMENSION(nb),INTENT(IN) :: VarIn
    INTEGER,DIMENSION(nb),INTENT(OUT) :: VarOut    
    INTEGER,INTENT(IN) :: nb
    INTEGER :: ierr
   
    IF (.not.is_ok_mpi) THEN
      VarOut(:)=VarIn(:)
      RETURN
    ENDIF


#ifdef CPP_PARA
    CALL MPI_REDUCE(VarIn,VarOut,nb,MPI_INTEGER,MPI_SUM,mpi_root_x,COMM_LMDZ_PHY,ierr)
#endif
            
  END SUBROUTINE reduce_sum_mpi_igen
  
  SUBROUTINE reduce_sum_mpi_rgen(VarIn,VarOut,nb)
    USE mod_phys_lmdz_mpi_data , mpi_root_x=>mpi_root
    USE mod_grid_phy_lmdz

    IMPLICIT NONE

#ifdef CPP_PARA
    INCLUDE 'mpif.h'
#endif
    
    REAL,DIMENSION(nb),INTENT(IN) :: VarIn
    REAL,DIMENSION(nb),INTENT(OUT) :: VarOut    
    INTEGER,INTENT(IN) :: nb
    INTEGER :: ierr
 
    IF (.not.is_ok_mpi) THEN
      VarOut(:)=VarIn(:)
      RETURN
    ENDIF
   
#ifdef CPP_PARA
    CALL MPI_REDUCE(VarIn,VarOut,nb,MPI_REAL_LMDZ,MPI_SUM,mpi_root_x,COMM_LMDZ_PHY,ierr)
#endif
        
  END SUBROUTINE reduce_sum_mpi_rgen



  SUBROUTINE grid1dTo2d_mpi_igen(VarIn,VarOut,dimsize)
    USE mod_phys_lmdz_mpi_data
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN) ,DIMENSION(klon_mpi,dimsize)       :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(nbp_lon*jj_nb,dimsize)  :: VarOut
    INTEGER :: i,ij,Offset

    
    VarOut(1:nbp_lon,:)=0
    VarOut(nbp_lon*(jj_nb-1)+1:nbp_lon*jj_nb,:)=0
    
    offset=ii_begin
    IF (is_north_pole) Offset=nbp_lon
    
    
    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij+offset-1,i)=VarIn(ij,i)
      ENDDO
    ENDDO
    
    
    IF (is_north_pole) THEN 
      DO i=1,dimsize
        DO ij=1,nbp_lon
         VarOut(ij,i)=VarIn(1,i)
	ENDDO
      ENDDO
    ENDIF
    
    IF (is_south_pole) THEN 
      DO i=1,dimsize
        DO ij=nbp_lon*(jj_nb-1)+1,nbp_lon*jj_nb
         VarOut(ij,i)=VarIn(klon_mpi,i)
	ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE grid1dTo2d_mpi_igen   


  SUBROUTINE grid1dTo2d_mpi_rgen(VarIn,VarOut,dimsize)
    USE mod_phys_lmdz_mpi_data
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN) ,DIMENSION(klon_mpi,dimsize)       :: VarIn
    REAL,INTENT(OUT),DIMENSION(nbp_lon*jj_nb,dimsize)  :: VarOut
    INTEGER :: i,ij,Offset

    
    VarOut(1:nbp_lon,:)=0
    VarOut(nbp_lon*(jj_nb-1)+1:nbp_lon*jj_nb,:)=0
    
    offset=ii_begin
    IF (is_north_pole) Offset=nbp_lon
    
    
    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij+offset-1,i)=VarIn(ij,i)
      ENDDO
    ENDDO
    
    
    IF (is_north_pole) THEN 
      DO i=1,dimsize
        DO ij=1,nbp_lon
         VarOut(ij,i)=VarIn(1,i)
	ENDDO
      ENDDO
    ENDIF
    
    IF (is_south_pole) THEN 
      DO i=1,dimsize
        DO ij=nbp_lon*(jj_nb-1)+1,nbp_lon*jj_nb
         VarOut(ij,i)=VarIn(klon_mpi,i)
	ENDDO
      ENDDO
    ENDIF

   END SUBROUTINE grid1dTo2d_mpi_rgen   



  SUBROUTINE grid1dTo2d_mpi_lgen(VarIn,VarOut,dimsize)
    USE mod_phys_lmdz_mpi_data
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN) ,DIMENSION(klon_mpi,dimsize)       :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(nbp_lon*jj_nb,dimsize)  :: VarOut
    INTEGER :: i,ij,Offset

    
    VarOut(1:nbp_lon,:)=.FALSE.
    VarOut(nbp_lon*(jj_nb-1)+1:nbp_lon*jj_nb,:)=.FALSE.
    
    offset=ii_begin
    IF (is_north_pole) Offset=nbp_lon
    
    
    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij+offset-1,i)=VarIn(ij,i)
      ENDDO
    ENDDO
    
    
    IF (is_north_pole) THEN 
      DO i=1,dimsize
        DO ij=1,nbp_lon
         VarOut(ij,i)=VarIn(1,i)
	ENDDO
      ENDDO
    ENDIF
    
    IF (is_south_pole) THEN 
      DO i=1,dimsize
        DO ij=nbp_lon*(jj_nb-1)+1,nbp_lon*jj_nb
         VarOut(ij,i)=VarIn(klon_mpi,i)
	ENDDO
      ENDDO
    ENDIF

   END SUBROUTINE grid1dTo2d_mpi_lgen   

  


  SUBROUTINE grid2dTo1d_mpi_igen(VarIn,VarOut,dimsize)
    USE mod_phys_lmdz_mpi_data
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN) ,DIMENSION(nbp_lon*jj_nb,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(klon_mpi,dimsize)      :: VarOut
    INTEGER :: i,ij,offset

    offset=ii_begin
    IF (is_north_pole) offset=nbp_lon

    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij,i)=VarIn(ij+offset-1,i)
      ENDDO
    ENDDO

    IF (is_north_pole) THEN 
      DO i=1,dimsize
        VarOut(1,i)=VarIn(1,i)
      ENDDO
    ENDIF
    
    
  END SUBROUTINE grid2dTo1d_mpi_igen   



  SUBROUTINE grid2dTo1d_mpi_rgen(VarIn,VarOut,dimsize)
    USE mod_phys_lmdz_mpi_data
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN) ,DIMENSION(nbp_lon*jj_nb,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize)      :: VarOut
    INTEGER :: i,ij,offset

    offset=ii_begin
    IF (is_north_pole) offset=nbp_lon

    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij,i)=VarIn(ij+offset-1,i)
      ENDDO
    ENDDO

    IF (is_north_pole) THEN 
      DO i=1,dimsize
         VarOut(1,i)=VarIn(1,i)
      ENDDO
    ENDIF
    
    
  END SUBROUTINE grid2dTo1d_mpi_rgen   
  

  SUBROUTINE grid2dTo1d_mpi_lgen(VarIn,VarOut,dimsize)
    USE mod_phys_lmdz_mpi_data
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN) ,DIMENSION(nbp_lon*jj_nb,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize)      :: VarOut
    INTEGER :: i,ij,offset

    offset=ii_begin
    IF (is_north_pole) offset=nbp_lon

    DO i=1,dimsize
      DO ij=1,klon_mpi
        VarOut(ij,i)=VarIn(ij+offset-1,i)
      ENDDO
    ENDDO

    IF (is_north_pole) THEN 
      DO i=1,dimsize
        VarOut(1,i)=VarIn(1,i)
      ENDDO
    ENDIF
    
    
  END SUBROUTINE grid2dTo1d_mpi_lgen   

