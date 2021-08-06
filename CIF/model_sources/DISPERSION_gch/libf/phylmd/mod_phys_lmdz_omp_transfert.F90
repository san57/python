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

MODULE mod_phys_lmdz_omp_transfert

  INTEGER,PARAMETER :: omp_buffer_size = 1024*1024*16
  INTEGER,SAVE,DIMENSION(omp_buffer_size) :: omp_buffer
  
  INTERFACE bcast_omp
    MODULE PROCEDURE bcast_omp_c,                                                     &
                     bcast_omp_i,bcast_omp_i1,bcast_omp_i2,bcast_omp_i3,bcast_omp_i4, &
                     bcast_omp_r,bcast_omp_r1,bcast_omp_r2,bcast_omp_r3,bcast_omp_r4, &
		     bcast_omp_l,bcast_omp_l1,bcast_omp_l2,bcast_omp_l3,bcast_omp_l4
  END INTERFACE

  INTERFACE scatter_omp
    MODULE PROCEDURE scatter_omp_i,scatter_omp_i1,scatter_omp_i2,scatter_omp_i3, &
                     scatter_omp_r,scatter_omp_r1,scatter_omp_r2,scatter_omp_r3, &
		     scatter_omp_l,scatter_omp_l1,scatter_omp_l2,scatter_omp_l3
  END INTERFACE

  
  INTERFACE gather_omp
    MODULE PROCEDURE gather_omp_i,gather_omp_i1,gather_omp_i2,gather_omp_i3, &
                     gather_omp_r,gather_omp_r1,gather_omp_r2,gather_omp_r3, &
		     gather_omp_l,gather_omp_l1,gather_omp_l2,gather_omp_l3  
  END INTERFACE
  
  
  INTERFACE reduce_sum_omp
    MODULE PROCEDURE reduce_sum_omp_i,reduce_sum_omp_i1,reduce_sum_omp_i2,reduce_sum_omp_i3,reduce_sum_omp_i4, &
                     reduce_sum_omp_r,reduce_sum_omp_r1,reduce_sum_omp_r2,reduce_sum_omp_r3,reduce_sum_omp_r4
  END INTERFACE 

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Broadcast --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! -- Les chaine de charactère -- !!

  SUBROUTINE bcast_omp_c(var)
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(INOUT) :: Var

    CALL bcast_omp_cgen(Var,len(Var),omp_buffer)
    
  END SUBROUTINE bcast_omp_c

!! -- Les entiers -- !!
  
  SUBROUTINE bcast_omp_i(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var

    CALL bcast_omp_igen(Var,1,omp_buffer)

  END SUBROUTINE bcast_omp_i


  SUBROUTINE bcast_omp_i1(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:)
   
    CALL bcast_omp_igen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_i1


  SUBROUTINE bcast_omp_i2(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:)
   
    CALL bcast_omp_igen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_i2


  SUBROUTINE bcast_omp_i3(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:)

    CALL bcast_omp_igen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_i3


  SUBROUTINE bcast_omp_i4(var)
  IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: Var(:,:,:,:)
   
    CALL bcast_omp_igen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_i4


!! -- Les reels -- !!

  SUBROUTINE bcast_omp_r(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var

    CALL bcast_omp_rgen(Var,1,omp_buffer)

  END SUBROUTINE bcast_omp_r


  SUBROUTINE bcast_omp_r1(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:)
   
    CALL bcast_omp_rgen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_r1


  SUBROUTINE bcast_omp_r2(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:)
   
    CALL bcast_omp_rgen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_r2


  SUBROUTINE bcast_omp_r3(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:)

    CALL bcast_omp_igen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_r3


  SUBROUTINE bcast_omp_r4(var)
  IMPLICIT NONE
    REAL,INTENT(INOUT) :: Var(:,:,:,:)
   
    CALL bcast_omp_rgen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_r4

  
!! -- Les booleans -- !!

  SUBROUTINE bcast_omp_l(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var

    CALL bcast_omp_lgen(Var,1,omp_buffer)

  END SUBROUTINE bcast_omp_l


  SUBROUTINE bcast_omp_l1(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:)
   
    CALL bcast_omp_lgen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_l1


  SUBROUTINE bcast_omp_l2(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:)
   
    CALL bcast_omp_lgen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_l2


  SUBROUTINE bcast_omp_l3(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:)

    CALL bcast_omp_lgen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_l3


  SUBROUTINE bcast_omp_l4(var)
  IMPLICIT NONE
    LOGICAL,INTENT(INOUT) :: Var(:,:,:,:)
   
    CALL bcast_omp_lgen(Var,size(Var),omp_buffer)

  END SUBROUTINE bcast_omp_l4



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Definition des Scatter   --> 4D   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE scatter_omp_i(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    INTEGER :: dummy


     IF (is_omp_root) THEN
      CALL scatter_omp_igen(VarIn,Varout,1,omp_buffer)
     ELSE
      CALL scatter_omp_igen(dummy,Varout,1,omp_buffer)
    ENDIF
    
  END SUBROUTINE scatter_omp_i


  SUBROUTINE scatter_omp_i1(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER :: dummy

    IF (is_omp_root) THEN
      CALL scatter_omp_igen(VarIn,Varout,Size(VarOut,2),omp_buffer)
    ELSE
      CALL scatter_omp_igen(dummy,Varout,Size(VarOut,2),omp_buffer)
    ENDIF
    
  END SUBROUTINE scatter_omp_i1
  
  
  SUBROUTINE scatter_omp_i2(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER :: dummy
    
    IF (is_omp_root) THEN
      CALL scatter_omp_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3),omp_buffer)
    ELSE
      CALL scatter_omp_igen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3),omp_buffer)
    ENDIF

  END SUBROUTINE scatter_omp_i2


  SUBROUTINE scatter_omp_i3(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    INTEGER :: dummy
    
    IF (is_omp_root) THEN
      CALL scatter_omp_igen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),omp_buffer)
    ELSE
      CALL scatter_omp_igen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),omp_buffer)
    ENDIF
  
  END SUBROUTINE scatter_omp_i3




  SUBROUTINE scatter_omp_r(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    REAL :: dummy


     IF (is_omp_root) THEN
      CALL scatter_omp_rgen(VarIn,Varout,1,omp_buffer)
     ELSE
      CALL scatter_omp_rgen(dummy,Varout,1,omp_buffer)
    ENDIF
    
  END SUBROUTINE scatter_omp_r


  SUBROUTINE scatter_omp_r1(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL :: dummy

    IF (is_omp_root) THEN
      CALL scatter_omp_rgen(VarIn,Varout,Size(VarOut,2),omp_buffer)
    ELSE
      CALL scatter_omp_rgen(dummy,Varout,Size(VarOut,2),omp_buffer)
    ENDIF
    
  END SUBROUTINE scatter_omp_r1
  
  
  SUBROUTINE scatter_omp_r2(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL :: dummy
    
    IF (is_omp_root) THEN
      CALL scatter_omp_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3),omp_buffer)
    ELSE
      CALL scatter_omp_rgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3),omp_buffer)
    ENDIF

  END SUBROUTINE scatter_omp_r2


  SUBROUTINE scatter_omp_r3(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    REAL :: dummy
    
    IF (is_omp_root) THEN
      CALL scatter_omp_rgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),omp_buffer)
    ELSE
      CALL scatter_omp_rgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),omp_buffer)
    ENDIF
  
  END SUBROUTINE scatter_omp_r3
  


  SUBROUTINE scatter_omp_l(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    LOGICAL :: dummy


     IF (is_omp_root) THEN
      CALL scatter_omp_lgen(VarIn,Varout,1,omp_buffer)
     ELSE
      CALL scatter_omp_lgen(dummy,Varout,1,omp_buffer)
    ENDIF
    
  END SUBROUTINE scatter_omp_l


  SUBROUTINE scatter_omp_l1(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    LOGICAL :: dummy

    IF (is_omp_root) THEN
      CALL scatter_omp_lgen(VarIn,Varout,Size(VarOut,2),omp_buffer)
    ELSE
      CALL scatter_omp_lgen(dummy,Varout,Size(VarOut,2),omp_buffer)
    ENDIF
    
  END SUBROUTINE scatter_omp_l1
  
  
  SUBROUTINE scatter_omp_l2(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    LOGICAL :: dummy
    
    IF (is_omp_root) THEN
      CALL scatter_omp_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3),omp_buffer)
    ELSE
      CALL scatter_omp_lgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3),omp_buffer)
    ENDIF

  END SUBROUTINE scatter_omp_l2


  SUBROUTINE scatter_omp_l3(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    LOGICAL :: dummy
    
    IF (is_omp_root) THEN
      CALL scatter_omp_lgen(VarIn,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),omp_buffer)
    ELSE
      CALL scatter_omp_lgen(dummy,Varout,Size(VarOut,2)*Size(VarOut,3)*Size(VarOut,4),omp_buffer)
    ENDIF
  
  END SUBROUTINE scatter_omp_l3  
  

  SUBROUTINE gather_omp_i(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut

    INTEGER :: dummy


     IF (is_omp_root) THEN
      CALL gather_omp_igen(VarIn,Varout,1,omp_buffer)
     ELSE
      CALL gather_omp_igen(dummy,Varout,1,omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_i


  SUBROUTINE gather_omp_i1(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    INTEGER :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_igen(VarIn,Varout,Size(VarIn,2),omp_buffer)
    ELSE
      CALL gather_omp_igen(VarIn,dummy,Size(VarIn,2),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_i1


  SUBROUTINE gather_omp_i2(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    INTEGER :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_igen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3),omp_buffer)
    ELSE
      CALL gather_omp_igen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_i2
  

  SUBROUTINE gather_omp_i3(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    INTEGER :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_igen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),omp_buffer)
    ELSE
      CALL gather_omp_igen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_i3



  SUBROUTINE gather_omp_r(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut

    REAL :: dummy


     IF (is_omp_root) THEN
      CALL gather_omp_rgen(VarIn,Varout,1,omp_buffer)
     ELSE
      CALL gather_omp_rgen(VarIn,dummy,1,omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_r


  SUBROUTINE gather_omp_r1(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    REAL :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_rgen(VarIn,Varout,Size(VarIn,2),omp_buffer)
    ELSE
      CALL gather_omp_rgen(VarIn,dummy,Size(VarIn,2),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_r1


  SUBROUTINE gather_omp_r2(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    REAL :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_rgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3),omp_buffer)
    ELSE
      CALL gather_omp_rgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_r2
  

  SUBROUTINE gather_omp_r3(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    REAL :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_rgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),omp_buffer)
    ELSE
      CALL gather_omp_rgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_r3


  SUBROUTINE gather_omp_l(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:) :: VarOut

    LOGICAL :: dummy


     IF (is_omp_root) THEN
      CALL gather_omp_lgen(VarIn,Varout,1,omp_buffer)
     ELSE
      CALL gather_omp_lgen(VarIn,dummy,1,omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_l


  SUBROUTINE gather_omp_l1(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    LOGICAL :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_lgen(VarIn,Varout,Size(VarIn,2),omp_buffer)
    ELSE
      CALL gather_omp_lgen(VarIn,dummy,Size(VarIn,2),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_l1


  SUBROUTINE gather_omp_l2(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    LOGICAL :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_lgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3),omp_buffer)
    ELSE
      CALL gather_omp_lgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_l2
  

  SUBROUTINE gather_omp_l3(VarIn, VarOut)
    USE mod_phys_lmdz_omp_data, ONLY : is_omp_root
    IMPLICIT NONE
  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
    
    LOGICAL :: dummy

    IF (is_omp_root) THEN
      CALL gather_omp_lgen(VarIn,Varout,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),omp_buffer)
    ELSE
      CALL gather_omp_lgen(VarIn,dummy,Size(VarIn,2)*Size(VarIn,3)*Size(VarIn,4),omp_buffer)
    ENDIF
    
  END SUBROUTINE gather_omp_l3




  SUBROUTINE reduce_sum_omp_i(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN)  :: VarIn
    INTEGER,INTENT(OUT) :: VarOut
    
    CALL reduce_sum_omp_igen(VarIn,Varout,1,omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_i

  SUBROUTINE reduce_sum_omp_i1(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:) :: VarOut
    
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),omp_buffer)
   
  END SUBROUTINE reduce_sum_omp_i1
  
  
  SUBROUTINE reduce_sum_omp_i2(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_i2


  SUBROUTINE reduce_sum_omp_i3(VarIn, VarOut)
    IMPLICIT NONE
  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_i3


  SUBROUTINE reduce_sum_omp_i4(VarIn, VarOut)
    IMPLICIT NONE

    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
  
    CALL reduce_sum_omp_igen(VarIn,Varout,Size(VarIn),omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_i4


  SUBROUTINE reduce_sum_omp_r(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN)  :: VarIn
    REAL,INTENT(OUT) :: VarOut
    
    CALL reduce_sum_omp_rgen(VarIn,Varout,1,omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_r

  SUBROUTINE reduce_sum_omp_r1(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:) :: VarOut
    
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),omp_buffer)
   
  END SUBROUTINE reduce_sum_omp_r1
  
  
  SUBROUTINE reduce_sum_omp_r2(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:) :: VarOut
    
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_r2


  SUBROUTINE reduce_sum_omp_r3(VarIn, VarOut)
    IMPLICIT NONE
  
    REAL,INTENT(IN),DIMENSION(:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:) :: VarOut
    
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_r3


  SUBROUTINE reduce_sum_omp_r4(VarIn, VarOut)
    IMPLICIT NONE

    REAL,INTENT(IN),DIMENSION(:,:,:,:)  :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:) :: VarOut
  
    CALL reduce_sum_omp_rgen(VarIn,Varout,Size(VarIn),omp_buffer)
  
  END SUBROUTINE reduce_sum_omp_r4



END MODULE mod_phys_lmdz_omp_transfert





SUBROUTINE bcast_omp_cgen(Var,Nb,Buff)
  IMPLICIT NONE
    
  CHARACTER(LEN=*),INTENT(INOUT) :: Var
  CHARACTER(LEN=*),INTENT(INOUT) :: Buff
  INTEGER,INTENT(IN) :: Nb
    
  INTEGER :: i
  
!$OMP MASTER
    Buff=Var
!$OMP END MASTER
!$OMP BARRIER

  DO i=1,Nb
    Var=Buff
  ENDDO
!$OMP BARRIER      
  
END SUBROUTINE bcast_omp_cgen


      
SUBROUTINE bcast_omp_igen(Var,Nb,Buff)
  IMPLICIT NONE
    
  INTEGER,DIMENSION(Nb),INTENT(INOUT) :: Var
  INTEGER,DIMENSION(Nb),INTENT(INOUT) :: Buff
  INTEGER,INTENT(IN) :: Nb

  INTEGER :: i
    
!$OMP MASTER
  DO i=1,Nb
    Buff(i)=Var(i)
  ENDDO
!$OMP END MASTER
!$OMP BARRIER

  DO i=1,Nb
    Var(i)=Buff(i)
  ENDDO
!$OMP BARRIER        

END SUBROUTINE bcast_omp_igen


SUBROUTINE bcast_omp_rgen(Var,Nb,Buff)
  IMPLICIT NONE
    
  REAL,DIMENSION(Nb),INTENT(INOUT) :: Var
  REAL,DIMENSION(Nb),INTENT(INOUT) :: Buff
  INTEGER,INTENT(IN) :: Nb

  INTEGER :: i
    
!$OMP MASTER
  DO i=1,Nb
    Buff(i)=Var(i)
  ENDDO
!$OMP END MASTER
!$OMP BARRIER

  DO i=1,Nb
    Var(i)=Buff(i)
  ENDDO
!$OMP BARRIER        

END SUBROUTINE bcast_omp_rgen

SUBROUTINE bcast_omp_lgen(Var,Nb,Buff)
  IMPLICIT NONE
    
  LOGICAL,DIMENSION(Nb),INTENT(INOUT) :: Var
  LOGICAL,DIMENSION(Nb),INTENT(INOUT) :: Buff
  INTEGER,INTENT(IN) :: Nb

  INTEGER :: i
    
!$OMP MASTER
  DO i=1,Nb
    Buff(i)=Var(i)
  ENDDO
!$OMP END MASTER
!$OMP BARRIER

  DO i=1,Nb
    Var(i)=Buff(i)
  ENDDO
!$OMP BARRIER        

END SUBROUTINE bcast_omp_lgen





SUBROUTINE scatter_omp_igen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  INTEGER,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
  INTEGER,INTENT(OUT),DIMENSION(klon_omp,dimsize) :: VarOut
  INTEGER,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

  INTEGER :: i,ij
    
!$OMP MASTER
  DO i=1,dimsize
    DO ij=1,klon_mpi
      Buff(ij,i)=VarIn(ij,i)
    ENDDO
  ENDDO  
!$OMP END MASTER
!$OMP BARRIER

  DO i=1,dimsize
    DO ij=1,klon_omp
      VarOut(ij,i)=Buff(klon_omp_begin-1+ij,i)
    ENDDO
  ENDDO
!$OMP BARRIER  

END SUBROUTINE scatter_omp_igen


SUBROUTINE scatter_omp_rgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  REAL,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
  REAL,INTENT(OUT),DIMENSION(klon_omp,dimsize) :: VarOut
  REAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

  INTEGER :: i,ij
    
!$OMP MASTER
  DO i=1,dimsize
    DO ij=1,klon_mpi
      Buff(ij,i)=VarIn(ij,i)
    ENDDO
  ENDDO  
!$OMP END MASTER
!$OMP BARRIER

  DO i=1,dimsize
    DO ij=1,klon_omp
      VarOut(ij,i)=Buff(klon_omp_begin-1+ij,i)
    ENDDO
  ENDDO
!$OMP BARRIER  

END SUBROUTINE scatter_omp_rgen


SUBROUTINE scatter_omp_lgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  LOGICAL,INTENT(IN),DIMENSION(klon_mpi,dimsize) :: VarIn
  LOGICAL,INTENT(OUT),DIMENSION(klon_omp,dimsize) :: VarOut
  LOGICAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

  INTEGER :: i,ij
    
!$OMP MASTER
  DO i=1,dimsize
    DO ij=1,klon_mpi
      Buff(ij,i)=VarIn(ij,i)
    ENDDO
  ENDDO  
!$OMP END MASTER
!$OMP BARRIER

  DO i=1,dimsize
    DO ij=1,klon_omp
      VarOut(ij,i)=Buff(klon_omp_begin-1+ij,i)
    ENDDO
  ENDDO
!$OMP BARRIER  

END SUBROUTINE scatter_omp_lgen





SUBROUTINE gather_omp_igen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  INTEGER,INTENT(IN),DIMENSION(klon_omp,dimsize) :: VarIn
  INTEGER,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
  INTEGER,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

  INTEGER :: i,ij
    
  DO i=1,dimsize
    DO ij=1,klon_omp
      Buff(klon_omp_begin-1+ij,i)=VarIn(ij,i)
    ENDDO
  ENDDO
!$OMP BARRIER  


!$OMP MASTER
  DO i=1,dimsize
    DO ij=1,klon_mpi
      VarOut(ij,i)=Buff(ij,i)
    ENDDO
  ENDDO  
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE gather_omp_igen


SUBROUTINE gather_omp_rgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  REAL,INTENT(IN),DIMENSION(klon_omp,dimsize) :: VarIn
  REAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
  REAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

  INTEGER :: i,ij
    
  DO i=1,dimsize
    DO ij=1,klon_omp
      Buff(klon_omp_begin-1+ij,i)=VarIn(ij,i)
    ENDDO
  ENDDO
!$OMP BARRIER  


!$OMP MASTER
  DO i=1,dimsize
    DO ij=1,klon_mpi
      VarOut(ij,i)=Buff(ij,i)
    ENDDO
  ENDDO  
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE gather_omp_rgen


SUBROUTINE gather_omp_lgen(VarIn,VarOut,dimsize,Buff)
  USE mod_phys_lmdz_omp_data
  USE mod_phys_lmdz_mpi_data, ONLY : klon_mpi 
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  LOGICAL,INTENT(IN),DIMENSION(klon_omp,dimsize) :: VarIn
  LOGICAL,INTENT(OUT),DIMENSION(klon_mpi,dimsize) :: VarOut
  LOGICAL,INTENT(INOUT),DIMENSION(klon_mpi,dimsize) :: Buff

  INTEGER :: i,ij
    
  DO i=1,dimsize
    DO ij=1,klon_omp
      Buff(klon_omp_begin-1+ij,i)=VarIn(ij,i)
    ENDDO
  ENDDO
!$OMP BARRIER  


!$OMP MASTER
  DO i=1,dimsize
    DO ij=1,klon_mpi
      VarOut(ij,i)=Buff(ij,i)
    ENDDO
  ENDDO  
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE gather_omp_lgen


SUBROUTINE reduce_sum_omp_igen(VarIn,VarOut,dimsize,Buff)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  INTEGER,INTENT(IN),DIMENSION(dimsize) :: VarIn
  INTEGER,INTENT(OUT),DIMENSION(dimsize) :: VarOut
  INTEGER,INTENT(INOUT),DIMENSION(dimsize) :: Buff

  INTEGER :: i

!$OMP MASTER
  Buff(:)=0
!$OMP END MASTER
!$OMP BARRIER

!$OMP CRITICAL     
  DO i=1,dimsize
    Buff(i)=Buff(i)+VarIn(i)
  ENDDO
!$OMP END CRITICAL
!$OMP BARRIER  

!$OMP MASTER
  DO i=1,dimsize
    VarOut(i)=Buff(i)
  ENDDO
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE reduce_sum_omp_igen

SUBROUTINE reduce_sum_omp_rgen(VarIn,VarOut,dimsize,Buff)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: dimsize
  REAL,INTENT(IN),DIMENSION(dimsize) :: VarIn
  REAL,INTENT(OUT),DIMENSION(dimsize) :: VarOut
  REAL,INTENT(INOUT),DIMENSION(dimsize) :: Buff

  INTEGER :: i

!$OMP MASTER
  Buff(:)=0
!$OMP END MASTER
!$OMP BARRIER

!$OMP CRITICAL     
  DO i=1,dimsize
    Buff(i)=Buff(i)+VarIn(i)
  ENDDO
!$OMP END CRITICAL
!$OMP BARRIER  

!$OMP MASTER
  DO i=1,dimsize
    VarOut(i)=Buff(i)
  ENDDO
!$OMP END MASTER
!$OMP BARRIER

END SUBROUTINE reduce_sum_omp_rgen
