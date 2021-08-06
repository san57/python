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

module mod_Hallo
USE parallel
implicit none

  integer, parameter :: MaxRequest=200
  integer, parameter :: MaxProc=80
  integer, parameter :: MaxBufferSize=1024*1024*16
  integer, parameter :: ListSize=1000
  
  integer,save       :: MaxBufferSize_Used
 
    real,save,pointer,dimension(:) :: Buffer

   integer,dimension(Listsize) :: Buffer_Pos
   integer :: Index_Pos
   
  type Hallo
    real, dimension(:,:),pointer :: Field
    integer :: offset
    integer :: size
    integer :: NbLevel
    integer :: Stride
  end type Hallo
  
  type request_SR
    integer :: NbRequest=0
    integer :: Pos
    integer :: Index 
    type(Hallo),dimension(MaxRequest) :: Hallo
    integer :: MSG_Request
  end type request_SR

  type request
    type(request_SR),dimension(0:MaxProc-1) :: RequestSend
    type(request_SR),dimension(0:MaxProc-1) :: RequestRecv
    integer :: tag=1
  end type request
  
    
  contains

  subroutine Init_mod_hallo
    implicit none

    Index_Pos=1
    Buffer_Pos(Index_Pos)=1
    MaxBufferSize_Used=0

    CALL create_global_mpi_buffer
    
  end subroutine init_mod_hallo


  SUBROUTINE create_standard_mpi_buffer
  IMPLICIT NONE
    
    ALLOCATE(Buffer(MaxBufferSize))
    
  END SUBROUTINE create_standard_mpi_buffer
  

  SUBROUTINE create_global_mpi_buffer
  IMPLICIT NONE
  INCLUDE 'mpif.h'  
    POINTER (Pbuffer,MPI_Buffer(MaxBufferSize))
    REAL :: MPI_Buffer
    INTEGER(KIND=MPI_ADDRESS_KIND) :: BS 
    INTEGER :: i,ierr


      Bs=8*MaxBufferSize
      CALL MPI_ALLOC_MEM(BS,MPI_INFO_NULL,Pbuffer,ierr)
      DO i=1,MaxBufferSize
	MPI_Buffer(i)=i
      ENDDO
     
      CALL  Associate_buffer(MPI_Buffer)
      
  CONTAINS
     
     SUBROUTINE Associate_buffer(MPI_Buffer)
     IMPLICIT NONE
       REAL,DIMENSION(:),target :: MPI_Buffer  

         Buffer=>MPI_Buffer
 
      END SUBROUTINE  Associate_buffer
                                      
  END SUBROUTINE create_global_mpi_buffer


      
  subroutine allocate_buffer(Size,Index,Pos)
  implicit none
    integer :: Size
    integer :: Index
    integer :: Pos

    if (Buffer_pos(Index_pos)+Size>MaxBufferSize_Used) MaxBufferSize_Used=Buffer_pos(Index_pos)+Size  
    if (Buffer_pos(Index_pos)+Size>MaxBufferSize) then
      print *,'STOP :: La taille de MaxBufferSize dans mod_hallo.F90 est trop petite !!!!'
      stop
    endif
    
    if (Index_pos>=ListSize) then
      print *,'STOP :: La taille de ListSize dans mod_hallo.F90 est trop petite !!!!'
      stop
    endif
     
    Pos=Buffer_Pos(Index_Pos)
    Buffer_Pos(Index_pos+1)=Buffer_Pos(Index_Pos)+Size
    Index_Pos=Index_Pos+1
    Index=Index_Pos
    
  end subroutine allocate_buffer
     
  subroutine deallocate_buffer(Index)
  implicit none
    integer :: Index
    
    Buffer_Pos(Index)=-1
    
    do while (Buffer_Pos(Index_Pos)==-1 .and. Index_Pos>1)
      Index_Pos=Index_Pos-1
    end do

  end subroutine deallocate_buffer  
  
  subroutine SetTag(a_request,tag)
  implicit none
    type(request):: a_request
    integer :: tag
    
    a_request%tag=tag
  end subroutine SetTag
  
  
  subroutine Init_Hallo(Field,Stride,NbLevel,offset,size,NewHallo)
    integer :: Stride
    integer :: NbLevel
    integer :: size
    integer :: offset
    real, dimension(Stride,NbLevel),target :: Field
    type(Hallo) :: NewHallo
    
    NewHallo%Field=>Field
    NewHallo%Stride=Stride
    NewHallo%NbLevel=NbLevel
    NewHallo%size=size
    NewHallo%offset=offset
    
    
  end subroutine Init_Hallo
  
  subroutine Register_SendField(Field,ij,ll,offset,size,target,a_request)
  implicit none

    include "dimensions.h"
    include "paramet.h"    
    include 'mpif.h'
    
      INTEGER :: ij,ll,offset,size,target
      REAL, dimension(ij,ll) :: Field
      type(request),target :: a_request
      type(request_SR),pointer :: Ptr_request

      Ptr_Request=>a_request%RequestSend(target)
      Ptr_Request%NbRequest=Ptr_Request%NbRequest+1
      if (Ptr_Request%NbRequest>=MaxRequest) then
        print *,'STOP :: La taille de MaxRequest dans mod_hallo.F90 est trop petite !!!!'
        stop
      endif      
      call Init_Hallo(Field,ij,ll,offset,size,Ptr_request%Hallo(Ptr_Request%NbRequest))
      
   end subroutine Register_SendField      
      
  subroutine Register_RecvField(Field,ij,ll,offset,size,target,a_request)
  implicit none

 include "dimensions.h"
 include "paramet.h"    
 include 'mpif.h'
    
      INTEGER :: ij,ll,offset,size,target
      REAL, dimension(ij,ll) :: Field
      type(request),target :: a_request
      type(request_SR),pointer :: Ptr_request

      Ptr_Request=>a_request%RequestRecv(target)
      Ptr_Request%NbRequest=Ptr_Request%NbRequest+1
      
      if (Ptr_Request%NbRequest>=MaxRequest) then
        print *,'STOP :: La taille de MaxRequest dans mod_hallo.F90 est trop petite !!!!'
        stop
      endif   
            
      call Init_Hallo(Field,ij,ll,offset,size,Ptr_request%Hallo(Ptr_Request%NbRequest))

      
   end subroutine Register_RecvField      
  
  subroutine Register_SwapField(FieldS,FieldR,ij,ll,jj_Nb_New,a_request)
  
      implicit none
 include "dimensions.h"
 include "paramet.h"    
 include 'mpif.h'
    
    INTEGER :: ij,ll
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    type(request) :: a_request
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    do i=0,MPI_Size-1
      if (i /= MPI_Rank) then
        jjb=max(jj_begin_new(i),jj_begin)
        jje=min(jj_end_new(i),jj_end)
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_SendField(FieldS,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
        jjb=max(jj_begin_new(MPI_Rank),jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),jj_end_Para(i))
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_RecvField(FieldR,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
      endif
    enddo
    
  end subroutine Register_SwapField    
  
  
    subroutine Register_SwapFieldHallo(FieldS,FieldR,ij,ll,jj_Nb_New,Up,Down,a_request)
  
      implicit none
 include "dimensions.h"
 include "paramet.h"    
 include 'mpif.h'
    
    INTEGER :: ij,ll,Up,Down
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    type(request) :: a_request
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    do i=0,MPI_Size-1
      jj_begin_New(i)=max(1,jj_begin_New(i)-Up)
      jj_end_New(i)=min(jjp1,jj_end_new(i)+Down)
    enddo
   
    do i=0,MPI_Size-1
      if (i /= MPI_Rank) then
        jjb=max(jj_begin_new(i),jj_begin)
        jje=min(jj_end_new(i),jj_end)
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_SendField(FieldS,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
        jjb=max(jj_begin_new(MPI_Rank),jj_begin_Para(i))
        jje=min(jj_end_new(MPI_Rank),jj_end_Para(i))
        
        if (ij==ip1jm .and. jje==jjp1) jje=jjm
        
        if (jje >= jjb) then
          call Register_RecvField(FieldR,ij,ll,jjb,jje-jjb+1,i,a_request) 
        endif
        
      endif
    enddo
    
  end subroutine Register_SwapFieldHallo
  
  subroutine Register_Hallo(Field,ij,ll,RUp,Rdown,SUp,SDown,a_request)
  
      implicit none
 include "dimensions.h"
 include "paramet.h"    
 include 'mpif.h'
    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: Sup,Sdown,rup,rdown
      type(request) :: a_request
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
   
 
      SendUp=.TRUE.
      SendDown=.TRUE.
      RecvUp=.TRUE.
      RecvDown=.TRUE.
        
      IF (north_pole) THEN
        SendUp=.FALSE.
        RecvUp=.FALSE.
      ENDIF
  
      IF (south_pole) THEN
        SendDown=.FALSE.
        RecvDown=.FALSE.
      ENDIF
      
      if (Sup.eq.0) then
        SendUp=.FALSE.
       endif
      
      if (Sdown.eq.0) then
        SendDown=.FALSE.
      endif

      if (Rup.eq.0) then
        RecvUp=.FALSE.
      endif
      
      if (Rdown.eq.0) then
        RecvDown=.FALSE.
      endif
      
      IF (SendUp) THEN
        call Register_SendField(Field,ij,ll,jj_begin,SUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (SendDown) THEN
        call Register_SendField(Field,ij,ll,jj_end-SDown+1,SDown,MPI_Rank+1,a_request)
      ENDIF
    
  
      IF (RecvUp) THEN
        call Register_RecvField(Field,ij,ll,jj_begin-Rup,RUp,MPI_Rank-1,a_request)
      ENDIF
  
      IF (RecvDown) THEN
        call Register_RecvField(Field,ij,ll,jj_end+1,RDown,MPI_Rank+1,a_request)
      ENDIF
  
    end subroutine Register_Hallo
    
    subroutine SendRequest(a_Request)
      implicit none

 include "dimensions.h"
 include "paramet.h"
 include 'mpif.h'

      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer :: SizeBuffer
      integer :: i,rank,l,ij,Pos,ierr
      integer :: offset
!      real,dimension(:),pointer :: Buffer
      real,dimension(:,:),pointer :: Field
      integer :: Nb
       
      do rank=0,MPI_SIZE-1
      
        Req=>a_Request%RequestSend(rank)
        
        SizeBuffer=0
        do i=1,Req%NbRequest
          PtrHallo=>Req%Hallo(i)
          SizeBuffer=SizeBuffer+PtrHallo%size*PtrHallo%NbLevel*iip1
        enddo
      
        if (SizeBuffer>0) then
       
!          allocate(Req%Buffer(SizeBuffer))
          call allocate_buffer(SizeBuffer,Req%Index,Req%pos)

          Pos=Req%Pos
!          Buffer=>req%Buffer
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
            Nb=iip1*PtrHallo%size-1
            Field=>PtrHallo%Field
            
            do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        Buffer(Pos+ij)=Field(Offset+ij,l)
	      enddo
!              Buffer(Pos:Pos+Nb)=Field(offset:offset+Nb,l)
              
              Pos=Pos+Nb+1
            enddo
            
          enddo
    
!         print *, 'process',MPI_RANK,'ISSEND: requette ',a_request%tag,'au process',rank,'de taille',SizeBuffer
!         call MPI_ISSEND(Req%Buffer,SizeBuffer,MPI_REAL8,rank,a_request%tag,     &
!                         COMM_LMDZ,Req%MSG_Request,ierr)
         call MPI_ISSEND(Buffer(req%Pos),SizeBuffer,MPI_REAL8,rank,a_request%tag,     &
                         COMM_LMDZ,Req%MSG_Request,ierr)

        endif

    enddo
   
           
      do rank=0,MPI_SIZE-1
         
          Req=>a_Request%RequestRecv(rank)
          SizeBuffer=0
          
	  do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            SizeBuffer=SizeBuffer+PtrHallo%size*PtrHallo%NbLevel*iip1
          enddo
        
          if (SizeBuffer>0) then
!            allocate(Req%Buffer(SizeBuffer))
             call allocate_buffer(SizeBuffer,Req%Index,Req%Pos)
!            print *, 'process',MPI_RANK,'IRECV: requette ',a_request%tag,'au process',rank,'de taille',SizeBuffer
            
!	    call MPI_IRECV(Req%Buffer,SizeBuffer,MPI_REAL8,rank,a_request%tag,     &
!                           COMM_LMDZ,Req%MSG_Request,ierr)
	    call MPI_IRECV(Buffer(Req%Pos),SizeBuffer,MPI_REAL8,rank,a_request%tag,     &
                           COMM_LMDZ,Req%MSG_Request,ierr)

          endif
      
      enddo
                        
   end subroutine SendRequest 
   
   subroutine WaitRequest(a_Request)
   implicit none
   
 include "dimensions.h"
 include "paramet.h"
 include 'mpif.h'   
      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(2*mpi_size) :: TabRequest
      integer, dimension(MPI_STATUS_SIZE,2*mpi_size) :: TabStatus
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset
      integer :: Nb
      
      
      NbRequest=0
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
     
      if (NbRequest>0) call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
      
      do rank=0,MPI_Size-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
	    Nb=iip1*PtrHallo%size-1
            
	    do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        PtrHallo%Field(offset+ij,l)=Buffer(Pos+ij)
	      enddo
!              PtrHallo%Field(offset:offset+Nb,l)=Buffer(Pos:Pos+Nb)
!	      do ij=offset,offset+iip1*PtrHallo%size-1
!                PtrHallo%Field(ij,l)=Buffer(Pos)
!                Pos=Pos+1
!              enddo
              Pos=Pos+Nb+1
	    enddo
	    
          enddo
        endif
      enddo
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
              
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
     
      a_request%tag=1
    end subroutine WaitRequest
     
   subroutine WaitSendRequest(a_Request)
   implicit none
   
 include "dimensions.h"
 include "paramet.h"
 include 'mpif.h'   
      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      integer, dimension(mpi_size) :: TabRequest
      integer, dimension(MPI_STATUS_SIZE,mpi_size) :: TabStatus
      integer :: NbRequest
      integer :: rank,ierr
      
      NbRequest=0
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
      

      if (NbRequest>0) call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
      
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestSend(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
              
      a_request%tag=1
    end subroutine WaitSendRequest
    
   subroutine WaitRecvRequest(a_Request)
   implicit none
   
 include "dimensions.h"
 include "paramet.h"
 include 'mpif.h'   
      
      type(request),target :: a_request
      type(request_SR),pointer :: Req
      type(Hallo),pointer :: PtrHallo
      integer, dimension(mpi_size) :: TabRequest
      integer, dimension(MPI_STATUS_SIZE,mpi_size) :: TabStatus
      integer :: NbRequest
      integer :: i,rank,pos,ij,l,ierr
      integer :: offset,Nb
      
      
      NbRequest=0
      
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          NbRequest=NbRequest+1
          TabRequest(NbRequest)=Req%MSG_Request
        endif
      enddo
     
      
      if (NbRequest>0) call MPI_WAITALL(NbRequest,TabRequest,TabStatus,ierr)
      
      do rank=0,MPI_Size-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          Pos=Req%Pos
          do i=1,Req%NbRequest
            PtrHallo=>Req%Hallo(i)
            offset=(PtrHallo%offset-1)*iip1+1
	    Nb=iip1*PtrHallo%size-1
            
	    do l=1,PtrHallo%NbLevel
!cdir NODEP
              do ij=0,Nb
	        PtrHallo%Field(offset+ij,l)=Buffer(Pos+ij)
	      enddo
                 Pos=Pos+Nb+1
            enddo
          enddo
        endif
      enddo
      
           
      do rank=0,MPI_SIZE-1
        Req=>a_request%RequestRecv(rank)
        if (Req%NbRequest>0) then
          call deallocate_buffer(Req%Index)
          Req%NbRequest=0 
        endif
      enddo
     
      a_request%tag=1
    end subroutine WaitRecvRequest
    
    
    
    subroutine CopyField(FieldS,FieldR,ij,ll,jj_Nb_New)
  
      implicit none
    include "dimensions.h"
    include "paramet.h"    
    include 'mpif.h'
    
    INTEGER :: ij,ll
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New
    
    integer ::i,jje,jjb,ijb,ije
    
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo
    
    jjb=max(jj_begin,jj_begin_new(MPI_Rank))
    jje=min(jj_end,jj_end_new(MPI_Rank))
    if (ij==ip1jm) jje=min(jje,jjm)

    if (jje >= jjb) then
      ijb=(jjb-1)*iip1+1
      ije=jje*iip1
      FieldR(ijb:ije,1:ll)=FieldS(ijb:ije,1:ll)
    endif

  end subroutine CopyField    

  subroutine CopyFieldHallo(FieldS,FieldR,ij,ll,jj_Nb_New,Up,Down)
  
      implicit none
    include "dimensions.h"
    include "paramet.h"    
    include 'mpif.h'
    
    INTEGER :: ij,ll,Up,Down
    REAL, dimension(ij,ll) :: FieldS
    REAL, dimension(ij,ll) :: FieldR
    integer,dimension(0:MPI_Size-1) :: jj_Nb_New   
    integer,dimension(0:MPI_Size-1) :: jj_Begin_New,jj_End_New

    integer ::i,jje,jjb,ijb,ije

     
    jj_begin_New(0)=1
    jj_End_New(0)=jj_Nb_New(0)
    do i=1,MPI_Size-1
      jj_begin_New(i)=jj_end_New(i-1)+1
      jj_end_New(i)=jj_begin_new(i)+jj_Nb_New(i)-1
    enddo

        
    jjb=max(jj_begin,jj_begin_new(MPI_Rank)-Up)
    jje=min(jj_end,jj_end_new(MPI_Rank)+Down)
    if (ij==ip1jm) jje=min(jje,jjm)
    
    
    if (jje >= jjb) then
      ijb=(jjb-1)*iip1+1
      ije=jje*iip1
      FieldR(ijb:ije,1:ll)=FieldS(ijb:ije,1:ll)
    endif
   end subroutine CopyFieldHallo        
          
end module mod_Hallo 
   
