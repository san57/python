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

  module parallel
  USE mod_const_mpi
   
    integer, save :: mpi_size
    integer, save :: mpi_rank
    integer, save :: jj_begin
    integer, save :: jj_end
    integer, save :: jj_nb
    integer, save :: ij_begin
    integer, save :: ij_end
    logical, save :: north_pole
    logical, save :: south_pole
    
    integer, allocatable, save, dimension(:) :: jj_begin_para
    integer, allocatable, save, dimension(:) :: jj_end_para
    integer, allocatable, save, dimension(:) :: jj_nb_para
    integer, save :: OMP_CHUNK
    
 contains
 
    subroutine init_parallel
    USE vampir
#ifdef CPP_COUPLE
    USE mod_prism_proto
#endif
    implicit none
    
      integer :: ierr
      integer :: i,j
      
      
      include 'mpif.h'
      include "dimensions.h"
      include "paramet.h"

      call InitVampir
      call MPI_COMM_SIZE(COMM_LMDZ,mpi_size,ierr)
      call MPI_COMM_RANK(COMM_LMDZ,mpi_rank,ierr)
  
      
      allocate(jj_begin_para(0:mpi_size-1))
      allocate(jj_end_para(0:mpi_size-1))
      allocate(jj_nb_para(0:mpi_size-1))
      
      do i=0,mpi_size-1
        jj_nb_para(i)=(jjm+1)/mpi_size
        if ( i < MOD((jjm+1),mpi_size) ) jj_nb_para(i)=jj_nb_para(i)+1
        
        if (jj_nb_para(i) <= 2 ) then
          
          print *,"Arret : le nombre de bande de lattitude par process est trop faible (<2)."
          print *," ---> diminuez le nombre de CPU ou augmentez la taille en lattitude"
          
          call MPI_ABORT(COMM_LMDZ,-1, ierr)
          
        endif
        
      enddo
      
!      jj_nb_para(0)=11
!      jj_nb_para(1)=25
!      jj_nb_para(2)=25
!      jj_nb_para(3)=12

      j=1
      
      do i=0,mpi_size-1
        
        jj_begin_para(i)=j
        jj_end_para(i)=j+jj_Nb_para(i)-1
        j=j+jj_Nb_para(i)
      
      enddo
      
      jj_begin = jj_begin_para(mpi_rank)
      jj_end   = jj_end_para(mpi_rank)
      jj_nb    = jj_nb_para(mpi_rank)
      
      ij_begin=(jj_begin-1)*iip1+1
      ij_end=jj_end*iip1
      
      if (mpi_rank.eq.0) then
        north_pole=.TRUE.
      else
        north_pole=.FALSE.
      endif
      
      if (mpi_rank.eq.mpi_size-1) then
        south_pole=.TRUE.
      else
        south_pole=.FALSE.
      endif
      
      print *,"jj_begin",jj_begin
      print *,"jj_end",jj_end
      print *,"ij_begin",ij_begin
      print *,"ij_end",ij_end
    
    
    end subroutine init_parallel

    
    subroutine SetDistrib(jj_Nb_New)
    implicit none

    include "dimensions.h"
    include "paramet.h"

      INTEGER,dimension(0:MPI_Size-1) :: jj_Nb_New
      INTEGER :: i
  
      jj_Nb_Para=jj_Nb_New
      
      jj_begin_para(0)=1
      jj_end_para(0)=jj_Nb_Para(0)
      
      do i=1,mpi_size-1
        
        jj_begin_para(i)=jj_end_para(i-1)+1
        jj_end_para(i)=jj_begin_para(i)+jj_Nb_para(i)-1
      
      enddo
      
      jj_begin = jj_begin_para(mpi_rank)
      jj_end   = jj_end_para(mpi_rank)
      jj_nb    = jj_nb_para(mpi_rank)
      
      ij_begin=(jj_begin-1)*iip1+1
      ij_end=jj_end*iip1

    end subroutine SetDistrib



    
    subroutine Finalize_parallel
#ifdef CPP_COUPLE
    use mod_prism_proto
#endif
    implicit none

    include "dimensions.h"
    include "paramet.h"
    include 'mpif.h'
      integer :: ierr
      
      deallocate(jj_begin_para)
      deallocate(jj_end_para)
      deallocate(jj_nb_para)

#ifdef CPP_COUPLE
     call prism_terminate_proto(ierr)
     IF (ierr .ne. PRISM_Ok) THEN
       WRITE(*,*) ' Problem in prism_terminate_proto '
       STOP
     endif
#else
      call MPI_FINALIZE(ierr)
#endif
    
    end subroutine Finalize_parallel
    
    subroutine Pack_Data(Field,ij,ll,row,Buffer)
    implicit none

    include "dimensions.h"
    include "paramet.h"

      integer, intent(in) :: ij,ll,row
      real,dimension(ij,ll),intent(in) ::Field
      real,dimension(ll*iip1*row), intent(out) :: Buffer
      
      integer :: Pos
      integer :: i,l
      
      Pos=0
      do l=1,ll
        do i=1,row*iip1
          Pos=Pos+1
          Buffer(Pos)=Field(i,l)
        enddo
      enddo
      
    end subroutine Pack_data
    
    subroutine Unpack_Data(Field,ij,ll,row,Buffer)
    implicit none

    include "dimensions.h"
    include "paramet.h"

      integer, intent(in) :: ij,ll,row
      real,dimension(ij,ll),intent(out) ::Field
      real,dimension(ll*iip1*row), intent(in) :: Buffer
      
      integer :: Pos
      integer :: i,l
      
      Pos=0
      
      do l=1,ll
        do i=1,row*iip1
          Pos=Pos+1
          Field(i,l)=Buffer(Pos)
        enddo
      enddo
      
    end subroutine UnPack_data
    
    subroutine exchange_hallo(Field,ij,ll,up,down)
    USE Vampir
    implicit none
    include "dimensions.h"
    include "paramet.h"
    include 'mpif.h'
    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: up,down
      
      INTEGER :: ierr
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
      INTEGER, DIMENSION(4) :: Request
      INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: Status
      INTEGER :: NbRequest
      REAL, dimension(:),allocatable :: Buffer_Send_up,Buffer_Send_down
      REAL, dimension(:),allocatable :: Buffer_Recv_up,Buffer_Recv_down
      INTEGER :: Buffer_size
      
      call MPI_Barrier(COMM_LMDZ,ierr)
      call VTb(VThallo)
      
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
      
      if (up.eq.0) then
        SendDown=.FALSE.
        RecvUp=.FALSE.
      endif
      
      if (down.eq.0) then
        SendUp=.FALSE.
        RecvDown=.FALSE.
      endif
      
      NbRequest=0
  
      IF (SendUp) THEN
        NbRequest=NbRequest+1
        buffer_size=down*iip1*ll
        allocate(Buffer_Send_up(Buffer_size))
        call PACK_Data(Field(ij_begin,1),ij,ll,down,Buffer_Send_up)
        call MPI_ISSEND(Buffer_send_up,Buffer_Size,MPI_REAL8,MPI_Rank-1,1,     &
                        COMM_LMDZ,Request(NbRequest),ierr)
      ENDIF
  
      IF (SendDown) THEN
        NbRequest=NbRequest+1
        
        buffer_size=up*iip1*ll
        allocate(Buffer_Send_down(Buffer_size))
        call PACK_Data(Field(ij_end+1-up*iip1,1),ij,ll,up,Buffer_send_down)
        
        call MPI_ISSEND(Buffer_send_down,Buffer_Size,MPI_REAL8,MPI_Rank+1,1,     &
                        COMM_LMDZ,Request(NbRequest),ierr)
      ENDIF
    
  
      IF (RecvUp) THEN
        NbRequest=NbRequest+1
        buffer_size=up*iip1*ll
        allocate(Buffer_recv_up(Buffer_size))
        
        call MPI_IRECV(Buffer_recv_up,Buffer_size,MPI_REAL8,MPI_Rank-1,1,  &
                        COMM_LMDZ,Request(NbRequest),ierr)
     
       
      ENDIF
  
      IF (RecvDown) THEN
        NbRequest=NbRequest+1
        buffer_size=down*iip1*ll
        allocate(Buffer_recv_down(Buffer_size))
        
        call MPI_IRECV(Buffer_recv_down,Buffer_size,MPI_REAL8,MPI_Rank+1,1,     &
                        COMM_LMDZ,Request(NbRequest),ierr)
      
        
      ENDIF
  
      if (NbRequest > 0) call MPI_WAITALL(NbRequest,Request,Status,ierr)
      IF (RecvUp)  call Unpack_Data(Field(ij_begin-up*iip1,1),ij,ll,up,Buffer_Recv_up)
      IF (RecvDown) call Unpack_Data(Field(ij_end+1,1),ij,ll,down,Buffer_Recv_down)

      call VTe(VThallo)
      call MPI_Barrier(COMM_LMDZ,ierr)
      RETURN
      
    end subroutine exchange_Hallo
    
    
    subroutine Gather_Field(Field,ij,ll,rank)
      implicit none
      include "dimensions.h"
      include "paramet.h"
      include 'mpif.h'
      
      INTEGER :: ij,ll,rank
      REAL, dimension(ij,ll) :: Field
      REAL, dimension(:),allocatable :: Buffer_send
      REAL, dimension(:),allocatable :: Buffer_Recv
      INTEGER, dimension(0:MPI_Size-1) :: Recv_count, displ
      INTEGER :: ierr
      INTEGER ::i
      
      if (ij==ip1jmp1) then
        allocate(Buffer_send(iip1*ll*(jj_end-jj_begin+1)))
        call Pack_Data(Field(ij_begin,1),ij,ll,jj_end-jj_begin+1,Buffer_send)
      else if (ij==ip1jm) then
        allocate(Buffer_send(iip1*ll*(min(jj_end,jjm)-jj_begin+1)))
        call Pack_Data(Field(ij_begin,1),ij,ll,min(jj_end,jjm)-jj_begin+1,Buffer_send)
      else
        print *,ij
        stop 'erreur dans Gather_Field'
      endif
      
      ALLOCATE(Buffer_Recv(ij*ll))
      IF (MPI_Rank==rank) THEN
        do i=0,MPI_Size-1
          
          if (ij==ip1jmp1) then
            Recv_count(i)=(jj_end_para(i)-jj_begin_para(i)+1)*ll*iip1
          else if (ij==ip1jm) then
            Recv_count(i)=(min(jj_end_para(i),jjm)-jj_begin_para(i)+1)*ll*iip1
          else
            stop 'erreur dans Gather_Field'
          endif
          
          if (i==0) then
            displ(i)=0
          else
            displ(i)=displ(i-1)+Recv_count(i-1)
          endif
          
        enddo
        
      endif
      
      call MPI_GATHERV(Buffer_send,(min(ij_end,ij)-ij_begin+1)*ll,MPI_REAL8,   &
                        Buffer_Recv,Recv_count,displ,MPI_REAL8,rank,COMM_LMDZ,ierr)
      
      if (MPI_Rank==rank) then
      
        if (ij==ip1jmp1) then
          do i=0,MPI_Size-1
            call Unpack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                 &
                             jj_end_para(i)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
          enddo
        else if (ij==ip1jm) then
          do i=0,MPI_Size-1
             call Unpack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                       &
                             min(jj_end_para(i),jjm)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
          enddo
        endif
      
      endif
      
    end subroutine Gather_Field

  subroutine Scatter_Field(Field,ij,ll,rank)
    ! distribue le champ Field de dimension (ij, ll) 
    ! du processus "rank" aux processus autres que rank
    ! en suivant le decoupage utilise a ce moment
    implicit none
    include "dimensions.h"
    include "paramet.h"
    include 'mpif.h'

    INTEGER :: ij,ll,rank
    REAL, dimension(ij,ll) :: Field

    ! Les tableaux Buffer_send et Buffer_Recv n'ont pas le role indique par leur nom

    REAL, dimension(iip1*ll*(jjm-jj_begin+1)) :: Buffer_send
    REAL, dimension(ij*ll) :: Buffer_Recv
    INTEGER, dimension(0:MPI_Size-1) :: Recv_count, displ
    INTEGER :: ierr
    INTEGER ::i


    ! on calcule dans "displ" l'offset dans le tableau "Buffer_Recv" des tranches 
    ! reordonnees et devant etre envoyees au processus "i"
    if (MPI_Rank==rank) then
       do i=0,MPI_Size-1

          if (ij==ip1jmp1) then
             Recv_count(i)=(jj_end_para(i)-jj_begin_para(i)+1)*ll*iip1
          else if (ij==ip1jm) then
             Recv_count(i)=(min(jj_end_para(i),jjm)-jj_begin_para(i)+1)*ll*iip1
          else
             stop 'erreur dans Scatter_Field'
          endif

          if (i==0) then
             displ(i)=0
          else
             displ(i)=displ(i-1)+Recv_count(i-1)
          endif

       enddo

    endif

    ! on reordonne les donnees du tableau Field dans le tableau Buffer_Recv avant expedition
    if (MPI_Rank==rank) then

       if (ij==ip1jmp1) then
          do i=0,MPI_Size-1
             call Pack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                 &
                  jj_end_para(i)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
          enddo
       else if (ij==ip1jm) then
          do i=0,MPI_Size-1
             call Pack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                       &
                  min(jj_end_para(i),jjm)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
          enddo
       endif

    endif

    call MPI_SCATTERV(Buffer_Recv, Recv_count, displ, MPI_REAL8,   &
         Buffer_send, (min(ij_end,ij)-ij_begin+1)*ll, MPI_REAL8, rank, COMM_LMDZ, ierr)

    if (ij==ip1jmp1) then
       call Unpack_Data(Field(ij_begin,1),ij,ll,jj_end-jj_begin+1,Buffer_send)
    else if (ij==ip1jm) then
       call Unpack_Data(Field(ij_begin,1),ij,ll,min(jj_end,jjm)-jj_begin+1,Buffer_send)
    else
       print *,ij
       stop 'erreur dans Scatter_Field'
    endif

  end subroutine Scatter_Field
 
    
    subroutine AllGather_Field(Field,ij,ll)
    implicit none
    include "dimensions.h"
    include "paramet.h"
    include 'mpif.h'
    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: ierr
      
      call Gather_Field(Field,ij,ll,0)
      call MPI_BCAST(Field,ij*ll,MPI_REAL8,0,COMM_LMDZ,ierr)
      
    end subroutine AllGather_Field
    
   subroutine Broadcast_Field(Field,ij,ll,rank)
    implicit none
    include "dimensions.h"
    include "paramet.h"
    include 'mpif.h'
    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: rank
      INTEGER :: ierr
      
      call MPI_BCAST(Field,ij*ll,MPI_REAL8,rank,COMM_LMDZ,ierr)
      
    end subroutine Broadcast_Field
    
   
    /*
  Subroutine verif_hallo(Field,ij,ll,up,down)
    implicit none
    include "dimensions.h"
    include "paramet.h"
    include 'mpif.h'
    
      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: up,down
      
      REAL,dimension(ij,ll): NewField
      
      NewField=0
      
      ijb=ij_begin
      ije=ij_end
      if (north_pole)
      NewField(ij_be
*/
  end module parallel
