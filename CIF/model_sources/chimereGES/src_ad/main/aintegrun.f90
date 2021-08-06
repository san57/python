subroutine aintegrun
  
  use chimere_consts
  use chimere_common
  use message_defs
  use master_message_subs
  
  implicit none
  
  !****************************************************************************
  integer :: nh,np,ir
  integer :: ksens
  integer :: ierr
  integer :: iprint
  ! IP
  real(kind=8), allocatable, dimension(:,:):: tabobs
  integer nb,i,ns,ip,ims,ime,izs,izo,nbw,nl
  character(len=16) :: obsname
  ! fin IP
  !****************************************************************************
  ! pour bien commencer l'adjoint
  ihourrun=ihourrun-1
  print*,'IHOURRUN ADJ',ihourrun
  
  call amaster_send_once ! initial aemis
  
  ! IP
  !  fichier sortie valeurs modele
  open(52,file='mod2.txt',form='formatted',action='write')
  call mpi_barrier(mpi_comm_world,ierr)
  ! fin IP
  
  ! send initial aconc and aconco
  call amaster_send_aconc
  
  call mpi_barrier(mpi_comm_world,ierr)
  
  do nh=nhourrun,1,-1
    
    call mpi_barrier(mpi_comm_world,ierr)
    
    ksens = 1
    call areadhour(ksens)
    call checkcfl
    thour = 0.5d0/nphour
    do np=1,nphour
      djul = djul - dun/(nphour*nhourpday)
    enddo
    call master_send_hourly
    call mpi_barrier(mpi_comm_world,ierr)
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print*,' DIRECT INTEGRATION master',ihourrun
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do np=1,nphour
      
      call mpi_barrier(mpi_comm_world,ierr)
      
      call master_locvalues
      
      call mpi_barrier(mpi_comm_world,ierr)
      call master_send_frac_hourly
      
      call mpi_barrier(mpi_comm_world,ierr)
      
      !  Timing update
      thour = thour + dun/nphour
      djul = djul + dun/(nphour*nhourpday)
    
    enddo ! np=1,nphour
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print*,' Adjoint integration master'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do np=nphour,1,-1
      
      !  Timing update
      thour = thour - dun/nphour ! catherine boucle a l envers  
      djul = djul - dun/(nphour*nhourpday)
      
      call mpi_barrier(mpi_comm_world,ierr)
      
      call master_locvalues
      
      call mpi_barrier(mpi_comm_world,ierr)
      call master_send_frac_hourly
      
      call mpi_barrier(mpi_comm_world,ierr)
      
      ! ici dans adjinteg, aoutprint = aconc pour les familles + outprint
      ! mais necessite passages avec workers donc pas fait pour le moment
      call mpi_barrier(mpi_comm_world,ierr)
      call amaster_recv_aconc_bounds
      call mpi_barrier(mpi_comm_world,ierr)
      call amaster_alocvalues
    enddo ! np=nphour,1,-1
    
    
    if(nh.gt.1)ihourrun=ihourrun-1 ! sinon, on n'avance pas!
    
    call arenewhour
  
  enddo ! nh=nhourrun,1,-1  
  
  call mpi_barrier(mpi_comm_world,ierr)
  
  ! Save of all aconcentrations
  call amaster_recv_aconc
  call awrite_aconcs
  
  ! adjoint for iniconc: here to simplify transfers
  call ainiconc
  
  call mpi_barrier(mpi_comm_world,ierr)
  
  ! Receiving and merging tabobs
  do ip = 1, ndoms
    allocate ( tabobs(nobs(ip),13) )
    call mpi_recv(tabobs,nobs(ip)*13,mpi_double_precision,ip,  &
            iar_tabobs,mpi_comm_world,mpi_status_ignore,ierr)
    do i = 1, nobs(ip)
      nb = tabobs(i, 8)
      tabobs_glo(nb,1:6) = tabobs(i,1:6)
      tabobs_glo(nb,8:11) = tabobs(i,9:12)
      tabobs_glo(nb,12) = tabobs(i,7)
    enddo
    deallocate(tabobs)
  enddo
  
  ! writing forward for check
  do i=1,nobs_glo
    write(52,1001)int(tabobs_glo(i,1:5)),species(int(tabobs_glo(i,6)))%name,tabobs_glo(i,7:12)
  enddo
  close(52)
  
  1001 format(5(i,' '),a16,' ',6(e64.55E3))
  
  ! receive and merge aemis
  call amaster_recv_once
  
  call awrite_adj ! from aadjinteg, awrtadjnhours


end subroutine aintegrun




  
