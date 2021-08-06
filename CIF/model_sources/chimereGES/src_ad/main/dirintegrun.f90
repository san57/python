subroutine dirintegrun
  
  !  Integration of the model                                             
  
  use chimere_consts
  use chimere_common
  use message_defs
  use master_message_subs
  
  implicit none
  
  !include 'mpif.h'
  
  
  !****************************************************************************
  integer :: nh,np,ir
  integer :: ksens
  integer :: ierr
  integer :: iprint
  ! IP
  real(kind=8), allocatable, dimension(:,:):: tabobs
  integer nb,i,ns,ip,ims,ime,izs,ize,nbw
  character(len=16) :: obsname
  ! end IP
  !****************************************************************************
  
  iprint = 0
  ihourrun = 0
  thour = 0.5d0/nphour
  ksens = 2
  call mpi_barrier(mpi_comm_world,ierr)
  
  !print*,' receive initial locvalues from workers for printout'
  call master_recv_locvalues
  call mpi_barrier(mpi_comm_world,ierr)
  
  ! printout of the initial concentrations and meteo parameters
  call outprint(iprint,.true.,1,NetCDF_output, NetCDF_parout)
  call mpi_barrier(mpi_comm_world,ierr)
  
  ! Loop on physical time steps
  do nh=1,nhourrun
    
    call mpi_barrier(mpi_comm_world,ierr)
    ksens = 2
    call readhour(ksens)
    call checkcfl
    thour = 0.5d0/nphour
    call master_send_hourly
    call mpi_barrier(mpi_comm_world,ierr)
    
    do np=1,nphour
      
      call mpi_barrier(mpi_comm_world,ierr)
      
      call master_locvalues ! also sends conc boundaries
      
      call mpi_barrier(mpi_comm_world,ierr)
      call master_send_frac_hourly
      call mpi_barrier(mpi_comm_world,ierr)
      
      !  Timing update
      thour = thour + dun/nphour
      djul = djul + dun/(nphour*nhourpday)
      if(np.eq.nphour) ihourrun = ihourrun + 1
      
      ! Printout. I need locvalues and conc
      if(np.eq.nphour) then
        call master_recv_locvalues
        iprint = iprint+1
        call outprint(iprint,.false.,nh,NetCDF_output, NetCDF_parout) ! also receives conc from workers
      endif
    enddo ! np=1,nphour
    
    call renewhour
    
    !print*,'NOT NEEDED IN ADJ: printing of depos',nh
  
  enddo ! nh=1,nhourrun!
  
  !IP receiving and merging obs
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
  
  ! writing pour verification de la partie directe
  do i=1,nobs_glo
    write(2,1001)int(tabobs_glo(i,1:5)),species(int(tabobs_glo(i,6)))%name,tabobs_glo(i,7:12)
  enddo
  close(2)
  1001 format(5(i,' '),a16,' ',6(e64.55E3))
  ! fin IP
  
  if(ad.eq.0)deallocate(tabobs_glo)
  ! end IP

end subroutine dirintegrun
