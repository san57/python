subroutine inidoms
  
  use chimere_consts
  use chimere_common
  use message_defs
  use master_message_subs
  
  implicit none
  
  
  !****************************************************************************
  integer :: ierr
  ! IP
  logical :: file_exists
  real(kind = 8), allocatable, dimension(:, :) :: tabobs
  integer nb, i, ns, ip, ims, ime, izs, ize, nbw
  character(len = 16) :: obsname
  ! fin IP
  !****************************************************************************
  
  ndoms = nzdoms * nmdoms
  
  if(ndoms==0) then
    print *
    print *, ' *** Error : NZDOMS or NMDOMS should not be null ! ***'
    print *, ' *** Check your configuration script. Exiting.     ***'
    call mpi_finalize(ierr)
    call exit1('Exiting')
  end if
  
  call master_init_mpi
  
  call master_send_once
  
  !print*,' IP reading obs.bin'
  ! XX RQ: pas tres efficace en cas de satellite??XX
  inquire(file = 'obs.txt', exist = file_exists)
  if (file_exists) then
    open(1, file = 'obs.txt', form = 'formatted', action = 'read')
    read(1, *) nb, nobs_glo
  else
    nobs_glo = 0
  end if
  
  if(fwd.eq.1)allocate (tabobs_glo(nobs_glo, 11))
  if(fwd.eq.0)allocate (tabobs_glo(nobs_glo, 12))
  tabobs_glo(:, :) = 0.
  
  print*, 'Total number of obs=', nobs_glo, rank
  do i = 1, nobs_glo
    if(ad.eq.0)read(1, *)tabobs_glo(i, 1:5), obsname
    if(ad.eq.1)read(1, *)tabobs_glo(i, 1:5), obsname, tabobs_glo(i, 7)
    do ns = 1, nspec
      if(species(ns)%name==obsname)nb = ns
    enddo
    tabobs_glo(i, 6) = dble(nb)
  enddo
  close(1)
  !
  ! distributing on the tiles
  allocate(nobs(ndoms))
  nobs(:) = 0
  do i = 1, nobs_glo
    do ip = 1, ndoms
      ims = dom(ip)%imstart
      ime = dom(ip)%imend
      izs = dom(ip)%izstart
      ize = dom(ip)%izend
      if(tabobs_glo(i, 3)>=ims .and. tabobs_glo(i, 3)<=ime) then
        if(tabobs_glo(i, 4)>=izs .and. tabobs_glo(i, 4)<=ize) then
          nobs(ip) = nobs(ip) + 1
        endif
      endif
    enddo
  enddo
  do ip = 1, ndoms
    print*, 'tile number', ip, ' number of obs=', nobs(ip)
    if(fwd.eq.1)allocate (tabobs(nobs(ip), 12))
    if(fwd.eq.0)allocate (tabobs(nobs(ip), 13))
    nb = 0
    do i = 1, nobs_glo
      ims = dom(ip)%imstart
      ime = dom(ip)%imend
      izs = dom(ip)%izstart
      ize = dom(ip)%izend
      if(tabobs_glo(i, 3)>=ims .and. tabobs_glo(i, 3)<=ime) then
        if(tabobs_glo(i, 4)>=izs .and. tabobs_glo(i, 4)<=ize) then
          nb = nb + 1
          ! First columns for obs info
          tabobs(nb, 1:6) = tabobs_glo(i, 1:6)
          
          ! 7th column is for simulation or increment
          tabobs(nb, 7) = tabobs_glo(i, 7)
          
          ! 8th column is for reconstitution at the end, stores index
          tabobs(nb, 8) = dble(i)
          
          ! last columns stores fwd for tl and adj, not needed for fwd
          tabobs(nb, 9:12) = 0.
          if(fwd.eq.0) then
            tabobs(nb, 13) = tabobs_glo(i, 7)
          end if
        endif
      endif
    enddo
    ! send to workers
    call mpi_send(nobs(ip), 1, mpi_integer, ip, ias_nobs, mpi_comm_world, ierr)
    if(fwd.eq.1)call mpi_send(tabobs, nobs(ip) * 12, mpi_double_precision, ip, &
        ias_tabobs, mpi_comm_world, ierr)
    if(fwd.eq.0)call mpi_send(tabobs, nobs(ip) * 13, mpi_double_precision, ip, &
        ias_tabobs, mpi_comm_world, ierr)
    deallocate(tabobs)
  
  enddo
  ! output file for model values
  open(2, file = 'mod.txt', form = 'formatted', action = 'write')
  ! fin IP
  
  call mpi_barrier(mpi_comm_world, ierr)

end subroutine inidoms
