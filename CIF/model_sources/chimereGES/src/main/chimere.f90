program chimere

!  Main program for simulation with CHIMERE

  use chimere_common
  use master_message_subs
  implicit none
  integer :: ierr
  
  !print*,rank,'  All readings and initializations'
  
  call init_mpi
  
  if(rank==0)  then
    call inichimere
  else
    !print*,'I, rank',rank,'am a worker'
    call iniworker
  endif

  if(usechemistry.eq.0 .and. rank==0) print*,'do NOT use CHEMISTRY'
  if(usedepos.eq.0 .and. rank==0)     print*,'do NOT use DEPOSITION'
  if(useemissions.eq.0 .and. rank==0) print*,'do NOT use EMISSIONS'
  if(usetransmix.eq.0 .and. rank==0)  print*,'do NOT use TRANSMIX'
  if(usewetdepos.eq.0 .and. rank==0)  print*,'do NOT use WET DEPOSITION'
  if(useabsclipconc.eq.1 .and. rank==0) print*,'ALLOW negative concentrations'

  !print*,rank,'  Run integration'
  if(rank==0) then
    call integrun
  else
    call worker
  endif
  !print*,rank,'fin run integration'

!  print*,rank,'Final Backup of all concentrations and clean-up'!
  if (rank==0) then
    call endchimere
    call master_mpi_finalize
  else
    call worker_mpi_finalize
  endif

  ! End properly
  OPEN(1, file='all_good')
  CLOSE(1)
  
end program chimere
