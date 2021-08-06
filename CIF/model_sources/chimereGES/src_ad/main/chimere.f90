program chimere

! main program for adjoint calculation with Chimere

  use chimere_common
  use master_message_subs
  implicit none
  integer :: ierr

  print*,rank,'  All readings and initializations'
  
  call init_mpi
  
  if(rank==0) then
    call inichimere
  else
    call ainiworker
  endif
  
  if(usechemistry.eq.0 .and. rank==0) print*,'do NOT use CHEMISTRY'
  if(usedepos.eq.0 .and. rank==0)     print*,'do NOT use DEPOSITION'
  if(useemissions.eq.0 .and. rank==0) print*,'do NOT use EMISSIONS'
  if(usetransmix.eq.0 .and. rank==0)  print*,'do NOT use TRANSMIX'
  if(usewetdepos.eq.0 .and. rank==0)  print*,'do NOT use WET DEPOSITION'
  if(useabsclipconc.eq.1 .and. rank==0) print*,'ALLOW negative concentrations'
  
  print*,rank,'  Run integration'

  if(rank==0) then
    call dirintegrun 
  else
    call worker
  endif
  
  print*,rank,' Initialize adjoint'
  
  if(rank==0)call ainichimere
  
  if(rank==0) then
    print*,rank,' INTEGRUN'
    call aintegrun
  else
    print*,rank,' WORKER'
    call aworker
  endif
 
  print*,rank,'Final Backup of all concentrations and clean-up'

  if (rank==0) then
    call endchimere
    call aendchimere
    call master_mpi_finalize
  else
    call worker_mpi_finalize
  endif

  ! End properly
  OPEN(1, file='all_good')
  CLOSE(1)
  
end program chimere

