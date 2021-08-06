subroutine inichimere


  use chimere_consts
  use chimere_common
  use master_message_subs

  !  Initialization of the model
  
  !print*,'  1a: Reading filenames'
  call iniread
 
  !print*,' initial time-step values'
  dtr=3600d0/(nphour_ref*ichemstep)
  dtr2=2d0*dtr/3d0
  nphour=nphour_ref

  ! Now that we know dynamic dimensions, we can allocate storage
  !print*,'debut master allocall'
  call master_allocall
  !print*,'fin master allocall'
  call master_allocwhole
  if(ad.eq.1) call  master_allocwhole_ad

  !print*,'  1aa: Initialization of input netcdf files that remain open during all the run'
  call inicdf

  !print*,'  2: Initialization of geometry'
  call inigeom
  psurf = bvcoord(1)*100000. ! Psurf in Pa

  !print*,'  3a: Initialization of chemistry'
  call inichem

  !print*,'  4: Initialization of deposition parameters'
  call inidepo

  !print*,'  5: Initialization of physics'
  call iniphys

  !print*,'  6: Initialization of emissions'
  call iniemis

  !print*,'  7: Initialization of immissions'
  call iniboun

  !print*,'  7bis: Initialization of MPI domains'
  call inidoms

  !print*,'  8: Initial of concentrations'
  call iniconc

  if (ad.eq.0) then
    !print*,'  9: Initialization of netCDF output files',NetCDF_output
    if (NetCDF_output.eq.1) then
      print*,'Creating netCDF output files'
      call iniout
      call iniparam
    endif
    !print*,'10: call ending procedures'
    call iniend
    call inidepout
  endif
  
end subroutine inichimere
