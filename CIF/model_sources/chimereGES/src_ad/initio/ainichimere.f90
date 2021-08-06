subroutine ainichimere
  
  
  ! 1 : initialisation of physics
  call ainiphys
  
  ! 2 : initialisation of emissions
  call ainiemis
  
  ! 3 : Initialization of adj emissions
  call ainiaemis
  
  ! 3 bis: Initialization of adj immissions
  call ainiaboun
  
  ! 4 : Initial of adj concentrations + adj ini conc
  call ainiaconc
  
  ! 5: Initialization of netCDF output files
  call ainiout
  call ainiend

end subroutine ainichimere
