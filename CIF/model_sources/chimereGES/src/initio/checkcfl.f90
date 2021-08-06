subroutine checkcfl
!-----------------------------------------------------------------------
! For each run hour, the CFL number is estimated over the whole domain
! and in diagmet. The minimum nphour value is written in METEO.nc.
! If its value is greater than 1, the time step for the run
! decreases accordingly i.e the 'nphour' value increases.
! - nphour_ref is the nphour chosen by the user
! - nphourm is the nphour to force a CFL < 1
! - nphour is the value really used during the run
!-----------------------------------------------------------------------
  use netcdf
  use chimere_consts
  use chimere_common
  implicit none
!----------------------------------------------------------------
  if(nphour_ref.lt.nphourm(ihourrun))then
     nphour=nphourm(ihourrun)
  else
     nphour=nphour_ref
  endif
  dtr=3600./(nphour*ichemstep)
  dtr2=2d0*dtr/3d0
!----------------------------------------------------------------
end subroutine checkcfl
