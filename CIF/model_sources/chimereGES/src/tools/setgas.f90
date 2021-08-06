subroutine setgas(tgas,pgas,gasmw,denair,viscos,freemp)

  ! Compute air properties

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! INPUT
  !  tgas : temperature (K)
  !  pgas : pressure (Pa)
  ! OUTPUT
  !  denair : air density (Kg/m3)
  !  viscos : dynamic air viscosity (Kg/m/s)
  !  freemp : mean free path in air (m)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  implicit none

  ! subroutine arguments
  real(kind=8) :: tgas,pgas,gasmw,denair,viscos,freemp


  denair=1.21d-4*pgas*gasmw/tgas
  viscos=3.661d-03*tgas
  viscos=6.6164d-03*viscos*sqrt(viscos)/(tgas+1.14d+02)
  freemp=viscos/denair*sqrt(1.89d-04*gasmw/tgas)

end subroutine setgas
