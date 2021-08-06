subroutine arenewhour

  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

  !*****************************************************************************************

  !*****************************************************************************************

  !  Renewing temporary read variables at the hourly time step

  !  1: Meteo variables

    winz(:,:,:,2) = winz(:,:,:,1)
    winm(:,:,:,2) = winm(:,:,:,1)
    airm(:,:,:,2) = airm(:,:,:,1)
    clwc(:,:,:,2) = clwc(:,:,:,1)
    ! lmbb deepconv
    dpeu(:,:,:,2) = dpeu(:,:,:,1)
    dped(:,:,:,2) = dped(:,:,:,1)
    dpdu(:,:,:,2) = dpdu(:,:,:,1)
    dpdd(:,:,:,2) = dpdd(:,:,:,1)
    
    sphu(:,:,:,2) = sphu(:,:,:,1)
    temp(:,:,:,2) = temp(:,:,:,1)
    hlay(:,:,:,2) = hlay(:,:,:,1)
    kzzz(:,:,:,2) = kzzz(:,:,:,1)
    tem2(:,:,2) = tem2(:,:,1)
    sreh(:,:,2) = sreh(:,:,1)
    atte(:,:,2) = atte(:,:,1)
    hght(:,:,2) = hght(:,:,1)
    usta(:,:,2) = usta(:,:,1)
    aerr(:,:,2) = aerr(:,:,1)
    obuk(:,:,2) = obuk(:,:,1)
    wsta(:,:,2) = wsta(:,:,1)
    topc(:,:,2) = topc(:,:,1)

  !  2: Biogenic emissions

    emisb(:,:,:,2) = emisb(:,:,:,1)

  !  3: Boundary conditions

    boundlat(:,:,2) = boundlat(:,:,1)
    boundtop(:,:,:,2) = boundtop(:,:,:,1)

end subroutine arenewhour
