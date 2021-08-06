subroutine renewhour_tl

  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

  !*****************************************************************************************

  !*****************************************************************************************

  !  Renewing temporary read variables at the hourly time step

  !  1: Meteo variables

    winz(:,:,:,1) = winz(:,:,:,2)
    winm(:,:,:,1) = winm(:,:,:,2)
    airm(:,:,:,1) = airm(:,:,:,2)
    clwc(:,:,:,1) = clwc(:,:,:,2)
    ! lmbb deepconv
    dpeu(:,:,:,1) = dpeu(:,:,:,2)
    dped(:,:,:,1) = dped(:,:,:,2)
    dpdu(:,:,:,1) = dpdu(:,:,:,2)
    dpdd(:,:,:,1) = dpdd(:,:,:,2)
     
    sphu(:,:,:,1) = sphu(:,:,:,2)
    temp(:,:,:,1) = temp(:,:,:,2)
    hlay(:,:,:,1) = hlay(:,:,:,2)
    kzzz(:,:,:,1) = kzzz(:,:,:,2)

    tem2(:,:,1) = tem2(:,:,2)
    sreh(:,:,1) = sreh(:,:,2)
    atte(:,:,1) = atte(:,:,2)
    hght(:,:,1) = hght(:,:,2)
    usta(:,:,1) = usta(:,:,2)
    aerr(:,:,1) = aerr(:,:,2)
    obuk(:,:,1) = obuk(:,:,2)
    wsta(:,:,1) = wsta(:,:,2)
    topc(:,:,1) = topc(:,:,2)

  !  2: Biogenic emissions
    emisb_tl(:,:,:,1) = emisb_tl(:,:,:,2)
    emisb(:,:,:,1) = emisb(:,:,:,2)

  !  3: Boundary conditions
    boundlat_tl(:,:,1) = boundlat_tl(:,:,2)
    boundlat(:,:,1) = boundlat(:,:,2)
    boundtop_tl(:,:,:,1) = boundtop_tl(:,:,:,2) 
    boundtop(:,:,:,1) = boundtop(:,:,:,2)

end subroutine renewhour_tl
