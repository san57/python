subroutine mixing

    !  Calculates instantaneous mixing parameters. These are then used
    !  in VTRANSPORT.
    !
    !  INPUT :  HGHTLOC   Current height of the mixed layer
    !           HLAYLOC   Current layer top altitudes
    !           THLAYLOC  Thicknesses of the layers
    !           KZZZLOC   Current Diffusivity
    !  OUTPUT:  WINXLOC   Velocity-equivalent of mixing exchange rate

    use worker_common

    implicit none


    !*****************************************************************************
    integer :: izo, ime, ivert
    real(kind = 8) :: thlayer, thupper, thmid


    !*****************************************************************************
    do ivert = 1, nverti
        do ime = 1, nmerid
            do izo = 1, nzonal

                !  THUPPER is the upper layer thickness. For the top layer, TMPL is kept
                !  as the thickness of the top layer
                !  THLAYER is the layer thickness
                !  THMID   is the thickness around interface

                thlayer = thlayloc(izo, ime, ivert)
                if(ivert.lt.nverti) then
                    thupper = thlayloc(izo, ime, ivert + 1)
                else
                    thupper = thlayer
                endif
                thmid = 5d-1 * (thlayer + thupper)
                winxloc(izo, ime, ivert) = kzzzloc(izo, ime, ivert) / thmid
            end do
        end do
    end do

end subroutine mixing
