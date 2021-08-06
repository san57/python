subroutine htransport

    !  Calculation of horizontal first-order UPWIND flux rates for every mod
    !  cell. Using the first-order scheme, the transport terms, in the gover
    !  equation, are of the form:
    !  dCi/dt = SUM ( AijCj ) - BiCi
    !  The AijCj terms represent incoming fluxes, and the sum runs over some
    !  Of the adjacent boxes j of box i. The factors Aij are called the Flux
    !  rates, and are calculated in this subroutine, as well as the outgoing
    !  flux rates Bi.
    !  This routine is called every physical time step.
    !  INPUT : WINZLOC    Zonal current wind field
    !          WINMLOC    Meridian current wind field
    !          THLAYLOC   Current thicknesses of the model layers
    !          XSIZE      Longitudinal size of cells
    !          YSIZE      Latitudinal  size of cells
    !          XBASX      Normalized vector in the dir. of cell west  side:
    !          XBASY      Normalized vector in the dir. of cell west  side:
    !          YBASX      Normalized vector in the dir. of cell south side:
    !          YBASY      Normalized vector in the dir. of cell south side:
    !  OUTPUT: FLUXW      Flux rate at the western boundary
    !          FLUXE      Flux rate at the eastern boundary
    !          FLUXS      Flux rate at the southern boundary
    !          FLUXN      Flux rate at the northern boundary
    !          UWEST      Zonal wind component at the western cell side
    !          UEAST      Zonal wind component at the eastern cell side
    !          UNORTH     Meridional wind component at the northern cell sid
    !          USOUTH     Meridional wind component at the southern cell sid
    !          HWEST      Layer thicknesses at the western cell side
    !          HEAST      Layer thicknesses at the eastern cell side
    !          HNORTH     Layer thicknesses at the northern cell side
    !          HSOUTH     Layer thicknesses at the southern cell side

    use chimere_consts
    use worker_common

    implicit none

    !*****************************************************************************************
    integer :: izo, ime, ivert

    !*****************************************************************************************

    do ivert = 1, nverti
        do ime = 1, nmerid
            do izo = 1, nzonal
                !  Calculations of wind and layer thicknesses on cell sides
                !  These are just linearly interpolated

                if((izo.le.1).and.((dom_i==1))) then
                    uwest(izo, ime, ivert) = 1.5d0 * (winzloc(izo, ime, ivert) * xbasx(izo, ime)    &
                            + winmloc(izo, ime, ivert) * xbasy(izo, ime))                 &
                            - 0.5d0 * (winzloc(izo + 1, ime, ivert) * xbasx(izo + 1, ime)                &
                                    + winmloc(izo + 1, ime, ivert) * xbasy(izo + 1, ime))
                    hwest(izo, ime, ivert) = 1.5d0 * thlayloc(izo, ime, ivert) - 0.5d0 * thlayloc(izo + 1, ime, ivert)
                    swest(izo, ime, ivert) = 1.5d0 * ysize(izo, ime) - 0.5d0 * ysize(izo + 1, ime)
                else
                    uwest(izo, ime, ivert) = 0.5d0 * (winzloc(izo, ime, ivert) * xbasx(izo, ime)    &
                            + winmloc(izo, ime, ivert) * xbasy(izo, ime)                  &
                            + winzloc(izo - 1, ime, ivert) * xbasx(izo - 1, ime)                &
                            + winmloc(izo - 1, ime, ivert) * xbasy(izo - 1, ime))
                    hwest(izo, ime, ivert) = 0.5d0 * (thlayloc(izo, ime, ivert) + thlayloc(izo - 1, ime, ivert))
                    swest(izo, ime, ivert) = 0.5d0 * (ysize(izo, ime) + ysize(izo - 1, ime))
                endif

                if((izo.ge.nzonal).and.((dom_i==nzdoms))) then
                    ueast(izo, ime, ivert) = 1.5d0 * (winzloc(izo, ime, ivert) * xbasx(izo, ime)    &
                            + winmloc(izo, ime, ivert) * xbasy(izo, ime))                 &
                            - 0.5d0 * (winzloc(izo - 1, ime, ivert) * xbasx(izo - 1, ime)                &
                                    + winmloc(izo - 1, ime, ivert) * xbasy(izo - 1, ime))
                    heast(izo, ime, ivert) = 1.5d0 * thlayloc(izo, ime, ivert) - 0.5d0 * thlayloc(izo - 1, ime, ivert)
                    seast(izo, ime, ivert) = 1.5d0 * ysize(izo, ime) - 0.5d0 * ysize(izo - 1, ime)
                else
                    ueast(izo, ime, ivert) = 0.5d0 * (winzloc(izo, ime, ivert) * xbasx(izo, ime)    &
                            + winmloc(izo, ime, ivert) * xbasy(izo, ime)                  &
                            + winzloc(izo + 1, ime, ivert) * xbasx(izo + 1, ime)                &
                            + winmloc(izo + 1, ime, ivert) * xbasy(izo + 1, ime))
                    heast(izo, ime, ivert) = 0.5d0 * (thlayloc(izo, ime, ivert) + thlayloc(izo + 1, ime, ivert))
                    seast(izo, ime, ivert) = 0.5d0 * (ysize(izo, ime) + ysize(izo + 1, ime))
                endif

                if ((ime.le.1).and.((dom_j==1)))  then
                    usouth(izo, ime, ivert) = 1.5d0 * (winzloc(izo, ime, ivert) * ybasx(izo, ime)   &
                            + winmloc(izo, ime, ivert) * ybasy(izo, ime))                 &
                            - 0.5d0 * (winzloc(izo, ime + 1, ivert) * ybasx(izo, ime + 1)                &
                                    + winmloc(izo, ime + 1, ivert) * ybasy(izo, ime + 1))
                    hsouth(izo, ime, ivert) = 1.5d0 * thlayloc(izo, ime, ivert) - 0.5d0 * thlayloc(izo, ime + 1, ivert)
                    ssouth(izo, ime, ivert) = 1.5d0 * xsize(izo, ime) - 0.5d0 * xsize(izo, ime + 1)
                else
                    usouth(izo, ime, ivert) = 0.5d0 * (winzloc(izo, ime, ivert) * ybasx(izo, ime)   &
                            + winmloc(izo, ime, ivert) * ybasy(izo, ime)                  &
                            + winzloc(izo, ime - 1, ivert) * ybasx(izo, ime - 1)                &
                            + winmloc(izo, ime - 1, ivert) * ybasy(izo, ime - 1))
                    hsouth(izo, ime, ivert) = 0.5d0 * (thlayloc(izo, ime, ivert) + thlayloc(izo, ime - 1, ivert))
                    ssouth(izo, ime, ivert) = 0.5d0 * (xsize(izo, ime) + xsize(izo, ime - 1))
                endif

                if ((ime.ge.nmerid).and.((dom_j==nmdoms)))  then
                    unorth(izo, ime, ivert) = 1.5d0 * (winzloc(izo, ime, ivert) * ybasx(izo, ime)    &
                            + winmloc(izo, ime, ivert) * ybasy(izo, ime))                  &
                            - 0.5d0 * (winzloc(izo, ime - 1, ivert) * ybasx(izo, ime - 1)                 &
                                    + winmloc(izo, ime - 1, ivert) * ybasy(izo, ime - 1))
                    hnorth(izo, ime, ivert) = 1.5d0 * thlayloc(izo, ime, ivert) - 0.5d0 * thlayloc(izo, ime - 1, ivert)
                    snorth(izo, ime, ivert) = 1.5d0 * xsize(izo, ime) - 0.5d0 * xsize(izo, ime - 1)
                else
                    unorth(izo, ime, ivert) = 0.5d0 * (winzloc(izo, ime, ivert) * ybasx(izo, ime)    &
                            + winmloc(izo, ime, ivert) * ybasy(izo, ime)                   &
                            + winzloc(izo, ime + 1, ivert) * ybasx(izo, ime + 1)                 &
                            + winmloc(izo, ime + 1, ivert) * ybasy(izo, ime + 1))
                    hnorth(izo, ime, ivert) = 0.5d0 * (thlayloc(izo, ime, ivert) + thlayloc(izo, ime + 1, ivert))
                    snorth(izo, ime, ivert) = 0.5d0 * (xsize(izo, ime) + xsize(izo, ime + 1))
                endif

                !  Calculation of flux rates (independent of sign) on each cell side

                ! Western cell side
                fluxw(izo, ime, ivert) = uwest(izo, ime, ivert) / xsize(izo, ime)                 &
                        * hwest(izo, ime, ivert) / thlayloc(izo, ime, ivert)                        &
                        * swest(izo, ime, ivert) / ysize(izo, ime)
                ! Eastern cell side
                fluxe(izo, ime, ivert) = ueast(izo, ime, ivert) / xsize(izo, ime)                 &
                        * heast(izo, ime, ivert) / thlayloc(izo, ime, ivert)                        &
                        * seast(izo, ime, ivert) / ysize(izo, ime)
                ! Southern cell side
                fluxs(izo, ime, ivert) = usouth(izo, ime, ivert) / ysize(izo, ime)                &
                        * hsouth(izo, ime, ivert) / thlayloc(izo, ime, ivert)                       &
                        * ssouth(izo, ime, ivert) / xsize(izo, ime)
                ! Northern cell side
                fluxn(izo, ime, ivert) = unorth(izo, ime, ivert) / ysize(izo, ime)                &
                        * hnorth(izo, ime, ivert) / thlayloc(izo, ime, ivert)                       &
                        * snorth(izo, ime, ivert) / xsize(izo, ime)

            enddo
        enddo
    enddo

end subroutine htransport
