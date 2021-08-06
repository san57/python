subroutine atransmix(ns, izo, ime, ivert, atrpr, atrlo)

    !!!!!!!!!!!ADJOINT of transmix !!!!!!!!!!!!!!!!!!!!!!

    !  Production/loss due to Transport and mixing
    !  INPUT : NS           Current species number
    !          IZO,IME,IVERT           Current cell coordinates
    !          IBWEST       Western cell
    !          IBEAST       Eastern cell
    !          IBSOUTH      Southern cell
    !          IBNORTH      Northern cell
    !          CONC         Current concentration array
    !          VFLUXO       Outgoing vertical fluxes
    !          VFLUXI       Incoming vertical fluxes
    !          UWEST        Western zonal wind
    !          UEAST        Eastern zonal wind
    !          USOUTH       Southern meridional wind
    !          UNORTH       Northern meridional wind
    !          FLUXW        Western flux rate
    !          FLUXE        Eastern flux rate
    !          FLUXN        Northern flux rate
    !          FLUXS        Southern flux rate
    !  OUTPUT: TRPR         Transport production fluxes
    !          TRLO         Transport loss fluxes

    use chimere_consts
    use worker_common

    implicit none

    include 'mpif.h'


    !***************************************************************************
    ! subroutine arguments
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: atrpr
    real(kind = 8) :: atrlo

    ! local variables
    integer :: inp, inpv
    integer :: nv
    real(kind = 8) :: cup, cupup, cdo, cprime
    real(kind = 8) :: concu, concd, flxe
    real(kind = 8) :: rhoup, rhodo, dxupup, dxup, dxdo, dens_bound

    real(kind = 8) :: acup, acupup, acdo, acprime
    real(kind = 8) :: aconcu, aconcd, aflxe

    ! external functions
    real(kind = 8) :: vanleer, vanleer_nonunif
    real(kind = 8) :: dvanleer, dvanleer_nonunif

    !***************************************************************************
    inp = species(ns)%transp
    inpv = species(ns)%transpv

    !!!!!!!!!!!!!!!! ADJOINT - horizontal transport !!!!!!!!!!!!!!!!
    if(ivert==nverti) then
        aflxeconc(ns, izo, ime, :) = dzero
        aflxuconc(ns, izo, ime, :) = dzero
    endif
    aflxdconc(ns, izo, ime, :) = dzero

    !  Horizontal transport
    !  4: Northern side

    if(unorth(izo, ime, ivert).lt.dzero) then
        if(inp.eq.1) then
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo, ime + 2, ivert) / airmloc(izo, ime + 2, ivert)
            cup = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime =vanleer(cupup,cup,cdo,ysize(izo,ime+1),-unorth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime + 2), &
                    &      ysize(izo, ime + 1), ysize(izo, ime), -unorth(izo, ime, ivert), dtr)

            acprime = -atrpr * fluxn(izo, ime, ivert) * airmloc(izo, ime + 1, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,ysize(izo,ime+1),-unorth(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,ysize(izo,ime+1),-unorth(izo,ime,ivert),dtr)
            !acupup =acprime* dvanleer(1,cupup,cup,cdo,ysize(izo,ime+1),-unorth(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, ysize(izo, ime + 2), &
                    &       ysize(izo, ime + 1), ysize(izo, ime), -unorth(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, ysize(izo, ime + 2), &
                    &       ysize(izo, ime + 1), ysize(izo, ime), -unorth(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, ysize(izo, ime + 2), &
                    &       ysize(izo, ime + 1), ysize(izo, ime), -unorth(izo, ime, ivert), dtr)

            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &          acdo / airmloc(izo, ime, ivert)
            aconc(ns, izo, ime + 1, ivert) = aconc(ns, izo, ime + 1, ivert) + &
                    &          acup / airmloc(izo, ime + 1, ivert)
            aconc(ns, izo, ime + 2, ivert) = aconc(ns, izo, ime + 2, ivert) + &
                    &          acupup / airmloc(izo, ime + 2, ivert)

        else
            ! UPWIND
            aconc(ns, izo, ime + 1, ivert) = aconc(ns, izo, ime + 1, ivert) &
                    & - fluxn(izo, ime, ivert) * atrpr

        endif
    else
        if(inp.eq.1) then
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            !cprime = vanleer(cupup,cup,cdo,ysize(izo,ime),unorth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime - 1), &
                    &     ysize(izo, ime), ysize(izo, ime + 1), unorth(izo, ime, ivert), dtr)
            acprime = atrlo * fluxn(izo, ime, ivert) * airmloc(izo, ime, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,ysize(izo,ime),unorth(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,ysize(izo,ime),unorth(izo,ime,ivert),dtr)
            !acupup = acprime*dvanleer(1,cupup,cup,cdo,ysize(izo,ime),unorth(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, ysize(izo, ime - 1), &
                    & ysize(izo, ime), ysize(izo, ime + 1), unorth(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, ysize(izo, ime - 1), &
                    & ysize(izo, ime), ysize(izo, ime + 1), unorth(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, ysize(izo, ime - 1), &
                    & ysize(izo, ime), ysize(izo, ime + 1), unorth(izo, ime, ivert), dtr)

            aconc(ns, izo, ime + 1, ivert) = aconc(ns, izo, ime + 1, ivert) + &
                    &         acdo / airmloc(izo, ime + 1, ivert)
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &         acup / airmloc(izo, ime, ivert)
            aconc(ns, izo, ime - 1, ivert) = aconc(ns, izo, ime - 1, ivert) + &
                    &         acupup / airmloc(izo, ime - 1, ivert)

        else
            ! UPWIND
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &         fluxn(izo, ime, ivert) * atrlo

        endif
    endif


    !  3: Southern side

    if(usouth(izo, ime, ivert).gt.dzero) then
        if(inp.eq.1) then
            ! PPM
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo, ime - 2, ivert) / airmloc(izo, ime - 2, ivert)
            cup = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,ysize(izo,ime-1),usouth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime - 2), &
                    & ysize(izo, ime - 1), ysize(izo, ime), usouth(izo, ime, ivert), dtr)
            acprime = atrpr * fluxs(izo, ime, ivert) * airmloc(izo, ime - 1, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,ysize(izo,ime-1),usouth(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,ysize(izo,ime-1),usouth(izo,ime,ivert),dtr)
            !acupup = acprime*dvanleer(1,cupup,cup,cdo,ysize(izo,ime-1),usouth(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, ysize(izo, ime - 2), &
                    &    ysize(izo, ime - 1), ysize(izo, ime), usouth(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, ysize(izo, ime - 2), &
                    &    ysize(izo, ime - 1), ysize(izo, ime), usouth(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, ysize(izo, ime - 2), &
                    &    ysize(izo, ime - 1), ysize(izo, ime), usouth(izo, ime, ivert), dtr)

            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &         acdo / airmloc(izo, ime, ivert)
            aconc(ns, izo, ime - 1, ivert) = aconc(ns, izo, ime - 1, ivert) + &
                    &         acup / airmloc(izo, ime - 1, ivert)
            aconc(ns, izo, ime - 2, ivert) = aconc(ns, izo, ime - 2, ivert) + &
                    &         acupup / airmloc(izo, ime - 2, ivert)

        else
            ! UPWIND
            aconc(ns, izo, ime - 1, ivert) = aconc(ns, izo, ime - 1, ivert) + &
                    &         fluxs(izo, ime, ivert) * atrpr

        endif
    else
        if(inp.eq.1) then
            ! PPM
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            !cprime = vanleer(cupup,cup,cdo,ysize(izo,ime),-usouth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime + 1), &
                    & ysize(izo, ime), ysize(izo, ime - 1), -usouth(izo, ime, ivert), dtr)
            acprime = - atrlo * fluxs(izo, ime, ivert) * airmloc(izo, ime, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,ysize(izo,ime),-usouth(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,ysize(izo,ime),-usouth(izo,ime,ivert),dtr)
            !acupup = acprime*dvanleer(1,cupup,cup,cdo,ysize(izo,ime),-usouth(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, ysize(izo, ime + 1), &
                    &    ysize(izo, ime), ysize(izo, ime - 1), -usouth(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, ysize(izo, ime + 1), &
                    &    ysize(izo, ime), ysize(izo, ime - 1), -usouth(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, ysize(izo, ime + 1), &
                    &    ysize(izo, ime), ysize(izo, ime - 1), -usouth(izo, ime, ivert), dtr)

            aconc(ns, izo, ime - 1, ivert) = aconc(ns, izo, ime - 1, ivert) + &
                    &         acdo / airmloc(izo, ime - 1, ivert)
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &         acup / airmloc(izo, ime, ivert)
            aconc(ns, izo, ime + 1, ivert) = aconc(ns, izo, ime + 1, ivert) + &
                    &         acupup / airmloc(izo, ime + 1, ivert)

        else
            ! UPWIND
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                    & - fluxs(izo, ime, ivert) * atrlo
        endif
    endif

    !  2: Eastern side

    if(ueast(izo, ime, ivert).lt.dzero) then

        if(inp.eq.1) then
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo + 2, ime, ivert) / airmloc(izo + 2, ime, ivert)
            cup = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,xsize(izo+1,ime),-ueast(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo + 2, ime), &
                    &       xsize(izo + 1, ime), xsize(izo, ime), -ueast(izo, ime, ivert), dtr)
            acprime = -atrpr * fluxe(izo, ime, ivert) * airmloc(izo + 1, ime, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,xsize(izo+1,ime),-ueast(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,xsize(izo+1,ime),-ueast(izo,ime,ivert),dtr)
            !acupup = acprime*dvanleer(1,cupup,cup,cdo,xsize(izo+1,ime),-ueast(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, xsize(izo + 2, ime), &
                    &     xsize(izo + 1, ime), xsize(izo, ime), -ueast(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, xsize(izo + 2, ime), &
                    &     xsize(izo + 1, ime), xsize(izo, ime), -ueast(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, xsize(izo + 2, ime), &
                    &     xsize(izo + 1, ime), xsize(izo, ime), -ueast(izo, ime, ivert), dtr)

            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &          acdo / airmloc(izo, ime, ivert)
            aconc(ns, izo + 1, ime, ivert) = aconc(ns, izo + 1, ime, ivert) + &
                    &          acup / airmloc(izo + 1, ime, ivert)
            aconc(ns, izo + 2, ime, ivert) = aconc(ns, izo + 2, ime, ivert) + &
                    &          acupup / airmloc(izo + 2, ime, ivert)

        else
            ! UPWIND
            aconc(ns, izo + 1, ime, ivert) = aconc(ns, izo + 1, ime, ivert) &
                    & - fluxe(izo, ime, ivert) * atrpr

        endif
    else
        if(inp.eq.1) then
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,xsize(izo,ime),ueast(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo - 1, ime), &
                    &      xsize(izo, ime), xsize(izo + 1, ime), ueast(izo, ime, ivert), dtr)
            acprime = atrlo * fluxe(izo, ime, ivert) * airmloc(izo, ime, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,xsize(izo,ime),ueast(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,xsize(izo,ime),ueast(izo,ime,ivert),dtr)
            !acupup = acprime*dvanleer(1,cupup,cup,cdo,xsize(izo,ime),ueast(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, xsize(izo - 1, ime), &
                    &     xsize(izo, ime), xsize(izo + 1, ime), ueast(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, xsize(izo - 1, ime), &
                    &     xsize(izo, ime), xsize(izo + 1, ime), ueast(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, xsize(izo - 1, ime), &
                    &     xsize(izo, ime), xsize(izo + 1, ime), ueast(izo, ime, ivert), dtr)

            aconc(ns, izo + 1, ime, ivert) = aconc(ns, izo + 1, ime, ivert) + &
                    &          acdo / airmloc(izo + 1, ime, ivert)
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &          acup / airmloc(izo, ime, ivert)
            aconc(ns, izo - 1, ime, ivert) = aconc(ns, izo - 1, ime, ivert) + &
                    &          acupup / airmloc(izo - 1, ime, ivert)

        else
            ! UPWIND
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &         fluxe(izo, ime, ivert) * atrlo

        endif
    endif

    !  1: Western side
    if(uwest(izo, ime, ivert).gt.dzero) then
        if(inp.eq.1) then
            ! PPM
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo - 2, ime, ivert) / airmloc(izo - 2, ime, ivert)
            cup = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,xsize(izo-1,ime),uwest(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo - 2, ime), &
                    &         xsize(izo - 1, ime), xsize(izo, ime), uwest(izo, ime, ivert), dtr)
            acprime = atrpr * fluxw(izo, ime, ivert) * airmloc(izo - 1, ime, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,xsize(izo-1,ime),uwest(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,xsize(izo-1,ime),uwest(izo,ime,ivert),dtr)
            !acupup = acprime*dvanleer(1,cupup,cup,cdo,xsize(izo-1,ime),uwest(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, xsize(izo - 2, ime), &
                    &      xsize(izo - 1, ime), xsize(izo, ime), uwest(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, xsize(izo - 2, ime), &
                    &      xsize(izo - 1, ime), xsize(izo, ime), uwest(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, xsize(izo - 2, ime), &
                    &      xsize(izo - 1, ime), xsize(izo, ime), uwest(izo, ime, ivert), dtr)

            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &         acdo / airmloc(izo, ime, ivert)
            aconc(ns, izo - 1, ime, ivert) = aconc(ns, izo - 1, ime, ivert) + &
                    &         acup / airmloc(izo - 1, ime, ivert)
            aconc(ns, izo - 2, ime, ivert) = aconc(ns, izo - 2, ime, ivert) + &
                    &         acupup / airmloc(izo - 2, ime, ivert)

        else
            ! UPWIND
            aconc(ns, izo - 1, ime, ivert) = aconc(ns, izo - 1, ime, ivert) + &
                    &         fluxw(izo, ime, ivert) * atrpr
        endif
    else
        if(inp.eq.1) then
            ! PPM
            print*, '******* PPM scheme not adjointized **********'
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,xsize(izo,ime),-uwest(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo + 1, ime), &
                    &      xsize(izo, ime), xsize(izo - 1, ime), -uwest(izo, ime, ivert), dtr)
            acprime = -atrlo * fluxw(izo, ime, ivert) * airmloc(izo, ime, ivert)
            !acdo = acprime*dvanleer(3,cupup,cup,cdo,xsize(izo,ime),-uwest(izo,ime,ivert),dtr)
            !acup = acprime*dvanleer(2,cupup,cup,cdo,xsize(izo,ime),-uwest(izo,ime,ivert),dtr)
            !acupup = acprime*dvanleer(1,cupup,cup,cdo,xsize(izo,ime),-uwest(izo,ime,ivert),dtr)
            acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, xsize(izo + 1, ime), &
                    &    xsize(izo, ime), xsize(izo - 1, ime), -uwest(izo, ime, ivert), dtr)
            acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, xsize(izo + 1, ime), &
                    &    xsize(izo, ime), xsize(izo - 1, ime), -uwest(izo, ime, ivert), dtr)
            acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, xsize(izo + 1, ime), &
                    &    xsize(izo, ime), xsize(izo - 1, ime), -uwest(izo, ime, ivert), dtr)

            aconc(ns, izo - 1, ime, ivert) = aconc(ns, izo - 1, ime, ivert) + &
                    &         acdo / airmloc(izo - 1, ime, ivert)
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                    &         acup / airmloc(izo, ime, ivert)
            aconc(ns, izo + 1, ime, ivert) = aconc(ns, izo + 1, ime, ivert) + &
                    &         acupup / airmloc(izo + 1, ime, ivert)

        else
            ! UPWIND
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                    & - fluxw(izo, ime, ivert) * atrlo
        endif
    endif


    !  print*,' DEEP CONVECTION IMPLEMENTATION !!!!!!!! ADJOINT !!!!!!!!!
    if(ideepconv.ne.0.and.ideep(izo, ime).eq.1)then
        !********************************************************************************
        !!!!!!!! direct calculation needed for internal variables like concu,flxe,...!!!!
        !********************************************************************************

        ! UPDRAUGHT FLUX
        if((flxuloc(izo, ime, ivert) + dpduloc(izo, ime, ivert)).le.dzero) then
            concu = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
        else
            if(ivert.eq.1) then
                concu = dpeuloc(izo, ime, ivert) * conc(ns, izo, ime, ivert) &
                        / airmloc(izo, ime, ivert) &
                        / (flxuloc(izo, ime, ivert) + dpduloc(izo, ime, ivert))
            else
                concu = (flxuconc(ns, izo, ime, ivert - 1)   &
                        + dpeuloc(izo, ime, ivert) * conc(ns, izo, ime, ivert) &
                                / airmloc(izo, ime, ivert)) &
                        / (flxuloc(izo, ime, ivert) + dpduloc(izo, ime, ivert))
            endif
        endif
        flxuconc(ns, izo, ime, ivert) = flxuloc(izo, ime, ivert) * concu

        ! DOWNDRAUGHT FLUX
        flxdconc(ns, izo, ime, :) = dzero
        do nv = nverti, ivert, -1
            if(nv.gt.1) then
                if((-flxdloc(izo, ime, nv - 1) + dpddloc(izo, ime, nv)).le.dzero) then
                    concd = conc(ns, izo, ime, nv) / airmloc(izo, ime, nv)
                else
                    concd = (flxdconc(ns, izo, ime, nv) &
                            - dpedloc(izo, ime, nv) * conc(ns, izo, ime, nv) / airmloc(izo, ime, nv)) &
                            / (flxdloc(izo, ime, nv - 1) - dpddloc(izo, ime, nv))
                endif
                flxdconc(ns, izo, ime, nv - 1) = flxdloc(izo, ime, nv - 1) * concd
            endif
        enddo

        flxe = - flxuloc(izo, ime, ivert) - flxdloc(izo, ime, ivert)
        if(ivert.ne.nverti)then
            if(flxe.lt.dzero)then
                flxeconc(ns, izo, ime, ivert) = flxe * conc(ns, izo, ime, ivert + 1) / airmloc(izo, ime, ivert + 1)
            else
                flxeconc(ns, izo, ime, ivert) = flxe * conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            endif
        else
            flxeconc(ns, izo, ime, ivert) = dzero
        endif
        !!! end needed direct calculation

        !!!!!!!!!!!!!! adjoint

        ! SOLVING TRANSMIX PRODUCTION AND LOSS
        if(ivert.eq.1)then
            if(flxeconc(ns, izo, ime, ivert).lt.dzero)then
                aflxeconc(ns, izo, ime, ivert) = aflxeconc(ns, izo, ime, ivert) &
                        - atrpr / thlayloc(izo, ime, ivert)
            else
                aflxeconc(ns, izo, ime, ivert) = aflxeconc(ns, izo, ime, ivert) &
                        + atrlo / thlayloc(izo, ime, ivert)
            endif
            aflxuconc(ns, izo, ime, ivert) = aflxuconc(ns, izo, ime, ivert) + atrlo / thlayloc(izo, ime, ivert)
            aflxdconc(ns, izo, ime, ivert) = -atrpr / thlayloc(izo, ime, ivert)
        else
            if((flxeconc(ns, izo, ime, ivert) - flxeconc(ns, izo, ime, ivert - 1)).lt.dzero)then
                aflxeconc(ns, izo, ime, ivert) = aflxeconc(ns, izo, ime, ivert) &
                        - atrpr / thlayloc(izo, ime, ivert)
                aflxeconc(ns, izo, ime, ivert - 1) = aflxeconc(ns, izo, ime, ivert - 1) &
                        + atrpr / thlayloc(izo, ime, ivert)
            else
                aflxeconc(ns, izo, ime, ivert) = aflxeconc(ns, izo, ime, ivert) &
                        + atrlo / thlayloc(izo, ime, ivert)
                aflxeconc(ns, izo, ime, ivert - 1) = aflxeconc(ns, izo, ime, ivert - 1) &
                        - atrlo / thlayloc(izo, ime, ivert)
            endif
            if((flxdconc(ns, izo, ime, ivert) - flxdconc(ns, izo, ime, ivert - 1)).lt.dzero)then
                aflxdconc(ns, izo, ime, ivert) = aflxdconc(ns, izo, ime, ivert) &
                        - atrpr / thlayloc(izo, ime, ivert)
                aflxdconc(ns, izo, ime, ivert - 1) = aflxdconc(ns, izo, ime, ivert - 1) &
                        + atrpr / thlayloc(izo, ime, ivert)
            else
                aflxdconc(ns, izo, ime, ivert) = aflxdconc(ns, izo, ime, ivert) &
                        + atrlo / thlayloc(izo, ime, ivert)
                aflxdconc(ns, izo, ime, ivert - 1) = aflxdconc(ns, izo, ime, ivert - 1) &
                        - atrlo / thlayloc(izo, ime, ivert)
            endif
            if((flxuconc(ns, izo, ime, ivert) - flxuconc(ns, izo, ime, ivert - 1)).lt.dzero)then
                aflxuconc(ns, izo, ime, ivert) = aflxuconc(ns, izo, ime, ivert) &
                        - atrpr / thlayloc(izo, ime, ivert)
                aflxuconc(ns, izo, ime, ivert - 1) = aflxuconc(ns, izo, ime, ivert - 1) &
                        + atrpr / thlayloc(izo, ime, ivert)
            else
                aflxuconc(ns, izo, ime, ivert) = aflxuconc(ns, izo, ime, ivert) &
                        + atrlo / thlayloc(izo, ime, ivert)
                aflxuconc(ns, izo, ime, ivert - 1) = aflxuconc(ns, izo, ime, ivert - 1) &
                        - atrlo / thlayloc(izo, ime, ivert)
            endif
        endif ! ivert==1

        ! ENTRAINMENT FLUX
        if(ivert.ne.nverti)then
            if(flxe.lt.dzero)then
                aconc(ns, izo, ime, ivert + 1) = aconc(ns, izo, ime, ivert + 1) &
                        + flxe / airmloc(izo, ime, ivert + 1) * aflxeconc(ns, izo, ime, ivert)
            else
                aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                        + flxe / airmloc(izo, ime, ivert) * aflxeconc(ns, izo, ime, ivert)
            endif
        endif

        ! DOWNDRAUGHT FLUX
        do nv = ivert, nverti
            if(nv.gt.1) then
                aconcd = aflxdconc(ns, izo, ime, nv - 1) * flxdloc(izo, ime, nv - 1)
                if((-flxdloc(izo, ime, nv - 1) + dpddloc(izo, ime, nv)).le.dzero) then
                    aconc(ns, izo, ime, nv) = aconc(ns, izo, ime, nv) &
                            + aconcd / airmloc(izo, ime, nv)
                else
                    aconc(ns, izo, ime, nv) = aconc(ns, izo, ime, nv) &
                            - aconcd * dpedloc(izo, ime, nv) / airmloc(izo, ime, nv) &
                                    / (flxdloc(izo, ime, nv - 1) - dpddloc(izo, ime, nv))
                    aflxdconc(ns, izo, ime, nv) = aflxdconc(ns, izo, ime, nv) &
                            + aconcd / (flxdloc(izo, ime, nv - 1) - dpddloc(izo, ime, nv))
                endif
            endif
        enddo

        ! UPDRAUGHT FLUX
        aconcu = flxuloc(izo, ime, ivert) * aflxuconc(ns, izo, ime, ivert)
        if((flxuloc(izo, ime, ivert) + dpduloc(izo, ime, ivert)).le.dzero) then
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                    + aconcu / airmloc(izo, ime, ivert)
        else
            if (ivert.eq.1) then
                aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                        + aconcu * dpeuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) &
                                / (flxuloc(izo, ime, ivert) + dpduloc(izo, ime, ivert))
            else
                aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                        + aconcu * dpeuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) &
                                / (flxuloc(izo, ime, ivert) + dpduloc(izo, ime, ivert))
                aflxuconc(ns, izo, ime, ivert - 1) = aflxuconc(ns, izo, ime, ivert - 1) + &
                        aconcu / (flxuloc(izo, ime, ivert) + dpduloc(izo, ime, ivert))
            endif
        endif

    endif ! of ideepconv=1
    ! DEEP CONVECTION IMPLEMENTATION
    !  Vertical transport

    ! 3 : Transport from resolved vertical velocity
    if(inpv.eq.1) then
        ! UPWIND
        ! up side
        if(winvloc(izo, ime, ivert).gt.dzero)then
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                    + winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * atrlo
        endif
        if(winvloc(izo, ime, ivert).lt.dzero)then
            aconc(ns, izo, ime, ivert + 1) = aconc(ns, izo, ime, ivert + 1) &
                    - winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * atrpr
        endif
        !down side
        if(ivert.gt.1)then
            if(winvloc(izo, ime, ivert - 1).gt.dzero)then
                aconc(ns, izo, ime, ivert - 1) = aconc(ns, izo, ime, ivert - 1) &
                        + winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * atrpr
            endif
            if(winvloc(izo, ime, ivert - 1).lt.dzero)then
                aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                        - winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * atrlo
            endif
        endif
    else
        if(inpv.eq.2) then
            ! VAN LEER
            !upper side
            if(ivert.lt.nverti)then
                !vanleer
                if(winvloc(izo, ime, ivert).lt.dzero)then !incoming flux
                    if(ivert.lt.nverti - 1)then
                        cupup = conc(ns, izo, ime, ivert + 2) / airmloc(izo, ime, ivert + 2)
                        cup = conc(ns, izo, ime, ivert + 1) / airmloc(izo, ime, ivert + 1)
                        cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
                        rhoup = airmloc(izo, ime, ivert + 1)
                        rhodo = airmloc(izo, ime, ivert)
                        dxupup = thlayloc(izo, ime, ivert + 2)
                        dxup = thlayloc(izo, ime, ivert + 1)
                        dxdo = thlayloc(izo, ime, ivert)
                        cprime = vanleer_nonunif(cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert)), dtr)
                        dens_bound = (rhoup * dxdo + rhodo * dxup) / (dxup + dxdo)
                        acprime = -atrpr * winvloc(izo, ime, ivert) * dens_bound / thlayloc(izo, ime, ivert)
                        acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert)), dtr)
                        acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert)), dtr)
                        acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert)), dtr)
                        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                                acdo / airmloc(izo, ime, ivert)
                        aconc(ns, izo, ime, ivert + 1) = aconc(ns, izo, ime, ivert + 1) + &
                                acup / airmloc(izo, ime, ivert + 1)
                        aconc(ns, izo, ime, ivert + 2) = aconc(ns, izo, ime, ivert + 2) + &
                                acupup / airmloc(izo, ime, ivert + 2)

                    else  ! upwind, because only one cell
                        aconc(ns, izo, ime, ivert + 1) = aconc(ns, izo, ime, ivert + 1) &
                                - winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * atrpr
                    endif ! ivert.lt.nverti-1
                else !outgoing flux
                    if(ivert.gt.1)then
                        cupup = conc(ns, izo, ime, ivert - 1) / airmloc(izo, ime, ivert - 1)
                        cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
                        cdo = conc(ns, izo, ime, ivert + 1) / airmloc(izo, ime, ivert + 1)
                        rhoup = airmloc(izo, ime, ivert)
                        rhodo = airmloc(izo, ime, ivert + 1)
                        dxupup = thlayloc(izo, ime, ivert - 1)
                        dxup = thlayloc(izo, ime, ivert)
                        dxdo = thlayloc(izo, ime, ivert + 1)
                        cprime = vanleer_nonunif(cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert), dtr)
                        dens_bound = (rhoup * dxdo + rhodo * dxup) / (dxup + dxdo)
                        acprime = atrlo * winvloc(izo, ime, ivert) * dens_bound / thlayloc(izo, ime, ivert)
                        acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert), dtr)
                        acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert), dtr)
                        acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert), dtr)
                        aconc(ns, izo, ime, ivert + 1) = aconc(ns, izo, ime, ivert + 1) + &
                                acdo / airmloc(izo, ime, ivert + 1)
                        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                                acup / airmloc(izo, ime, ivert)
                        aconc(ns, izo, ime, ivert - 1) = aconc(ns, izo, ime, ivert - 1) + &
                                acupup / airmloc(izo, ime, ivert - 1)
                    else
                        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                                + winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * atrlo
                    endif !ivert.gt.1
                endif ! winvloc(izo,ime,ivert).lt.0
            else ! ivert.lt.nverti
                !upwind, because only one cell
                if(winvloc(izo, ime, ivert).gt.dzero) then
                    aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                            + winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * atrlo
                endif
                if(winvloc(izo, ime, ivert).lt.dzero) then
                    aconc(ns, izo, ime, ivert + 1) = aconc(ns, izo, ime, ivert + 1) &
                            - winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * atrpr
                endif
            endif !end treatment of upper side
            !lower side
            if(ivert.gt.1)then
                if(winvloc(izo, ime, ivert - 1).gt.dzero)then
                    !incoming flux from lower side
                    if(ivert.gt.2)then
                        !vanleer
                        cupup = conc(ns, izo, ime, ivert - 2) / airmloc(izo, ime, ivert - 2)
                        cup = conc(ns, izo, ime, ivert - 1) / airmloc(izo, ime, ivert - 1)
                        cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
                        rhoup = airmloc(izo, ime, ivert - 1)
                        rhodo = airmloc(izo, ime, ivert)
                        dxupup = thlayloc(izo, ime, ivert - 2)
                        dxup = thlayloc(izo, ime, ivert - 1)
                        dxdo = thlayloc(izo, ime, ivert)
                        cprime = vanleer_nonunif(cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert - 1), dtr)
                        dens_bound = (rhoup * dxdo + rhodo * dxup) / (dxup + dxdo)

                        acprime = atrpr * winvloc(izo, ime, ivert - 1) * dens_bound / thlayloc(izo, ime, ivert)
                        acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert - 1), dtr)
                        acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert - 1), dtr)
                        acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, dxupup, dxup, dxdo, winvloc(izo, ime, ivert - 1), dtr)
                        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                                acdo / airmloc(izo, ime, ivert)
                        aconc(ns, izo, ime, ivert - 1) = aconc(ns, izo, ime, ivert - 1) + &
                                acup / airmloc(izo, ime, ivert - 1)
                        aconc(ns, izo, ime, ivert - 2) = aconc(ns, izo, ime, ivert - 2) + &
                                acupup / airmloc(izo, ime, ivert - 2)
                    else ! upwind, because only one cell below
                        aconc(ns, izo, ime, ivert - 1) = aconc(ns, izo, ime, ivert - 1) &
                                + winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * atrpr
                    endif !of ivert.gt.2
                else
                    !outgoing flux on lower side
                    if(ivert.lt.nverti)then
                        !vanleer
                        cupup = conc(ns, izo, ime, ivert + 1) / airmloc(izo, ime, ivert + 1)
                        cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
                        cdo = conc(ns, izo, ime, ivert - 1) / airmloc(izo, ime, ivert - 1)
                        rhoup = airmloc(izo, ime, ivert)
                        rhodo = airmloc(izo, ime, ivert - 1)
                        dxupup = thlayloc(izo, ime, ivert + 1)
                        dxup = thlayloc(izo, ime, ivert)
                        dxdo = thlayloc(izo, ime, ivert - 1)
                        cprime = vanleer_nonunif(cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert - 1)), dtr)
                        dens_bound = (rhoup * dxdo + rhodo * dxup) / (dxup + dxdo)

                        acprime = -atrlo * winvloc(izo, ime, ivert - 1) * dens_bound / thlayloc(izo, ime, ivert)
                        acdo = acprime * dvanleer_nonunif(3, cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert - 1)), dtr)
                        acup = acprime * dvanleer_nonunif(2, cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert - 1)), dtr)
                        acupup = acprime * dvanleer_nonunif(1, cupup, cup, cdo, dxupup, dxup, dxdo, abs(winvloc(izo, ime, ivert - 1)), dtr)
                        aconc(ns, izo, ime, ivert - 1) = aconc(ns, izo, ime, ivert - 1) + &
                                acdo / airmloc(izo, ime, ivert - 1)
                        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                                acup / airmloc(izo, ime, ivert)
                        aconc(ns, izo, ime, ivert + 1) = aconc(ns, izo, ime, ivert + 1) + &
                                acupup / airmloc(izo, ime, ivert + 1)
                    else ! upwind, because only one cell

                        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                                - winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * atrlo
                    endif !of ivert.lt.nverti
                endif ! winvloc
            endif !ivert.gt.1
        else ! inpv
            print*, 'NOT READY NOT READY'
        endif
    endif !inpv

    !  2: Incoming mixing fluxes
    do nv = nverti + 1, 1, -1
        if (vfluxi(izo, ime, ivert, nv).ne.dzero) then
            aconc(ns, izo, ime, nv) = aconc(ns, izo, ime, nv) &
                    + vfluxi(izo, ime, ivert, nv) * atrpr
        endif
    enddo

    !  1: Outgoing mixing fluxes
    aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
            + vfluxo(izo, ime, ivert) * atrlo

END subroutine atransmix
