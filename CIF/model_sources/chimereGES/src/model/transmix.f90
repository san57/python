subroutine transmix(ns, izo, ime, ivert, trpr, trlo)

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
    real(kind = 8) :: trpr
    real(kind = 8) :: trlo

    ! local variables
    integer :: inp, inpv
    integer :: nv
    real(kind = 8) :: cup, cupup, cdo, cprime
    real(kind = 8) :: concu, concd, flxe
    real(kind = 8) :: rhoup, rhodo, dxupup, dxup, dxdo, dens_bound

    ! external functions
    real(kind = 8) :: vanleer, vanleer_nonunif

    !***************************************************************************
    !  Initializations

    trpr = dzero
    trlo = dzero
    ! lmbb
    inp = species(ns)%transp
    inpv = species(ns)%transpv

    !  Vertical transport
    !  1: Outgoing mixing fluxes

    trlo = trlo + vfluxo(izo, ime, ivert) * conc(ns, izo, ime, ivert)

    !  2: Incoming mixing fluxes
    do nv = 1, nverti + 1
        if (vfluxi(izo, ime, ivert, nv).ne.dzero) then
                trpr = trpr + vfluxi(izo, ime, ivert, nv) * conc(ns, izo, ime, nv)
        endif
    enddo

    !  3: Transport from resolved vertical velocity
    if(inpv.eq.1) then
        ! UPWIND
        !down side
        if(ivert.gt.1)then
            if(winvloc(izo, ime, ivert - 1).gt.dzero)then
                trpr = trpr + winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert - 1)
            endif
            if(winvloc(izo, ime, ivert - 1).lt.dzero)then
                trlo = trlo - winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert)
            endif
        endif
        ! up side
        if(winvloc(izo, ime, ivert).gt.dzero)then
            trlo = trlo + winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert)
        endif
        if(winvloc(izo, ime, ivert).lt.dzero)then
            trpr = trpr - winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert + 1)
        endif
    else
        if(inpv.eq.2) then
            ! VAN LEER
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
                        trpr = trpr + winvloc(izo, ime, ivert - 1) * dens_bound * cprime / thlayloc(izo, ime, ivert)
                    else ! upwind, because only one cell below
                        trpr = trpr + winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert - 1)
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
                        trlo = trlo - winvloc(izo, ime, ivert - 1) * dens_bound * cprime / thlayloc(izo, ime, ivert)
                    else ! upwind, because only one cell
                        trlo = trlo - winvloc(izo, ime, ivert - 1) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert)
                    endif !of ivert.lt.nverti
                endif ! winvloc
            endif !ivert.gt.1
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
                        trpr = trpr - winvloc(izo, ime, ivert) * dens_bound * cprime / thlayloc(izo, ime, ivert)
                    else  ! upwind, because only one cell
                        trpr = trpr - winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert + 1)
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
                        trlo = trlo + winvloc(izo, ime, ivert) * dens_bound * cprime / thlayloc(izo, ime, ivert)
                    else
                        trlo = trlo + winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert)
                    endif !ivert.gt.1
                endif ! winvloc(izo,ime,ivert).lt.0
            else ! ivert.lt.nverti
                !upwind, because only one cell
                if(winvloc(izo, ime, ivert).gt.dzero) then
                    trlo = trlo + winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert)
                endif
                if(winvloc(izo, ime, ivert).lt.dzero) then
                    trpr = trpr - winvloc(izo, ime, ivert) / thlayloc(izo, ime, ivert) * conc(ns, izo, ime, ivert + 1)
                endif
            endif !end treatment of upper side
        else ! inpv
            print*, 'NOT AN EXISTING SCHEME'
        endif
    endif !inpv

    ! DEEP CONVECTION IMPLEMENTATION
    if(ideepconv.ne.0.and.ideep(izo, ime).eq.1)then

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

        ! ENTRAINMENT FLUX
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

        ! SOLVING TRANSMIX PRODUCTION AND LOSS
        if(ivert.eq.1)then
            trpr = trpr - flxdconc(ns, izo, ime, ivert) / thlayloc(izo, ime, ivert)
            trlo = trlo + flxuconc(ns, izo, ime, ivert) / thlayloc(izo, ime, ivert)
            if(flxeconc(ns, izo, ime, ivert).lt.dzero)then
                trpr = trpr - flxeconc(ns, izo, ime, ivert) / thlayloc(izo, ime, ivert)
            else
                trlo = trlo + flxeconc(ns, izo, ime, ivert) / thlayloc(izo, ime, ivert)
            endif
        else
            if((flxuconc(ns, izo, ime, ivert) - flxuconc(ns, izo, ime, ivert - 1)).lt.dzero)then

                trpr = trpr - (flxuconc(ns, izo, ime, ivert) - flxuconc(ns, izo, ime, ivert - 1)) &
                        / thlayloc(izo, ime, ivert)
            else
                trlo = trlo + (flxuconc(ns, izo, ime, ivert) - flxuconc(ns, izo, ime, ivert - 1)) &
                        / thlayloc(izo, ime, ivert)
            endif
            if((flxdconc(ns, izo, ime, ivert) - flxdconc(ns, izo, ime, ivert - 1)).lt.dzero)then
                trpr = trpr - (flxdconc(ns, izo, ime, ivert) - flxdconc(ns, izo, ime, ivert - 1)) &
                        / thlayloc(izo, ime, ivert)
            else
                trlo = trlo + (flxdconc(ns, izo, ime, ivert) - flxdconc(ns, izo, ime, ivert - 1)) &
                        / thlayloc(izo, ime, ivert)
            endif
            if((flxeconc(ns, izo, ime, ivert) - flxeconc(ns, izo, ime, ivert - 1)).lt.dzero)then
                trpr = trpr - (flxeconc(ns, izo, ime, ivert) - flxeconc(ns, izo, ime, ivert - 1)) &
                        / thlayloc(izo, ime, ivert)
            else
                trlo = trlo + (flxeconc(ns, izo, ime, ivert) - flxeconc(ns, izo, ime, ivert - 1)) &
                        / thlayloc(izo, ime, ivert)
            endif
        endif
    endif ! of ideepconv=1
    ! DEEP CONVECTION IMPLEMENTATION

    !  Horizontal transport
    !  1: Western side
    if(uwest(izo, ime, ivert).gt.dzero) then
        if(inp.eq.1) then
            ! PPM
            call meanconce(cprime, ns, izo - 1, ime, ivert, uwest(izo, ime, ivert) * dtr)
            trpr = trpr + fluxw(izo, ime, ivert) * cprime * airmloc(izo - 1, ime, ivert)
        elseif(inp.eq.2) then
            ! Vanleer
            cupup = conc(ns, izo - 2, ime, ivert) / airmloc(izo - 2, ime, ivert)
            cup = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,xsize(izo-1,ime),uwest(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo - 2, ime), xsize(izo - 1, ime), xsize(izo, ime), uwest(izo, ime, ivert), dtr)
            trpr = trpr + fluxw(izo, ime, ivert) * cprime * airmloc(izo - 1, ime, ivert)
        else
            ! UPWIND
            trpr = trpr + fluxw(izo, ime, ivert) * conc(ns, izo - 1, ime, ivert)
        endif
    else
        if(inp.eq.1) then
            ! PPM
            call meanconcw(cprime, ns, izo, ime, ivert, -uwest(izo, ime, ivert) * dtr)
            trlo = trlo - fluxw(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        elseif(inp.eq.2) then
            cupup = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,xsize(izo,ime),-uwest(izo,ime,ivert),dtr) 
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo + 1, ime), xsize(izo, ime), xsize(izo - 1, ime), -uwest(izo, ime, ivert), dtr)
            trlo = trlo - fluxw(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        else
            ! UPWIND
            trlo = trlo - fluxw(izo, ime, ivert) * conc(ns, izo, ime, ivert)
        endif
    endif
     !  2: Eastern side
    if(ueast(izo, ime, ivert).lt.dzero) then
        if(inp.eq.1) then
            call meanconcw(cprime, ns, izo + 1, ime, ivert, -ueast(izo, ime, ivert) * dtr)
            trpr = trpr - fluxe(izo, ime, ivert) * cprime * airmloc(izo + 1, ime, ivert)
        elseif(inp.eq.2) then
            cupup = conc(ns, izo + 2, ime, ivert) / airmloc(izo + 2, ime, ivert)
            cup = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,xsize(izo+1,ime),-ueast(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo + 2, ime), xsize(izo + 1, ime), xsize(izo, ime), -ueast(izo, ime, ivert), dtr)
            trpr = trpr - fluxe(izo, ime, ivert) * cprime * airmloc(izo + 1, ime, ivert)
        else
            ! UPWIND
            trpr = trpr - fluxe(izo, ime, ivert) * conc(ns, izo + 1, ime, ivert)
        endif
    else
        if(inp.eq.1) then
            call meanconce(cprime, ns, izo, ime, ivert, ueast(izo, ime, ivert) * dtr)
            trlo = trlo + fluxe(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        elseif(inp.eq.2) then
            cupup = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)

            !cprime = vanleer(cupup,cup,cdo,xsize(izo,ime),ueast(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, xsize(izo - 1, ime), xsize(izo, ime), xsize(izo + 1, ime), ueast(izo, ime, ivert), dtr)
            trlo = trlo + fluxe(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        else
            ! UPWIND
            trlo = trlo + fluxe(izo, ime, ivert) * conc(ns, izo, ime, ivert)
        endif
    endif
    !  3: Southern side
    if(usouth(izo, ime, ivert).gt.dzero) then
        if(inp.eq.1) then
            ! PPM
            call meanconcn(cprime, ns, izo, ime - 1, ivert, usouth(izo, ime, ivert) * dtr)
            trpr = trpr + fluxs(izo, ime, ivert) * cprime * airmloc(izo, ime - 1, ivert)
        elseif(inp.eq.2) then
            cupup = conc(ns, izo, ime - 2, ivert) / airmloc(izo, ime - 2, ivert)
            cup = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime = vanleer(cupup,cup,cdo,ysize(izo,ime-1),usouth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime - 2), ysize(izo, ime - 1), ysize(izo, ime), usouth(izo, ime, ivert), dtr)
            trpr = trpr + fluxs(izo, ime, ivert) * cprime * airmloc(izo, ime - 1, ivert)
        else
            ! UPWIND
            trpr = trpr + fluxs(izo, ime, ivert) * conc(ns, izo, ime - 1, ivert)
        endif
    else
        if(inp.eq.1) then
            ! PPM
            call meanconcs(cprime, ns, izo, ime, ivert, -usouth(izo, ime, ivert) * dtr)
            trlo = trlo - fluxs(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        elseif(inp.eq.2) then
            cupup = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            !cprime = vanleer(cupup,cup,cdo,ysize(izo,ime),-usouth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime + 1), ysize(izo, ime), ysize(izo, ime - 1), -usouth(izo, ime, ivert), dtr)
            trlo = trlo - fluxs(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        else
            ! UPWIND
            trlo = trlo - fluxs(izo, ime, ivert) * conc(ns, izo, ime, ivert)
        endif
    endif
     !  4: Northern side

    if(unorth(izo, ime, ivert).lt.dzero) then
        if(inp.eq.1) then
            call meanconcs(cprime, ns, izo, ime + 1, ivert, -unorth(izo, ime, ivert) * dtr)
            trpr = trpr - fluxn(izo, ime, ivert) * cprime * airmloc(izo, ime + 1, ivert)
        elseif(inp.eq.2) then
            cupup = conc(ns, izo, ime + 2, ivert) / airmloc(izo, ime + 2, ivert)
            cup = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            cdo = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            !cprime =vanleer(cupup,cup,cdo,ysize(izo,ime+1),-unorth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime + 2), ysize(izo, ime + 1), ysize(izo, ime), -unorth(izo, ime, ivert), dtr)
            trpr = trpr - fluxn(izo, ime, ivert) * cprime * airmloc(izo, ime + 1, ivert)
        else
            ! UPWIND
            trpr = trpr - fluxn(izo, ime, ivert) * conc(ns, izo, ime + 1, ivert)
        endif
    else
        if(inp.eq.1) then
            call meanconcn(cprime, ns, izo, ime, ivert, unorth(izo, ime, ivert) * dtr)
            trlo = trlo + fluxn(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        elseif(inp.eq.2) then
            cupup = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            cup = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)
            cdo = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            !cprime = vanleer(cupup,cup,cdo,ysize(izo,ime),unorth(izo,ime,ivert),dtr)
            cprime = vanleer_nonunif(cupup, cup, cdo, ysize(izo, ime - 1), ysize(izo, ime), ysize(izo, ime + 1), unorth(izo, ime, ivert), dtr)
            trlo = trlo + fluxn(izo, ime, ivert) * cprime * airmloc(izo, ime, ivert)
        else
            ! UPWIND
            trlo = trlo + fluxn(izo, ime, ivert) * conc(ns, izo, ime, ivert)
        endif
    endif
!  if(ime==66.and.izo==57.and.ns==32) print 102,'TTTTTTTTTTTTTT ',trpr


END subroutine transmix
