function vanleer(cupup, cup, cdo, dx, u, dt)

    implicit none

    !****************************************************************
    ! function arguments
    real(kind = 8) :: cupup
    real(kind = 8) :: cup
    real(kind = 8) :: cdo
    real(kind = 8) :: dx
    real(kind = 8) :: u
    real(kind = 8) :: dt

    ! local variables
    real(kind = 8) :: cfl, vanleer, minmod

    !*****************************************************************
    cfl = u * dt / dx
    if(cfl.gt.1.) then
        print *, '*** WARNING: CFL > 1 IN VANLEER scheme :', cfl
    endif

    !  Nonmonotonic cases

    if(cup.lt.min(cupup, cdo)) then
        vanleer = cup
    elseif(cup.gt.max(cupup, cdo)) then
        vanleer = cup
    else

        !  Monotonic case

        if(abs(cdo - cup).lt.abs(cup - cupup)) then
            minmod = cdo - cup
        else
            minmod = cup - cupup
        endif
        vanleer = cup + 0.5 * max(1d0 - cfl, 0d0) * minmod
    endif

END function vanleer


function vanleer_nonunif(cupup, cup, cdo, dxupup, dx, dxdo, u, dt)

    implicit none

    !****************************************************************
    ! function arguments
    real(kind = 8) :: cupup
    real(kind = 8) :: cup
    real(kind = 8) :: cdo
    real(kind = 8) :: dx, dxupup, dxdo
    real(kind = 8) :: u
    real(kind = 8) :: dt

    ! local variables
    real(kind = 8) :: cfl, vanleer_nonunif, minmod

    !*****************************************************************

    !  Nonmonotonic cases

    if(cup.lt.min(cupup, cdo)) then
        vanleer_nonunif = cup
    elseif(cup.gt.max(cupup, cdo)) then
        vanleer_nonunif = cup
    else

        !  Monotonic case

        if(abs((cdo - cup) / (dx + dxdo)).lt.abs((cup - cupup) / (dx + dxupup))) then
            cfl = 2. * u * dt / (dxdo + dx)
            if(cfl.gt.1.) then
                !!!!!!print *, '*** WARNING: CFL > 1 IN VANLEER scheme :', cfl
                vanleer_nonunif = cup
            else
                vanleer_nonunif = (cup * dxdo + cdo * dx) / (dx + dxdo) - cfl / 2. * (cdo - cup)
            endif
        else
            cfl = 2. * u * dt / (dxupup + dx)
            if(cfl.gt.1.) then
                !!!!!print *, '*** WARNING: CFL > 1 IN VANLEER scheme :', cfl
                vanleer_nonunif = cup
            else
                vanleer_nonunif = cup + dx * (cup - cupup) / (dx + dxupup) - cfl / 2. * (cup - cupup)
            endif
        endif
    endif

END function vanleer_nonunif
