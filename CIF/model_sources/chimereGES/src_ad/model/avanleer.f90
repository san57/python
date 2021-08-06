function dvanleer(var, cupup, cup, cdo, dx, u, dt)

    implicit none

    !****************************************************************
    ! function arguments
    real(kind = 8) :: cupup
    real(kind = 8) :: cup
    real(kind = 8) :: cdo
    real(kind = 8) :: dx
    real(kind = 8) :: u
    real(kind = 8) :: dt
    integer :: var

    ! local variables
    real(kind = 8) :: cfl, dvanleer

    !*****************************************************************
    cfl = u * dt / dx
    !  if(cfl.gt.1.) then
    !     print *,'*** WARNING: CFL > 1 IN VANLEER scheme :',cfl
    !  endif

    !  Nonmonotonic cases

    if(cup.lt.min(cupup, cdo)) then

        ! vanleer = cup then the derivative against other variables = 0
        if (var.eq.2) then
            dvanleer = 1.
        else
            dvanleer = 0.
        endif

    elseif(cup.gt.max(cupup, cdo)) then

        if(var.eq.2)then
            dvanleer = 1.
        else
            dvanleer = 0.
        endif

    else

        !  Monotonic case

        if(abs(cdo - cup).lt.abs(cup - cupup)) then
            ! vanleer = cup + 0.5*(1-cfl)*(cdo - cup)
            if(var.eq.1)then
                dvanleer = 0.
            elseif(var.eq.2) then
                dvanleer = 1. - 0.5 * (1. - cfl)
            elseif(var.eq.3) then
                dvanleer = 0.5 * (1. - cfl)
            endif
        else
            ! vanleer = cup + 0.5*(1-cfl)*(cup - cupup)
            if(var.eq.1) then
                dvanleer = -0.5 * (1. - cfl)
            elseif(var.eq.2) then
                dvanleer = 1. + 0.5 * (1. - cfl)
            elseif(var.eq.3) then
                dvanleer = 0.
            endif
        endif
    endif

END function dvanleer

function dvanleer_nonunif(var, cupup, cup, cdo, dxupup, dx, dxdo, u, dt)

    implicit none

    !****************************************************************
    ! function arguments
    real(kind = 8) :: cupup
    real(kind = 8) :: cup
    real(kind = 8) :: cdo
    real(kind = 8) :: dx, dxupup, dxdo
    real(kind = 8) :: u
    real(kind = 8) :: dt
    integer :: var

    ! local variables
    real(kind = 8) :: cfl, dvanleer_nonunif, minmod

    !*****************************************************************

    !  Nonmonotonic cases

    if(cup.lt.min(cupup, cdo)) then
        !   vanleer_nonunif = cup then the derivative against other variables = 0
        if (var.eq.2) then
            dvanleer_nonunif = 1.
        else
            dvanleer_nonunif = 0.
        endif

    elseif(cup.gt.max(cupup, cdo)) then
        !    vanleer_nonunif = cup
        if (var.eq.2) then
            dvanleer_nonunif = 1.
        else
            dvanleer_nonunif = 0.
        endif

    else

        !  Monotonic case

        if(abs((cdo - cup) / (dx + dxdo)).lt.abs((cup - cupup) / (dx + dxupup))) then
            cfl = 2. * u * dt / (dxdo + dx)
            if(cfl.gt.1.) then
                !       vanleer_nonunif=cup
                if (var.eq.2) then
                    dvanleer_nonunif = 1.
                else
                    dvanleer_nonunif = 0.
                endif

            else
                !        vanleer_nonunif=(cup*dxdo+cdo*dx)/(dx+dxdo)-cfl/2.*(cdo-cup)
                if(var.eq.1)then
                    dvanleer_nonunif = 0.
                elseif(var.eq.2) then
                    dvanleer_nonunif = dxdo / (dx + dxdo) + cfl / 2.
                elseif(var.eq.3)then
                    dvanleer_nonunif = dx / (dx + dxdo) - cfl / 2.
                endif

            endif
        else
            cfl = 2. * u * dt / (dxupup + dx)
            if(cfl.gt.1.) then
                !       vanleer_nonunif=cup
                if (var.eq.2) then
                    dvanleer_nonunif = 1.
                else
                    dvanleer_nonunif = 0.
                endif

            else
                !        vanleer_nonunif=cup+dx*(cup-cupup)/(dx+dxupup)-cfl/2.*(cup-cupup)
                if(var.eq.1)then
                    dvanleer_nonunif = -dx * 1. / (dx + dxupup) + cfl / 2.
                elseif(var.eq.2) then
                    dvanleer_nonunif = 1. + dx / (dx + dxupup) - cfl / 2.
                elseif(var.eq.3)then
                    dvanleer_nonunif = 0.
                endif

            endif
        endif
    endif

END function dvanleer_nonunif

