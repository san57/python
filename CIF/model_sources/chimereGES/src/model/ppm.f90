!  This package of subroutines is the PPM for horizontal transport      


!*******************************************************************************************
subroutine zreconstruct(ns, izo, ime, ivert, y, cfleft, cfright, delc, parab, x)
    !*******************************************************************************************
    !     ******************************************
    !     Parabolic reconstruction (zonal direction)
    !     ******************************************


    ! Attention: only to be used when deltax constant in the direction of ad

    use chimere_consts
    use worker_common

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: y
    real(kind = 8) :: cfleft
    real(kind = 8) :: cfright
    real(kind = 8) :: delc
    real(kind = 8) :: parab
    real(kind = 8) :: x

    ! local variables
    real(kind = 8) :: cc, ccl, ccr, ccrr, ccll

    !*****************************************************************************************

    cc = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)

    !     Boundary

    if (ivert.gt.nverti .or. ((izo.le.1).and.(dom_i==1)) .or. ((izo.ge.nzonal).and.(dom_i==nzdoms))) then

        cfleft = cc
        cfright = cc
        delc = dzero
        parab = dzero
        x = dzero

    else

        !     Courant number
        x = y / xsize(izo, ime)


        !     Interpolations at the faces of the volume

        if ((izo.eq.2).and.(dom_i==1)) then
            ccl = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            ccr = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)
            ccrr = conc(ns, izo + 2, ime, ivert) / airmloc(izo + 2, ime, ivert)

            !           dxl =xsize(1,nm)
            !           dx  =xsize(2,nm)
            !           dxr =xsize(3,nm)
            !           dxrr=xsize(4,nm)

            !           cfleft  = ( ccl*dx + cc*dxl ) / (dxl+dx)
            cfleft = 5.d-1 * (ccl + cc)

            call interp(cfright, ccl, cc, ccr, ccrr)

        elseif ((izo.eq.nzonal - 1).and.(dom_i==nzdoms)) then

            ccll = conc(ns, izo - 2, ime, ivert) / airmloc(izo - 2, ime, ivert)
            ccl = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            ccr = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)

            !           dxll=xsize(nz-2,nm)
            !           dxl =xsize(nz-1,nm)
            !           dx  =xsize(nz,nm)
            !           dxr =xsize(nzonal,nm)

            call interp(cfleft, ccll, ccl, cc, ccr)

            !           cfright = ( cc*dxr + ccr*dx) / (dx+dxr)
            cfright = 5.d-1 * (cc + ccr)

        else

            ccll = conc(ns, izo - 2, ime, ivert) / airmloc(izo - 2, ime, ivert)
            ccl = conc(ns, izo - 1, ime, ivert) / airmloc(izo - 1, ime, ivert)
            ccr = conc(ns, izo + 1, ime, ivert) / airmloc(izo + 1, ime, ivert)
            ccrr = conc(ns, izo + 2, ime, ivert) / airmloc(izo + 2, ime, ivert)

            !           dxll=xsize(nz-2,nm)
            !           dxl =xsize(nz-1,nm)
            !           dx  =xsize(nz  ,nm)
            !           dxr =xsize(nz+1,nm)
            !           dxrr=xsize(nz+2,nm)

            call interp(cfleft, ccll, ccl, cc, ccr)
            call interp(cfright, ccl, cc, ccr, ccrr)

        endif

        !     Second slope limiter
        call filter(cc, cfleft, cfright, parab, delc)

    endif

END subroutine zreconstruct


!*******************************************************************************************
subroutine mreconstruct(ns, izo, ime, ivert, y, cfleft, cfright, delc, parab, x)
    !*******************************************************************************************
    !     *******************************************
    !     Parabolic reconstruction (merid. direction)
    !     *******************************************
    ! Attention: only to be used when deltax constant in the direction of ad

    use chimere_consts
    use worker_common

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: y
    real(kind = 8) :: cfleft
    real(kind = 8) :: cfright
    real(kind = 8) :: delc
    real(kind = 8) :: parab
    real(kind = 8) :: x

    ! local variables
    real(kind = 8) :: cc, ccl, ccr, ccrr, ccll

    !*****************************************************************************************

    cc = conc(ns, izo, ime, ivert) / airmloc(izo, ime, ivert)

    !     Boundary

    if (ivert.gt.nverti .or. ((ime.le.1).and.(dom_j==1)) .or. ((ime.ge.nmerid).and.(dom_j==nmdoms))) then

        cfleft = cc
        cfright = cc
        delc = dzero
        parab = dzero
        x = dzero

    else

        !     Courant number
        x = y / ysize(izo, ime)

        !     Interpolations at the faces of the volume

        if ((ime.eq.2).and.(dom_j==1)) then

            ccl = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            ccr = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            ccrr = conc(ns, izo, ime + 2, ivert) / airmloc(izo, ime + 2, ivert)

            !           dxl =ysize(nz,1)
            !           dx  =ysize(nz,2)
            !           dxr =ysize(nz,3)
            !           dxrr=ysize(nz,4)

            !           cfleft  = ( ccl*dx + cc*dxl ) / (dxl+dx)
            cfleft = 5.d-1 * (ccl + cc)

            call interp(cfright, ccl, cc, ccr, ccrr)

        elseif ((ime.eq.nmerid - 1).and.(dom_j==nmdoms)) then

            ccll = conc(ns, izo, ime - 2, ivert) / airmloc(izo, ime - 2, ivert)
            ccl = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            ccr = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)

            !           dxll=ysize(nz,nm-2)
            !           dxl =ysize(nz,nm-1)
            !           dx  =ysize(nz,nm)
            !           dxr =ysize(nz,nmerid)

            call interp(cfleft, ccll, ccl, cc, ccr)

            !           cfright = ( cc*dxr + ccr*dx) / (dx+dxr)
            cfright = 5.d-1 * (cc + ccr)

        else

            ccll = conc(ns, izo, ime - 2, ivert) / airmloc(izo, ime - 2, ivert)
            ccl = conc(ns, izo, ime - 1, ivert) / airmloc(izo, ime - 1, ivert)
            ccr = conc(ns, izo, ime + 1, ivert) / airmloc(izo, ime + 1, ivert)
            ccrr = conc(ns, izo, ime + 2, ivert) / airmloc(izo, ime + 2, ivert)

            !           dxll=ysize(nz,nm-2)
            !           dxl =ysize(nz,nm-1)
            !           dx  =ysize(nz,nm  )
            !           dxr =ysize(nz,nm+1)
            !           dxrr=ysize(nz,nm+2)

            call interp(cfleft, ccll, ccl, cc, ccr)
            call interp(cfright, ccl, cc, ccr, ccrr)

        endif

        !     Second slope limiter
        call filter(cc, cfleft, cfright, parab, delc)

    endif

END subroutine mreconstruct


!*******************************************************************************************
subroutine delta(cdelta, c1, c2, c3)
    !*******************************************************************************************
    !     *******************
    !     First slope limiter
    !     *******************

    ! Attention: only to be used when deltax constant in the direction of ad

    use chimere_consts
    implicit none

    !*****************************************************************************************
    ! function arguments
    real(kind = 8) :: cdelta
    real(kind = 8) :: c1
    real(kind = 8) :: c2
    real(kind = 8) :: c3

    ! local variables
    real(kind = 8) :: dfront, dbehind

    !*****************************************************************************************

    dfront = c3 - c2
    dbehind = c2 - c1
    cdelta = dzero

    if(dfront>dzero) then
        if(dbehind>dzero) cdelta = 2.d0 * min(2.5d-1 * (c3 - c1), dfront, dbehind)
    else
        if(dbehind<dzero) cdelta = 2.d0 * max(2.5d-1 * (c3 - c1), dfront, dbehind)
    end if

END subroutine delta


!*******************************************************************************************
subroutine interp(cinterp, c1, c2, c3, c4)
    !*******************************************************************************************
    !     ***********************************************
    !     Interpolation at the face between cells 2 and 3
    !     ***********************************************

    ! Attention: only to be used when deltax constant in the direction of ad

    implicit none

    !*****************************************************************************************
    !subroutine arguments
    real(kind = 8) :: cinterp
    real(kind = 8) :: c1
    real(kind = 8) :: c2
    real(kind = 8) :: c3
    real(kind = 8) :: c4

    ! local variables
    real(kind = 8) :: deltal, deltar

    !*****************************************************************************************

    call delta(deltal, c1, c2, c3)
    call delta(deltar, c2, c3, c4)

    cinterp = 5.d-1 * (c3 + c2) + (deltal - deltar) / 6.d0

END subroutine interp


!*******************************************************************************************
subroutine filter(cc, cfleft, cfright, parab, delc)
    !*******************************************************************************************
    !     *************************************************
    !     Second slope limiter for parabolic reconstruction
    !     *************************************************

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    real(kind = 8) :: cc
    real(kind = 8) :: cfleft
    real(kind = 8) :: cfright
    real(kind = 8) :: parab
    real(kind = 8) :: delc

    ! local variables
    real(kind = 8) :: test, dcpar, dc2

    !*****************************************************************************************

    delc = cfright - cfleft
    parab = 6.d0 * (cc - (cfleft + cfright) / 2.d0)
    test = (cfright - cc) * (cc - cfleft)

    if (test.le.0.) then
        cfleft = cc
        cfright = cc
        delc = 0d0
        parab = 0d0
    else
        dcpar = delc * parab
        dc2 = delc * delc
        if (dcpar.gt.dc2) then
            cfleft = 3.d0 * cc - 2.d0 * cfright
            delc = cfright - cfleft
            parab = 6.d0 * (cc - (cfleft + cfright) / 2.d0)
        end if
        if (dcpar.lt.-dc2) then
            cfright = 3.d0 * cc - 2.d0 * cfleft
            delc = cfright - cfleft
            parab = 6.d0 * (cc - (cfleft + cfright) / 2.d0)
        end if
    end if

END subroutine filter


!*******************************************************************************************
subroutine averleft(aver, cfleft, delc, parab, x)
    !*******************************************************************************************
    !     ************************
    !     Left hand side averaging
    !     ************************

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    real(kind = 8) :: aver
    real(kind = 8) :: cfleft
    real(kind = 8) :: delc
    real(kind = 8) :: parab
    real(kind = 8) :: x

    !*****************************************************************************************

    aver = cfleft + x / 2.d0 * (delc + (1.d0 - 2.d0 / 3.d0 * x) * parab)

END subroutine averleft


!*******************************************************************************************
subroutine averright(aver, cfright, delc, parab, x)
    !*******************************************************************************************
    !     *************************
    !     Right hand side averaging
    !     *************************

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    real(kind = 8) :: aver
    real(kind = 8) :: cfright
    real(kind = 8) :: delc
    real(kind = 8) :: parab
    real(kind = 8) :: x

    aver = cfright - x / 2.d0 * (delc - (1.d0 - 2.d0 / 3.d0 * x) * parab)

END subroutine averright


!*******************************************************************************************
subroutine meanconcw(aver, ns, izo, ime, ivert, y)
    !*******************************************************************************************
    !     *************************************************
    !     Mean concentration at the western face during dtr
    !     *************************************************

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    real(kind = 8) :: aver
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: y

    ! local variables
    real(kind = 8) :: cfleft, cfright, delc, parab, x
    !*****************************************************************************************

    !     Parabolic reconstruction
    call zreconstruct(ns, izo, ime, ivert, y, cfleft, cfright, delc, parab, x)

    !     Mean concentration at the western face
    call averleft(aver, cfleft, delc, parab, x)

END subroutine meanconcw


!*******************************************************************************************
subroutine meanconce(aver, ns, izo, ime, ivert, y)
    !*******************************************************************************************
    !     *************************************************
    !     Mean concentration at the eastern face during dtr
    !     *************************************************
    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    real(kind = 8) :: aver
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: y

    ! local variables
    real(kind = 8) :: cfleft, cfright, delc, parab, x
    !*****************************************************************************************

    !     Parabolic reconstruction
    call zreconstruct(ns, izo, ime, ivert, y, cfleft, cfright, delc, parab, x)

    !     Mean concentration at the eastern face
    call averright(aver, cfright, delc, parab, x)

END subroutine meanconce


!*******************************************************************************************
subroutine meanconcs(aver, ns, izo, ime, ivert, y)
    !*******************************************************************************************
    !     **************************************************
    !     Mean concentration at the southern face during dtr
    !     **************************************************

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    real(kind = 8) :: aver
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: y

    ! local variables
    real(kind = 8) :: cfleft, cfright, delc, parab, x
    !*****************************************************************************************

    !     Parabolic reconstruction
    call mreconstruct(ns, izo, ime, ivert, y, cfleft, cfright, delc, parab, x)

    !     Mean concentration at the southern face
    call averleft(aver, cfleft, delc, parab, x)

END subroutine meanconcs


!*******************************************************************************************
subroutine meanconcn(aver, ns, izo, ime, ivert, y)
    !*******************************************************************************************
    !     **************************************************
    !     Mean concentration at the northern face during dtr
    !     **************************************************

    implicit none

    !*****************************************************************************************
    ! subroutine arguments
    real(kind = 8) :: aver
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: y

    ! local variables
    real(kind = 8) :: cfleft, cfright, delc, parab, x
    !*****************************************************************************************

    !     Parabolic reconstruction
    call mreconstruct(ns, izo, ime, ivert, y, cfleft, cfright, delc, parab, x)

    !     Mean concentration at the northern face
    call averright(aver, cfright, delc, parab, x)

END subroutine meanconcn
