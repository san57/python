subroutine outprint_tl(iprint, print_title, nh, writeout, parout)
    
    !  Outputs on file and screen.
    
    use netcdf
    use chimere_consts
    use chimere_common
    use wholedomain_common
    use master_message_subs
    
    implicit none
    
#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'outprint.f90')

    !****************************************************************************
    integer :: iprint
    logical :: print_title
    integer :: nh
    integer :: writeout, parout
    
    integer :: nf, ne, no, ns
    integer :: nsp
    integer :: cl
    
    integer :: ncstat
    integer :: ioutspec
    integer :: iout, nout
    integer :: ispec, izo, ime, ivert
    integer :: tsaveconcs
    
    real(kind = 8) :: fppb, fdry
    real(kind = 4), allocatable, dimension(:, :, :) :: toprint
    real(kind = 4), allocatable, dimension(:, :, :) :: toprint_tl
    real(kind = 4), allocatable, dimension(:, :) :: toprint2d
    real(kind = 4), allocatable, dimension(:, :, :) :: toprint3d
    real(kind = 4), allocatable, dimension(:, :, :, :, :) :: toprintexc
    real, allocatable, dimension(:, :, :) :: toscreen
    real(kind = 8), allocatable, dimension(:, :, :) :: buf3d, buf3do, buf3d_tl, buf3do_tl
    
    character(len = dlen) :: datebuf
    
    ! Functions
    character(len = dlen) :: numeric2mm5date
    !***************************************************************************
    allocate(toprint(nzonal_domain, nmerid_domain, nverti))
    allocate(toprint_tl(nzonal_domain, nmerid_domain, nverti))
    allocate(toscreen(2, nzonal_domain, nmerid_domain))
    allocate(buf3d(nzonal_domain, nmerid_domain, nverti))
    allocate(buf3do(nzonal_domain, nmerid_domain, nverti))
    allocate(buf3d_tl(nzonal_domain, nmerid_domain, nverti))
    allocate(buf3do_tl(nzonal_domain, nmerid_domain, nverti))
    
    !  print*,'  Output to files "out"'
    
    ! TIMES
    datebuf = numeric2mm5date(idate(ihourrun))
    write(unit = datebuf(15:16), fmt = '(i2.2)')  mod(60 * iprint, 60)
    if(writeout.eq.1)then
        ncstat = nf90_put_var(out_ncid, out_times_varid, datebuf, (/1, iprint + 1/), (/dlen, 1/))
        NCERR(__LINE__)
    endif
    
    
    ! SPECIES
    do ioutspec = 1, noutspec
        call master_recv_toprint(toprint)
        if(writeout.eq.1)then
            ncstat = nf90_put_var(&
                    out_ncid, &
                    output_species(ioutspec)%varid, &
                    toprint(:, :, :), &
                    (/1, 1, 1, iprint + 1/), (/nzonal_domain, nmerid_domain, nivout, 1/) &
                    )
            NCERR(__LINE__)
        endif
        if (ioutspec==1) toscreen(1, :, :) = toprint(:, :, 1)
        call master_recv_toprint(toprint_tl)
        if(writeout.eq.1)then
            ncstat = nf90_put_var(&
                    out_ncid, &
                    output_species(ioutspec)%varid_tl, &
                    toprint_tl(:, :, :), &
                    (/1, 1, 1, iprint + 1/), (/nzonal_domain, nmerid_domain, nivout, 1/) &
                    )
            NCERR(__LINE__)
        endif
        if (ioutspec==1) toscreen(2, :, :) = toprint_tl(:, :, 1)
        if (ioutspec==2) toscreen(2, :, :) = toprint(:, :, 1)
    end do ! ioutspec=1,noutspec
    
   
    ! Synchronize disk to avoid data loss
    if(writeout.eq.1)then
        ncstat = nf90_sync(out_ncid)
        NCERR(__LINE__)
    endif
    
    ! Output of other gridded variables to file "par"
    
    ! TIMES
    if(parout.eq.1)then
        ncstat = nf90_put_var(par_ncid, par_times_varid, datebuf, (/1, iprint + 1/), (/dlen, 1/))
        NCERR(__LINE__)
    endif
    
    
    ! METEO FIELDS
    if(parout.eq.1)then
        allocate(toprint2d(nzonal_domain, nmerid_domain))
        do iout = 1, nparammax
            if(output_par(iout)%printit) then
                do ime = 1, nmerid_domain
                    do izo = 1, nzonal_domain
                        select case (output_par(iout)%name)
                        case('winhs')
                            toprint2d(izo, ime) = &
                                    1.e-2 * sqrt(winzloc(izo, ime, 1)**2 + winmloc(izo, ime, 1)**2)
                        case('winh2')
                            toprint2d(izo, ime) = &
                                    1.e-2 * sqrt(winzloc(izo, ime, 2)**2 + winmloc(izo, ime, 2)**2)
                        case('winvs')
                            toprint2d(izo, ime) = 1.e-2 * winvloc(izo, ime, 1)
                        case('winvt')
                            toprint2d(izo, ime) = 1.e-2 * winvloc(izo, ime, nverti)
                        case('sphu')
                            toprint2d(izo, ime) = sphuloc(izo, ime, 1) / (1.6 * airmloc(izo, ime, 1))
                        case('airm')
                            toprint2d(izo, ime) = airmloc(izo, ime, 1)
                        case('temp')
                            toprint2d(izo, ime) = temploc(izo, ime, 1) - t0k
                        case('winxs')
                            toprint2d(izo, ime) = 1.e-2 * winxloc(izo, ime, 1)
                        case('winxt')
                            toprint2d(izo, ime) = 1.e-2 * winxloc(izo, ime, nverti)
                        case('thlay')
                            toprint2d(izo, ime) = 1.e-2 * thlayloc(izo, ime, 1)
                        case('hght')
                            toprint2d(izo, ime) = 1.e-2 * hghtloc(izo, ime)
                        case('atte')
                            toprint2d(izo, ime) = atteloc(izo, ime)
                        case('zeni')
                            toprint2d(izo, ime) = zeniloc(izo, ime)
                        case('tem2')
                            toprint2d(izo, ime) = tem2loc(izo, ime) - t0k
                        case('usta')
                            toprint2d(izo, ime) = 1.e-2 * ustaloc(izo, ime)
                        case('aerr')
                            toprint2d(izo, ime) = 1.e+2 * aerrloc(izo, ime)
                        case('obuk')
                            toprint2d(izo, ime) = 1.e-2 * obukloc(izo, ime)
                        case('wsta')
                            toprint2d(izo, ime) = 1.e-2 * wstaloc(izo, ime)
                        case('depo')
                            toprint2d(izo, ime) = depoloc(1, izo, ime) * hlayloc(izo, ime, 1)
                        case('topc')
                            toprint2d(izo, ime) = topcloc(izo, ime)
                        case('clwc')
                            toprint2d(izo, ime) = clwcloc(izo, ime, 1) * an / (airmloc(izo, ime, 1) * 29.)
                        case default
                            call exit1('outprint.f90 : error in output parameter list')
                        end select
                    end do !izo=1,nzonal_domain
                end do !ime=1,nmerid_domain
                ncstat = nf90_put_var(&
                        par_ncid, &
                        output_par(iout)%varid, &
                        toprint2d, &
                        (/1, 1, iprint + 1/), (/nzonal_domain, nmerid_domain, 1/)          &
                        )
                NCERR(__LINE__)
            end if !(output_par(iout)%printit)
        enddo ! iout=1,nparammax
        deallocate (toprint2d)
    endif
    
    ! Synchronize disk to avoid data loss
    if(parout.eq.1)then
        ncstat = nf90_sync(par_ncid)
        NCERR(__LINE__)
    endif
    
    !print*,'  Output on screen'
    
    if (print_title) then
        
        print 105, 'DATE        '                                   &
                , output_species(1)%name(1:cl(output_species(1)%name)) &
                , output_species(1)%name(1:cl(output_species(1)%name)) // '_tl' &
                , ' Temp (2m)'                                         &
                , ' Surf Wind'                                         &
                , '  PBL Hght'                                         &
                , '  Rad Atte'                                         &
                , '   Kz(1)  ', '   Kz(2)  '
    endif
    
    izo = int(nzonal_domain / 2.) + 1
    ime = int(nmerid_domain / 2.) + 1
    print 104, idate(ihourrun), &
            toscreen(1, izo, ime) &
            !       ,toscreen(2,izo,ime) &
            , tem2loc(izo, ime) - t0k                                                    &
            , sqrt(winzloc(izo, ime, 1)**2 + winmloc(izo, ime, 1)**2) * 1e-2           &
            , hghtloc(izo, ime) * 1e-2                                                   &
            , atteloc(izo, ime)                                                        &
            , kzzzloc(izo, ime, 1) * 1e-4                                            &
            , kzzzloc(izo, ime, 2) * 1e-4
    !104 format(i10,f8.1,a5,f8.1,f10.1,f10.1,f10.0,f10.2,f10.2,f10.2)
    104 format(i10, f8.1, f8.1, f10.1, f10.1, f10.0, f10.2, f10.2, f10.2)
    !105 format(a9,a8,a9,3x,6(a10))
    105 format(a9, a8, 3x, 6(a10))
    if(mod(nh, nsaveconcs).eq.0) then
        ! TIMES
        datebuf = numeric2mm5date(idate(ihourrun))
        tsaveconcs = ihourrun / nsaveconcs
        ncstat = nf90_put_var(end_ncid, end_times_varid, datebuf, (/1, tsaveconcs/), (/dlen, 1/))
        NCERR(__LINE__)
        
        !print*,' CONC for "end" file'
        do ispec = 1, nspectot
            call master_recv_conc
            call master_recv_conc_tl
            do ivert = 1, nverti
                do ime = 1, nmerid_domain
                    do izo = 1, nzonal_domain
                        buf3d(izo, ime, ivert) = conc(izo, ime, ivert)
                        buf3do(izo, ime, ivert) = conco(izo, ime, ivert)
                        buf3d_tl(izo, ime, ivert) = conc_tl(izo, ime, ivert)
                        buf3do_tl(izo, ime, ivert) = conco_tl(izo, ime, ivert)
                    end do
                end do
            end do
            ncstat = nf90_put_var(&
                    end_ncid, &
                    species(ispec)%varid, &
                    buf3d, &
                    (/1, 1, 1, tsaveconcs/), (/nzonal_domain, nmerid_domain, nverti, 1/) &
                    )
            NCERR(__LINE__)
            ncstat = nf90_put_var(&
                    end_ncid, &
                    species(ispec)%varido, &
                    buf3do, &
                    (/1, 1, 1, tsaveconcs/), (/nzonal_domain, nmerid_domain, nverti, 1/) &
                    )
            NCERR(__LINE__)
            ncstat = nf90_put_var(&
                    end_ncid, &
                    species(ispec)%varid_tl, &
                    buf3d_tl, &
                    (/1, 1, 1, tsaveconcs/), (/nzonal_domain, nmerid_domain, nverti, 1/) &
                    )
            NCERR(__LINE__)
            ncstat = nf90_put_var(&
                    end_ncid, &
                    species(ispec)%varido_tl, &
                    buf3do_tl, &
                    (/1, 1, 1, tsaveconcs/), (/nzonal_domain, nmerid_domain, nverti, 1/) &
                    )
            NCERR(__LINE__)
        
        end do
        ! Synchronize disk to avoid data loss
        ncstat = nf90_sync(end_ncid)
        NCERR(__LINE__)
    end if
    
    deallocate(toprint)
    deallocate(toprint_tl)
    deallocate(toscreen)
    deallocate(buf3d)
    deallocate(buf3do)
    deallocate(buf3d_tl)
    deallocate(buf3do_tl)

end subroutine outprint_tl
