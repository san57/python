subroutine iniboun
    
    !  Initialization of boundary concentrations
    
    use netcdf
    use chimere_consts
    use chimere_common
    use wholedomain_common
    
    implicit none
    
#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'iniboun.f90')

    !****************************************************************************
    integer :: ns, ni, nl, ih, ihz
    integer :: ifnbounspec
    integer :: idr, idex
    integer :: izo, ime, ivert
    integer :: sn, we
    integer :: ntimes, hb, tb, bt, dl, sp
    character(len = dlen) :: datebuf
    character(len = splen), dimension(:), allocatable :: specbuf
    real(kind = 8), allocatable, dimension(:, :, :) :: buf3
    
    ! netCDF stuff
    integer :: ncstat                    ! return code for netCDF functions
    integer :: timedimid                 ! time dimension ID of conc files
    integer :: datedimid                 ! date string dimension ID of conc files
    integer :: hbdimid                   ! horizontal boundary dimension ID
    integer :: btdimid                   ! bottom_top dimension ID of bounconc file
    integer :: wedimid                   ! west_east dimension ID
    integer :: sndimid                   ! south_north dimension ID
    integer :: spdimid                   ! species dimnsion ID of bounconc file
    integer :: bounconc_species_varid
    
    integer, dimension(3) :: stvec3  ! start vectors for R/W functions
    integer, dimension(3) :: cntvec3 ! count vectors for R/W functions
    integer, dimension(4) :: stvec4  ! start vectors for R/W functions
    integer, dimension(4) :: cntvec4 ! count vectors for R/W functions
    
    
    ! External function
    integer :: mm5date2numeric
    
    character*15 cspec
    character*100000 usedval, notusedval
    
    !*****************************************************************************
    
    !!  List of lateral boundary species
    ! how many species ?
    ncstat = nf90_inq_dimid(bounconc_ncid, 'Species', spdimid)
NCERR(__LINE__)
    ncstat = nf90_inquire_dimension(bounconc_ncid, spdimid, len = nspecboun)
NCERR(__LINE__)
    
    
    allocate(specbuf(nspecboun))
    allocate(isboun(max(nspec, nspecboun)))
    
    ! read "species" variable
    ncstat = nf90_inq_varid(bounconc_ncid, 'species', bounconc_species_varid)
NCERR(__LINE__)
    ncstat = nf90_get_var(&
            bounconc_ncid, bounconc_species_varid, specbuf)
NCERR(__LINE__)
    
    do ns = 1, max(nspec, nspecboun)
        isboun(ns) = 0
    enddo
    usedval = ''
    notusedval = ''
    
    print*, '* BOUNDARY SPECIES: '
    nspecbounout = 0
    do ni = 1, nspecboun
        call findspec(specbuf(ni), ns)
        if((ns.gt.0).and.(ns.le.nspec)) then
            nspecbounout = nspecbounout + 1
            isboun(ni) = ns
            cspec = specbuf(ni)
            usedval = usedval(1:len_trim(usedval))  &
                    & // specbuf(ni)(1:len_trim(specbuf(ni))) // '; '
        else
            notusedval = notusedval(1:len_trim(notusedval))  &
                    & // specbuf(ni)(1:len_trim(specbuf(ni))) // '; '
        endif
    enddo
    deallocate(specbuf)
    print*, '  o used species:'
    if(len_trim(usedval).gt.1)    print*, 'u   ', usedval(1:len_trim(usedval) - 1)
    print*, '  o NOT used species:'
    if(len_trim(notusedval).gt.1) print*, ' nu   ', notusedval(1:len_trim(notusedval) - 1)
    
    ! Boundary concentrations file is already open
    
    ! Checking dimensions consistency between file and chimere parameters
    !   Getting dimensions IDs
    !     Are Boundary concentrations time-dependant ?
    ncstat = nf90_inq_dimid(bounconc_ncid, 'Time', timedimid)
    ncstat = nf90_inq_dimid(bounconc_ncid, 'DateStrLen', datedimid)
NCERR(__LINE__)
    ncstat = nf90_inq_dimid(bounconc_ncid, 'h_boundary', hbdimid)
NCERR(__LINE__)
    ncstat = nf90_inq_dimid(bounconc_ncid, 'bottom_top', btdimid)
NCERR(__LINE__)
    ncstat = nf90_inq_dimid(bounconc_ncid, 'Species', spdimid)
NCERR(__LINE__)
    ncstat = nf90_inq_dimid(bounconc_ncid, 'west_east', wedimid)
NCERR(__LINE__)
    ncstat = nf90_inq_dimid(bounconc_ncid, 'south_north', sndimid)
NCERR(__LINE__)
    
    !   Checking lengths
    ncstat = nf90_inquire_dimension(bounconc_ncid, timedimid, len = ntimes)
NCERR(__LINE__)
    ncstat = nf90_inquire_dimension(bounconc_ncid, datedimid, len = dl)
NCERR(__LINE__)
    ncstat = nf90_inquire_dimension(bounconc_ncid, hbdimid, len = hb)
NCERR(__LINE__)
    ncstat = nf90_inquire_dimension(bounconc_ncid, btdimid, len = bt)
NCERR(__LINE__)
    ncstat = nf90_inquire_dimension(bounconc_ncid, wedimid, len = we)
NCERR(__LINE__)
    ncstat = nf90_inquire_dimension(bounconc_ncid, sndimid, len = sn)
NCERR(__LINE__)
    if (ntimes<nhourrun) then
        print *, 'iniboun : Boundary concentrations file '
        print *, 'does not contain enough time slots. Exiting.'
        call exit1('Exiting')
    end if
    if (dl /= dlen)        call exit1('iniboun : date format error in Boundary concentrations file')
    if (hb /= nhbound_domain)     call exit1(&
            'iniboun : h_boundary dimension inconsistency in Boundary concentrations file')
    if (bt /= nverti)      call exit1(&
            'iniboun : bottom_top dimension inconsistency in Boundary concentrations file')
    if (we /= nzonal_domain)      call exit1(&
            'iniboun : west_east dimension inconsistency in Boundary concentrations file')
    if (sn /= nmerid_domain)      call exit1(&
            'iniboun : west_east dimension inconsistency in Boundary concentrations file')
    
    ! Checking time
    ncstat = nf90_inq_varid(bounconc_ncid, 'Times', bounconc_times_varid)
NCERR(__LINE__)
    ncstat = nf90_get_var(&
            bounconc_ncid, bounconc_times_varid, datebuf, (/1, 1/), (/dlen, 1/))
NCERR(__LINE__)
    idex = idate(0)
    idr = mm5date2numeric(datebuf)
    if(idr.ne.idex) then
        print *, '*** ERROR: WRONG EXPECTED DATE IN BOUN_CONCS FILE'
        print *, 'IHOURRUN=0,  EXPECTED= ', idex, ' BOUN_CONCS= ', idr
        call exit1('Exiting')
    endif
    
    !*** Lateral Boundary Concentrations
    
    ! Initialization
    do ns = 1, nspec
        do nl = 1, nlatbound_domain
            boundlat(nl, ns, 1) = 1d-17
            boundlat(nl, ns, 2) = 1d-17
        enddo
    enddo
    
    !' Reading lateral concentrations in ppb'
    ncstat = nf90_inq_varid(bounconc_ncid, 'lat_conc', latconc_conc_varid)
NCERR(__LINE__)
    stvec4 = (/1, 1, 1, 1/)
    cntvec4 = (/nspecboun, nhbound_domain, nverti, 1/)
    allocate(buf3(nspecboun, nhbound_domain, nverti))
    ncstat = nf90_get_var(bounconc_ncid, latconc_conc_varid, buf3, stvec4, cntvec4)
NCERR(__LINE__)
    
    nl = 0
    do ivert = 1, nverti
        do ih = 1, nhbound_domain
            nl = nl + 1
            do ns = 1, nspecboun
                if(isboun(ns).gt.0) then
                    boundlat(nl, isboun(ns), 1) = buf3(ns, ih, ivert)
                end if
            end do
        end do
    end do
    
    ncstat = nf90_get_var(&
            bounconc_ncid, bounconc_times_varid, datebuf, (/1, 2/), (/dlen, 1/))
NCERR(__LINE__)
    idex = idate(1)
    idr = mm5date2numeric(datebuf)
    if(idr.ne.idex) then
        print *, '*** ERROR: WRONG EXPECTED DATE IN BOUN_CONCS FILE'
        print *, 'IHOURRUN=1,  EXPECTED= ', idex, ' BOUN_CONCS= ', idr
        call exit1('Exiting')
    endif
    stvec4 = (/1, 1, 1, 2/)
    cntvec4 = (/nspecboun, nhbound_domain, nverti, 1/)
    ncstat = nf90_get_var(bounconc_ncid, latconc_conc_varid, buf3, stvec4, cntvec4)
NCERR(__LINE__)
    
    nl = 0
    do ivert = 1, nverti
        do ih = 1, nhbound_domain
            nl = nl + 1
            do ns = 1, nspecboun
                if(isboun(ns).gt.0) then
                    boundlat(nl, isboun(ns), 2) = buf3(ns, ih, ivert)
                end if
            end do
        end do
    end do
    
    deallocate(buf3)
    
    !print*,'*** Top Boundary Concentrations'
    
    ! Initialization
    do ns = 1, nspec
        do ime = 1, nmerid_domain
            do izo = 1, nzonal_domain
                boundtop(izo, ime, ns, 1) = 1d-17
                boundtop(izo, ime, ns, 2) = 1d-17
            end do
        end do
    end do
    
    !print*,'Reading top concentrations in ppb'
    ncstat = nf90_inq_varid(bounconc_ncid, 'top_conc', topconc_conc_varid)
NCERR(__LINE__)
    stvec4 = (/1, 1, 1, 1/)
    cntvec4 = (/nspecboun, nzonal_domain, nmerid_domain, 1/)
    allocate(buf3(nspecboun, nzonal_domain, nmerid_domain))
    ncstat = nf90_get_var(bounconc_ncid, topconc_conc_varid, buf3, stvec4, cntvec4)
NCERR(__LINE__)
    
    do ime = 1, nmerid_domain
        do izo = 1, nzonal_domain
            do ns = 1, nspecboun
                if(isboun(ns).gt.0) then
                    boundtop(izo, ime, isboun(ns), 1) = buf3(ns, izo, ime)
                end if
            end do
        end do
    end do
    
    stvec4 = (/1, 1, 1, 2/)
    cntvec4 = (/nspecboun, nzonal_domain, nmerid_domain, 1/)
    ncstat = nf90_get_var(bounconc_ncid, topconc_conc_varid, buf3, stvec4, cntvec4)
NCERR(__LINE__)
    
    do ime = 1, nmerid_domain
        do izo = 1, nzonal_domain
            do ns = 1, nspecboun
                if(isboun(ns).gt.0) then
                    boundtop(izo, ime, isboun(ns), 2) = buf3(ns, izo, ime)
                end if
            end do
        end do
    end do
    
    deallocate(buf3)

end subroutine iniboun
