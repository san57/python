subroutine iniconc_tl
  
  !  Initialization of species concentrations                             
  !  INPUT : ---                                                          
  !  OUTPUT: CONC   Array of concentrations. Only active concs are initial
  !          CONCO  Preceding array of concentrations (used by the TWOSTEP
  
  
  !  Default initialization : Interpolation of boundary values            
  !  Warning: These concentrations are in ppb. 
  
  use chimere_consts
  use wholedomain_common
  use chimere_common
  use master_message_subs
  use netcdf
  implicit none
  
#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'iniconc.f90')

  !******************************************************************************
  integer :: ns, nl
  integer :: iidate, ispec, ivert, ime, izo
  integer :: times, i_nmerid, i_nzonal, i_nverti, i_land, i_species, i_splen
  integer :: dl
  integer :: nbd0
  integer :: ierr
  real(kind = 8) :: fwest, feast, fsouth, fnorth
  real(kind = 8) :: cwest, ceast, csouth, cnorth
  real(kind = 8), allocatable, dimension(:, :, :) :: buf3d, buf3dincr
  real(kind = 8), allocatable, dimension(:, :, :, :) :: buf4d
  
  character(len = dlen) :: cur_date              ! current time
  character(len = 132) :: z_units
  character(len = 23), allocatable, dimension(:) :: buf1d_species
  integer :: start_time_idx                      ! netcdf index of the start time slot
  
  ! netCDF stuff
  integer :: ncid, ncidincr                ! ID of init file, incr file
  integer :: ncid1              ! ID of init file
  integer :: ncstat, ncstati             ! return code for netCDF functions
  integer :: date_dimid         ! Date string dimension IDs
  integer :: spstr_dimid        ! Species string dimension ID
  integer :: time_dimid         ! Time dimension IDs
  integer :: zonal_dimid
  integer :: merid_dimid
  integer :: bt_dimid
  integer :: land_dimid
  integer :: species_dimid
  
  integer :: times_varid        ! Times variable IDs
  integer :: species_varid      ! 
  integer :: varid              ! current variable ID
  
  ! functions
  character(len = dlen), external :: numeric2mm5date
  integer, external :: mm5date2numeric
  
  ! Conversion to dry air mole fraction
  real(kind = 8), allocatable, dimension(:, :, :) :: fdry
  
  
  !***************************************************************************
  conc_tl = dzero
  conc = 1d-17 !!!!!!!!!!!!!!!!!dzero
  
  !  Read in a file
  
  !print*,' open init file',fninit
  ncstat = nf90_open(fninit, NF90_NOWRITE, ncid)
NCERR(__LINE__)
  
  !  read dimensions from init file
  call read_dims
  
  ! allocate storage for a temporary buffer
  allocate(buf3d(nzonal_domain, nmerid_domain, nverti))
  allocate(fdry(nzonal_domain, nmerid_domain, nverti))
  
  ! Verify the presence of the required time slot in init file
  call check_time(start_time_idx)
  
  ! read species in the order CHIMERE waits for
  allocate(buf1d_species(i_species))
  ncstat = nf90_inq_varid(ncid, "species", species_varid)
NCERR(__LINE__)
  ncstat = nf90_get_var(&
      ncid, &
      species_varid, &
      buf1d_species, &
      (/1, 1/), &
      (/i_splen, i_species/) &
      )
NCERR(__LINE__)
  !   conc(1:nspectot,1:nzonal_domain,1:nmerid_domain,1:nverti)=1d-17
  !  Remark: one also initializes conco (used for TWOSTEP)
  !   conco = conc(:,1:nzonal_domain,1:nmerid_domain,1:nverti)
  do ispec = 1, i_species
    call findspec(buf1d_species(ispec), ns)
    if(ns.gt.0) then
      ncstat = nf90_inq_varid(ncid, species(ns)%name, varid)
NCERR(__LINE__)
      ncstat = nf90_get_att(&
          ncid, &
          varid, &
          "units", &
          z_units)
NCERR(__LINE__)
      ! read values for this species and this time slot
      ncstat = nf90_get_var(&
          ncid, &
          varid, &
          buf3d, &
          (/1, 1, 1, start_time_idx/), &
          (/nzonal_domain, nmerid_domain, nverti, 1/) &
          )
NCERR(__LINE__)
      ! copy to "conc" array
      if(z_units(1:3).eq."ppb") then
        conc(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d &
            * airm(1:nzonal_domain, 1:nmerid_domain, 1:nverti, 1) * 1d-9
        if(species(ns)%bounddry.eq.1) then
          fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = &
                  m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6) &
              / (m_air * sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6 &
                  + m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6))
          conc(1:nzonal_domain,1:nmerid_domain,1:nverti) = conc(1:nzonal_domain,1:nmerid_domain,1:nverti) * fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti)
        end if
      else
        conc(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d
      endif
      
      !conco
      ncstat = nf90_inq_varid(ncid, species(ns)%name&
          (1:len_trim(species(ns)%name)) // '_o', varid)
      if (ncstat.lt.0) then
        print*, 'Not a chained simulation, conco=conc for ', species(ns)%name
        conco(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = conc(1:nzonal_domain, 1:nmerid_domain, 1:nverti)
      else
        print*, 'Chained simulation, re-read conco for ', species(ns)%name
        ncstat = nf90_get_att(&
            ncid, &
            varid, &
            "units", &
            z_units)
NCERR(__LINE__)
        ! read values for this species and this time slot
        ncstat = nf90_get_var(&
            ncid, &
            varid, &
            buf3d, &
            (/1, 1, 1, start_time_idx/), &
            (/nzonal_domain, nmerid_domain, nverti, 1/) &
            )
NCERR(__LINE__)
        ! copy to "conco" array
        if(z_units(1:3).eq."ppb") then
          conco(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d &
              * airm(1:nzonal_domain, 1:nmerid_domain, 1:nverti, 1) * 1d-9
          if(species(ns)%bounddry.eq.1) then
            fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = &
                    m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6) &
                / (m_air * sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6 &
                    + m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6))
            conco(1:nzonal_domain,1:nmerid_domain,1:nverti) = conco(1:nzonal_domain,1:nmerid_domain,1:nverti) * fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti)
          end if
        else
          conco(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d
        endif
      endif ! chained simulation
      ! conc_tl
      ncstat = nf90_inq_varid(ncid, &
          species(ns)%name(1:len_trim(species(ns)%name)) // '_tl', varid)
      if (ncstat.lt.0) then
        print*, 'Open increment file', fninitincr
        ncstat = nf90_open(fninitincr, NF90_NOWRITE, ncidincr)
NCERR(__LINE__)
        allocate(buf3dincr(nzonal_domain, nmerid_domain, nverti))
        call check_time_incr(start_time_idx)
        
        ncstati = nf90_inq_varid(ncidincr, species(ns)%name, varid)
        if (ncstati.lt.0) then
          print*, 'No initialization for TL: conc_tl=0.', species(ns)%name
          conc_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = dzero
        else
          print*, 'Initialize TL from increments ', species(ns)%name
          ncstat = nf90_get_att(&
              ncidincr, &
              varid, &
              "units", &
              z_units)
NCERR(__LINE__)
          ncstat = nf90_get_var(&
              ncidincr, &
              varid, &
              buf3dincr, &
              (/1, 1, 1, start_time_idx/), &
              (/nzonal_domain, nmerid_domain, nverti, 1/) &
              )
NCERR(__LINE__)
          ! copy to "conc_tl" array
          if(z_units(1:3).eq."ppb") then
            conc_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3dincr &
                * airm(1:nzonal_domain, 1:nmerid_domain, 1:nverti, 1) * 1d-9
            if(species(ns)%bounddry.eq.1) then
              fdry(1:nzonal_domain,1:nmerid_domain,1:nverti) = m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6) &
                  / (m_air * sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6 &
                      + m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6))
              conc_tl(1:nzonal_domain,1:nmerid_domain,1:nverti) = conc_tl(1:nzonal_domain,1:nmerid_domain,1:nverti) * fdry(1:nzonal_domain,1:nmerid_domain,1:nverti)
            end if
          else
            conc_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3dincr
          endif
        endif ! incr given
        deallocate(buf3dincr)
      else
        print*, 'Initialize TL for a chained simulation', species(ns)%name
        ncstat = nf90_get_att(&
            ncid, &
            varid, &
            "units", &
            z_units)
NCERR(__LINE__)
        ! read values for this species and this time slot
        ncstat = nf90_get_var(&
            ncid, &
            varid, &
            buf3d, &
            (/1, 1, 1, start_time_idx/), &
            (/nzonal_domain, nmerid_domain, nverti, 1/) &
            )
NCERR(__LINE__)
        ! copy to "conc" array
        if(z_units(1:3).eq."ppb") then
          conc_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d &
              * airm(1:nzonal_domain, 1:nmerid_domain, 1:nverti, 1) * 1d-9
          if(species(ns)%bounddry.eq.1) then
            fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6) &
                / (m_air * sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6 &
                    + m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6))
            conc_tl(1:nzonal_domain,1:nmerid_domain,1:nverti) = conc_tl(1:nzonal_domain,1:nmerid_domain,1:nverti) * fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti)
          end if
        else
          conc_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d
        endif
      endif ! chained simu
      
      !conco_tl
      ncstat = nf90_inq_varid(ncid, &
          species(ns)%name(1:len_trim(species(ns)%name)) // '_o_tl', varid)
      if (ncstat.lt.0) then
        ncstati = nf90_inq_varid(ncidincr, species(ns)%name, varid)
        if (ncstati.lt.0) then
          print*, 'No initalization for TL: conco_tl=0. ', species(ns)%name
          conco_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = dzero
        else
          print*, 'Not a chained simulation, conco_tl=conc_tl for ', species(ns)%name
          conco_tl = conc_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti)
        endif
      else
        print*, 'Initialize TL for a chained simulation ', species(ns)%name
        ncstat = nf90_get_att(&
            ncid, &
            varid, &
            "units", &
            z_units)
NCERR(__LINE__)
        ! read values for this species and this time slot
        ncstat = nf90_get_var(&
            ncid, &
            varid, &
            buf3d, &
            (/1, 1, 1, start_time_idx/), &
            (/nzonal_domain, nmerid_domain, nverti, 1/) &
            )
NCERR(__LINE__)
        ! copy to "conco" array
        if(z_units(1:3).eq."ppb") then
          conco_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d &
              * airm(1:nzonal_domain, 1:nmerid_domain, 1:nverti, 1) * 1d-9
          if(species(ns)%bounddry.eq.1) then
            fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6) &
                / (m_air * sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6 &
                    + m_h2o * (1 - sphu(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / airm(1:nzonal_domain,1:nmerid_domain,1:nverti,1) / 1.6))
            conco_tl(1:nzonal_domain,1:nmerid_domain,1:nverti) = conc_tl(1:nzonal_domain,1:nmerid_domain,1:nverti) * fdry(1:nzonal_domain, 1:nmerid_domain, 1:nverti)
          end if
        else
          conco_tl(1:nzonal_domain, 1:nmerid_domain, 1:nverti) = buf3d
        endif
      
      endif
      call duplicate_boundary_tl
      call master_send_ns(ns)
      call master_send_conc
      call master_send_conc_tl
      call mpi_barrier(mpi_comm_world, ierr)
      print*, 'CHECK MPI barriers', rank, 'after barrier 2'
    
    endif
  end do
  
  
  !  Initialization of accumulated variables
  drydep_tl = dzero
  drydep = dzero
  wetdep_tl = dzero
  wetdep = dzero
  
  deallocate(buf3d)
  deallocate(buf1d_species)
  ncstat = nf90_close(ncid)
NCERR(__LINE__)
  call master_send_ns(10000)
  call mpi_barrier(mpi_comm_world, ierr)
  print*, 'CHECK MPI barriers', rank, 'after barrier 3'
contains
  !********************************************************************
  subroutine duplicate_boundary_tl
    ! Duplicate boundary lines and columns
    !     Corners keep their initial value of dzero. They are never used
    
    do ime = 1, nmerid_domain
      conc(0, ime, :) = conc(1, ime, :)
      conc(-1, ime, :) = conc(1, ime, :)
      conc(-2, ime, :) = conc(1, ime, :)
      conc(nzonal_domain + 1, ime, :) = conc(nzonal_domain, ime, :)
      conc(nzonal_domain + 2, ime, :) = conc(nzonal_domain, ime, :)
      conc(nzonal_domain + 3, ime, :) = conc(nzonal_domain, ime, :)
    end do
    do izo = 1, nzonal_domain
      conc(izo, 0, :) = conc(izo, 1, :)
      conc(izo, -1, :) = conc(izo, 1, :)
      conc(izo, -2, :) = conc(izo, 1, :)
      conc(izo, nmerid_domain + 1, :) = conc(izo, nmerid_domain, :)
      conc(izo, nmerid_domain + 2, :) = conc(izo, nmerid_domain, :)
      conc(izo, nmerid_domain + 3, :) = conc(izo, nmerid_domain, :)
    end do
    do ime = 1, nmerid_domain
      conc_tl(0, ime, :) = conc_tl(1, ime, :)
      conc_tl(-1, ime, :) = conc_tl(1, ime, :)
      conc_tl(-2, ime, :) = conc_tl(1, ime, :)
      conc_tl(nzonal_domain + 1, ime, :) = conc_tl(nzonal_domain, ime, :)
      conc_tl(nzonal_domain + 2, ime, :) = conc_tl(nzonal_domain, ime, :)
      conc_tl(nzonal_domain + 3, ime, :) = conc_tl(nzonal_domain, ime, :)
    end do
    do izo = 1, nzonal_domain
      conc_tl(izo, 0, :) = conc_tl(izo, 1, :)
      conc_tl(izo, -1, :) = conc_tl(izo, 1, :)
      conc_tl(izo, -2, :) = conc_tl(izo, 1, :)
      conc_tl(izo, nmerid_domain + 1, :) = conc_tl(izo, nmerid_domain, :)
      conc_tl(izo, nmerid_domain + 2, :) = conc_tl(izo, nmerid_domain, :)
      conc_tl(izo, nmerid_domain + 3, :) = conc_tl(izo, nmerid_domain, :)
    end do
  
  end subroutine duplicate_boundary_tl
  
  !***********
  subroutine read_dims
    ! dimensions IDs
    ncstat = nf90_inq_dimid(ncid, 'Time', time_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid, 'DateStrLen', date_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid, 'SpStrLen', spstr_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid, 'Species', species_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid, 'west_east', zonal_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid, 'south_north', merid_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid, 'bottom_top', bt_dimid)
NCERR(__LINE__)
    
    ! lengths
    ncstat = nf90_inquire_dimension(ncid, time_dimid, len = times)
NCERR(__LINE__)
    
    ncstat = nf90_inquire_dimension(ncid, date_dimid, len = dl)
NCERR(__LINE__)
    if (dl/=dlen) call exit1('iniconc : date format error in init file')
    
    ncstat = nf90_inquire_dimension(ncid, spstr_dimid, len = i_splen)
NCERR(__LINE__)
    
    ncstat = nf90_inquire_dimension(ncid, species_dimid, len = i_species)
NCERR(__LINE__)
    
    ncstat = nf90_inquire_dimension(ncid, zonal_dimid, len = i_nzonal)
NCERR(__LINE__)
    if (i_nzonal/=nzonal_domain) call exit1('iniconc : X dimension error in init file')
    
    ncstat = nf90_inquire_dimension(ncid, merid_dimid, len = i_nmerid)
NCERR(__LINE__)
    if (i_nmerid/=nmerid_domain) call exit1('iniconc : Y dimension error in init file')
    
    ncstat = nf90_inquire_dimension(ncid, bt_dimid, len = i_nverti)
NCERR(__LINE__)
    if (i_nverti/=nverti) call exit1('iniconc : Z dimension error in init file')
  
  end subroutine read_dims
  
  !***********
  subroutine check_time(idx)
    ! checks if start date/time for continued run is included in init file
    ! sets idx as the netcdf index of this start date/time
    
    integer, intent(out) :: idx
    logical :: found
    
    ncstat = nf90_inq_varid(ncid, 'Times', times_varid)
NCERR(__LINE__)
    do idx = 1, times
      ncstat = nf90_get_var(ncid, times_varid, cur_date, (/1, idx/), (/dlen, 1/))
NCERR(__LINE__)
      iidate = mm5date2numeric(cur_date)
      found = (iidate==idatestart)
      if (found) exit
    end do
    if (.not. found) call exit1('iniconc : init file does not contain the requested start time')
  end subroutine check_time
  
  !***********
  subroutine check_time_incr(idx)
    ! checks if start date/time for continued run is included in init file
    ! sets idx as the netcdf index of this start date/time
    
    integer, intent(out) :: idx
    logical :: found
    
    ncstat = nf90_inq_varid(ncidincr, 'Times', times_varid)
NCERR(__LINE__)
    do idx = 1, times
      ncstat = nf90_get_var(ncidincr, times_varid, cur_date, (/1, idx/), (/dlen, 1/))
NCERR(__LINE__)
      iidate = mm5date2numeric(cur_date)
      found = (iidate==idatestart)
      if (found) exit
    end do
    if (.not. found) call exit1('iniconc : init increment file does not contain the requested start time')
  end subroutine check_time_incr
  !***********
  subroutine check_time1(idx)
    ! checks if start date/time for continued run is included in init file
    ! sets idx as the netcdf index of this start date/time
    
    integer, intent(out) :: idx
    logical :: found
    
    ncstat = nf90_inq_varid(ncid1, 'Times', times_varid)
NCERR(__LINE__)
    do idx = 1, times
      ncstat = nf90_get_var(ncid1, times_varid, cur_date, (/1, idx/), (/dlen, 1/))
NCERR(__LINE__)
      iidate = mm5date2numeric(cur_date)
      found = (iidate==idatestart)
      if (found) exit
    end do
    if (.not. found) call exit1('iniconc : init file does not contain the requested start time')
  end subroutine check_time1
  
  subroutine read_dims1
    ! dimensions IDs
    ncstat = nf90_inq_dimid(ncid1, 'Time', time_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid1, 'DateStrLen', date_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid1, 'SpStrLen', spstr_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid1, 'west_east', zonal_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid1, 'south_north', merid_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid1, 'bottom_top', bt_dimid)
NCERR(__LINE__)
    
    ncstat = nf90_inq_dimid(ncid1, 'landuse', land_dimid)
NCERR(__LINE__)
    
    ! lengths
    ncstat = nf90_inquire_dimension(ncid1, time_dimid, len = times)
NCERR(__LINE__)
    
    ncstat = nf90_inquire_dimension(ncid1, date_dimid, len = dl)
NCERR(__LINE__)
    if (dl/=dlen) call exit1('iniconc : date format error in init file')
    
    ncstat = nf90_inquire_dimension(ncid1, zonal_dimid, len = i_nzonal)
NCERR(__LINE__)
    if (i_nzonal/=nzonal_domain) call exit1('iniconc : X dimension error in init file')
    
    ncstat = nf90_inquire_dimension(ncid1, merid_dimid, len = i_nmerid)
NCERR(__LINE__)
    if (i_nmerid/=nmerid_domain) call exit1('iniconc : Y dimension error in init file')
    
    ncstat = nf90_inquire_dimension(ncid1, bt_dimid, len = i_nverti)
NCERR(__LINE__)
    
    ncstat = nf90_inquire_dimension(ncid1, land_dimid, len = i_land)
NCERR(__LINE__)
    if (i_land/=nlduse) call exit1('iniconc : L dimension error in init file')
  
  end subroutine read_dims1


end subroutine iniconc_tl
