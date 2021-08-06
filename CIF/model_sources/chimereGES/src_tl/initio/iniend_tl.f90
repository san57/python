subroutine iniend_tl

  !  Initialization of concentration file

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'iniend.f90')

  !****** parameters
  integer,parameter :: spnamelen = 23

  !  Constant netCDF attributes
  character(len=*),parameter :: title=      'CHIMERE SUITE'
  character(len=*),parameter :: subtitle=   'Final concentrations file'
  character(len=*),parameter :: generator=  'Generated by chimere'
  character(len=*),parameter :: conventions='None'

  !*****************************************************************************************
  ! Local variables
  ! arrays
  real(kind=8),dimension(nzonal_domain,nmerid_domain) :: buf2d
  character(len=spnamelen),dimension(nspectot*4) :: bufspec
  ! strings
  character(len=1024) :: history 
  character(len=spnamelen) :: nametl

  ! scalars
  integer :: ispec
  integer :: ioutend
  integer :: ncstat
  integer :: end_time_dimid
  integer :: end_date_dimid
  integer :: end_splen_dimid
  integer :: end_species_dimid
  integer :: end_zonal_dimid
  integer :: end_merid_dimid
  integer :: end_verti_dimid

  ! System information stuff
  character(len=32) :: systime
  character(len=64) :: hname
  character(len=32) :: usrname
  character(len=255) :: cwd

  !*****************************************************************************************

  ! get various system informations
  call get_system(usrname,hname,systime,cwd)

  ! create "end" file
  ncstat=nf90_create(fnconcs,NF90_CLOBBER,end_ncid)   ; NCERR(__LINE__)

  ! create dimensions in "end" file
  call create_dims_end

  ! write attributes to "end" file
  call create_atts_end

  ! scatter the "conc" array into separate variables
  call create_vars_end

  ! end define mode
  ncstat=nf90_enddef(end_ncid)
  NCERR(__LINE__)

  ! record longitude and latitude
  buf2d=xlong
  ncstat=nf90_put_var(end_ncid,end_lon_varid,buf2d)
  NCERR(__LINE__)
  buf2d=xlati
  ncstat=nf90_put_var(end_ncid,end_lat_varid,buf2d)
  NCERR(__LINE__)

  do ispec=1,nspectot
    bufspec(ispec)=species(ispec)%name
  end do
  do ispec=1,nspectot
    nametl=species(ispec)%name(1:len_trim(species(ispec)%name))//'_tl'
    bufspec(nspectot+ispec)=nametl
  end do 
  do ispec=1,nspectot
    nametl=species(ispec)%name(1:len_trim(species(ispec)%name))//'_o'
    bufspec(nspectot*2+ispec)=nametl
  end do
  do ispec=1,nspectot
    nametl=species(ispec)%name(1:len_trim(species(ispec)%name))//'_o_tl'
    bufspec(nspectot*3+ispec)=nametl
  end do
  ! record species
  ncstat=nf90_put_var(end_ncid,end_species_varid,bufspec)
  NCERR(__LINE__)

  ! Synchronize disk to avoid data loss
  ncstat=nf90_sync(end_ncid)
  NCERR(__LINE__)

 contains

  !****************************************************
  subroutine create_dims_end
    implicit none

    ncstat=nf90_def_dim(end_ncid,'Time',      NF90_UNLIMITED,  end_time_dimid)
NCERR(__LINE__)
    ncstat=nf90_def_dim(end_ncid,'DateStrLen',          dlen,  end_date_dimid)
NCERR(__LINE__)
    ncstat=nf90_def_dim(end_ncid,'SpStrLen',     spnamelen,  end_splen_dimid)
NCERR(__LINE__)
   ! nspectot plus the TL species plus the conco
    ncstat=nf90_def_dim(end_ncid,'Species',         nspectot*4,  end_species_dimid)
NCERR(__LINE__)
    ncstat=nf90_def_dim(end_ncid,'west_east',         nzonal_domain,  end_zonal_dimid)
NCERR(__LINE__)
    ncstat=nf90_def_dim(end_ncid,'south_north',       nmerid_domain,  end_merid_dimid)
NCERR(__LINE__)
    ncstat=nf90_def_dim(end_ncid,'bottom_top',        nverti,  end_verti_dimid)
NCERR(__LINE__)
  end subroutine create_dims_end


  !****************************************************
  subroutine create_atts_end
    implicit none

    ncstat=nf90_put_att(end_ncid,NF90_GLOBAL,'Title',title)
NCERR(__LINE__)
    ncstat=nf90_put_att(end_ncid,NF90_GLOBAL,'Sub-title',subtitle)
NCERR(__LINE__)
    ncstat=nf90_put_att(end_ncid,NF90_GLOBAL,'Chimere_type','end')
NCERR(__LINE__)
    ncstat=nf90_put_att(end_ncid,NF90_GLOBAL,'Generating_process',generator)
NCERR(__LINE__)
    ncstat=nf90_put_att(end_ncid,NF90_GLOBAL,'Conventions',conventions)
NCERR(__LINE__)
    ncstat=nf90_put_att(end_ncid,NF90_GLOBAL,'Domain',domain)
NCERR(__LINE__)

#if defined(IFORT)

    history = &
         'File '//fnconcs(1:len_trim(fnconcs))   // &
         ' was generated on '                    // &
         systime(1:len_trim(systime))            // &
         ' by '                                  // &
         usrname(1:len_trim(usrname))            // &
         ' on '                                  // &
         hname(1:len_trim(hname))                //'\n'C
    history = history(1:len_trim(history)-1)

#elif (defined(G95)+defined(PGI))

    history =                                       &
         'File '//fnconcs(1:len_trim(fnconcs))   // &
         ' was generated on '                    // &
         systime(1:len_trim(systime))            // &
         ' by '                                  // &
         usrname(1:len_trim(usrname))            // &
         ' on '                                  // &
         hname(1:len_trim(hname))
    history = history(1:len_trim(history)-1)

#else

    history='unknown'

#endif

    ncstat=nf90_put_att(end_ncid,NF90_GLOBAL,'history',history(1:len_trim(history)))
NCERR(__LINE__)

  end subroutine create_atts_end

  !****************************************************
  subroutine create_vars_end
    implicit none

    ! LON
    ncstat=nf90_def_var(                      &
         end_ncid,                            &
         'lon',                               &
         NF90_FLOAT,                          &
         (/end_zonal_dimid,end_merid_dimid/), &
         end_lon_varid                        &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                      &
         end_ncid,                            &
         end_lon_varid,                       &
         'units',                             &
         'degrees_east'                       &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                      &
         end_ncid,                            &
         end_lon_varid,                       &
         'long_name',                         &
         'Longitude'                          &
         )
NCERR(__LINE__)

    ! LAT
    ncstat=nf90_def_var( &
         end_ncid, &
         'lat', &
         NF90_FLOAT, &
         (/end_zonal_dimid,end_merid_dimid/), &
         end_lat_varid                        &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                      &
         end_ncid,                            &
         end_lat_varid,                       &
         'units',                             &
         'degrees_north'                      &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                      &
         end_ncid,                            &
         end_lat_varid,                       &
         'long_name',                         &
         'Latitude'                           &
         )
NCERR(__LINE__)

    ! SPECIES
    ncstat=nf90_def_var( &
         end_ncid,       &
         'species', &
         NF90_CHAR,      &
         (/end_splen_dimid,end_species_dimid/),  &
         end_species_varid                       &
         )
NCERR(__LINE__)

    ! TIMES
    ncstat=nf90_def_var(                         &
         end_ncid,                               &
         'Times',                                &
         NF90_CHAR,                              &
         (/end_date_dimid,end_time_dimid/),      &
         end_times_varid                         &
         )
NCERR(__LINE__)

    ! CONCENTRATION FIELDS
    do ioutend=1,nspectot
       ncstat=nf90_def_var(                                                      &
            end_ncid,                                                            &
            species(ioutend)%name,                                               &
            NF90_DOUBLE,                                                          &
            (/end_zonal_dimid,end_merid_dimid,end_verti_dimid,end_time_dimid/),  &
            species(ioutend)%varid                                               &
            )
       NCERR(__LINE__)
       ncstat=nf90_put_att(                                                  &
            end_ncid,                                                        &
            species(ioutend)%varid,                                          &
            'units',                                                         &
            'molecules/cm3'                                                  &
            )
       NCERR(__LINE__)
       nametl=species(ioutend)%name(1:len_trim(species(ioutend)%name))//'_tl'
       ! TL FIELDS
        ncstat=nf90_def_var(                                                      &
              end_ncid,                                                            &
           nametl,                                        &
           NF90_DOUBLE,                                                          &
         (/end_zonal_dimid,end_merid_dimid,end_verti_dimid,end_time_dimid/),  &
            species(ioutend)%varid_tl                                            &
              )
       NCERR(__LINE__)
         ncstat=nf90_put_att(                                                  &
              end_ncid,                                                        &
              species(ioutend)%varid_tl,                                       &
              'units',                                                         &
              'molecules/cm3'                                                  &
              )
        NCERR(__LINE__)
	! les noms!
       nametl=species(ioutend)%name(1:len_trim(species(ioutend)%name))//'_o'
       !CONCO FIELDS
        ncstat=nf90_def_var(                                                      &
              end_ncid,                                                            &
           nametl,                                        &
           NF90_DOUBLE,                                                          &
         (/end_zonal_dimid,end_merid_dimid,end_verti_dimid,end_time_dimid/),  &
            species(ioutend)%varido                                            &
              )
       NCERR(__LINE__)
         ncstat=nf90_put_att(                                                  &
              end_ncid,                                                        &
              species(ioutend)%varido,                                       &
              'units',                                                         &
              'molecules/cm3'                                                  &
              )
        NCERR(__LINE__)
       ! les noms!
       nametl=species(ioutend)%name(1:len_trim(species(ioutend)%name))//'_o_tl'
       !CONCO FIELDS
        ncstat=nf90_def_var(                                                      &
              end_ncid,                                                            &
           nametl,                                        &
           NF90_DOUBLE,                                                          &
         (/end_zonal_dimid,end_merid_dimid,end_verti_dimid,end_time_dimid/),  &
            species(ioutend)%varido_tl                                         &
              )
       NCERR(__LINE__)
         ncstat=nf90_put_att(                                                  &
              end_ncid,                                                        &
              species(ioutend)%varido_tl,                                       &
              'units',                                                         &
              'molecules/cm3'                                                  &
              )
        NCERR(__LINE__)
    end do


  end subroutine create_vars_end

end subroutine iniend_tl