subroutine ainiout

  ! initializes the netcdf output files

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'ainiout.f90')

  !****** parameters
  integer,parameter :: spnamelen = 23
  !  Constant netCDF attributes
  character(len=*),parameter :: title=      'CHIMERE SUITE'
  character(len=*),parameter :: subtitle=   'Hourly avar file'
  character(len=*),parameter :: generator=  'Generated by chimere'
  character(len=*),parameter :: conventions='None'

  !*****************************************************************************************
  ! Local variables
  ! arrays
  real(kind=8),dimension(nzonal_domain,nmerid_domain) :: buf2d
  real(kind=8),dimension(nzonal_domain,nmerid_domain,nverti) :: buf3d
  character(len=spnamelen),dimension(nemisa) :: bufspec
  character(len=spnamelen),dimension(nemisb) :: bufspecb
  character(len=spnamelen),dimension(nspecbounout) :: bufspecbc
  ! strings
  character(len=1024) :: history
  character(len=15)   :: long_name
  ! scalars
  integer :: iospec,ispec,ns,ia,ib,ne

  integer :: ncstat
  integer,dimension(4) :: out_time_dimid
  integer,dimension(4) :: out_date_dimid
  integer ,dimension(4):: out_zonal_dimid
  integer ,dimension(4):: out_merid_dimid
  integer ,dimension(4):: out_layers_dimid
  integer,dimension(4) :: out_verti_dimid
  integer :: out_hb_dimid
  integer,dimension(4) :: out_splen_dimid
  integer,dimension(4) :: out_species_dimid
  
  ! System information stuff
  character(len=32) :: systime
  character(len=64) :: hname
  character(len=32) :: usrname
  character(len=255) :: cwd
  
  allocate(aboun_species(nspecbounout))

  !*****************************************************************************************

  ! 'get various system informations'
  call get_system(usrname,hname,systime,cwd)

  !print*,'create out files'
  print*,afnout
  ncstat=nf90_create(afnout(1:len_trim(afnout))//'aemis.nc',NF90_CLOBBER,aoutea_ncid)
  NCERR(__LINE__)
  if(optemisb.ne.0) then
  ncstat=nf90_create(afnout(1:len_trim(afnout))//'bemis.nc',NF90_CLOBBER,aouteb_ncid)
  NCERR(__LINE__)
  endif
  ncstat=nf90_create(afnout(1:len_trim(afnout))//'bc.nc',NF90_CLOBBER,aoutbc_ncid)
  NCERR(__LINE__)
  ncstat=nf90_create(afnout(1:len_trim(afnout))//'ini.nc',NF90_CLOBBER,aoutini_ncid)
  NCERR(__LINE__)
!  open(79,form='formatted',status='new',action='write')
  
  !print*,'create dimensions in out files'
  call create_dims_out_ea
  if(optemisb.ne.0) call create_dims_out_eb
  call create_dims_out_bc
  call create_dims_out_ini
  
  !print*,'create variables in out file'
  do ispec=1,nspec
    if (inemisa(ispec).ne.0) then
      emisa_species(inemisa(ispec))%name=species(ispec)%name
    endif
    if (inemisb(ispec).ne.0) then
      emisb_species(inemisb(ispec))%name=species(ispec)%name
    endif
  end do
  
  ne = 0
  do ispec=1,nspecboun
    if(isboun(ispec).gt.0) then
      ne = ne + 1
      aboun_species(ne)%name=species(isboun(ispec))%name
    end if
  end do
  
  call create_vars_out_ea
  if(optemisb.ne.0) call create_vars_out_eb
  call create_vars_out_bc
  call create_vars_out_ini
  
  
  !print*,'write attributes to out file'
  call create_atts_out_ea
  if(optemisb.ne.0) call create_atts_out_eb
  call create_atts_out_bc
  call create_atts_out_ini
  
!  print*,'end define mode'
  ncstat=nf90_enddef(aoutea_ncid)
  NCERR(__LINE__)
  if(optemisb.ne.0) then
    ncstat=nf90_enddef(aouteb_ncid)
NCERR(__LINE__)
  endif
  ncstat=nf90_enddef(aoutbc_ncid)
  NCERR(__LINE__)
  ncstat=nf90_enddef(aoutini_ncid)
  NCERR(__LINE__)
  
!  print*,'record longitude and latitude'
  buf2d=xlong
  ncstat=nf90_put_var(aoutea_ncid,out_lon_varid_ad(1),buf2d)
  NCERR(__LINE__)
  if(optemisb.ne.0) then
  ncstat=nf90_put_var(aouteb_ncid,out_lon_varid_ad(4),buf2d)
  NCERR(__LINE__)
  endif
  ncstat=nf90_put_var(aoutbc_ncid,out_lon_varid_ad(2),buf2d)
  NCERR(__LINE__)
  ncstat=nf90_put_var(aoutini_ncid,out_lon_varid_ad(3),buf2d)
  NCERR(__LINE__)
  buf2d=xlati
  ncstat=nf90_put_var(aoutea_ncid,out_lat_varid_ad(1),buf2d)
  NCERR(__LINE__)
  if(optemisb.ne.0) then
  ncstat=nf90_put_var(aouteb_ncid,out_lat_varid_ad(4),buf2d)
  NCERR(__LINE__)  
  endif
  ncstat=nf90_put_var(aoutbc_ncid,out_lat_varid_ad(2),buf2d)
  NCERR(__LINE__) 
  ncstat=nf90_put_var(aoutini_ncid,out_lat_varid_ad(3),buf2d)
  NCERR(__LINE__)
!   print*,' record vcoord coefficient'
  ncstat=nf90_put_var(aoutea_ncid,out_avcoord_varid,avcoord)
  NCERR(__LINE__)
  ncstat=nf90_put_var(aoutea_ncid,out_bvcoord_varid,bvcoord)
  NCERR(__LINE__)
  
!  print *,' record species'
  do ne=1,nemisa
    bufspec(ne)=emisa_species(ne)%name
  enddo  
  ncstat=nf90_put_var(aoutea_ncid,out_species_varid(1),bufspec)
  NCERR(__LINE__)
  if (optemisb.ne.0) then
    do ib=1,nemisb
      bufspecb(ib)=emisb_species(ib)%name
    enddo
    ncstat=nf90_put_var(aouteb_ncid,out_species_varid(4),bufspecb)
NCERR(__LINE__)
  endif
  
  do ispec=1,nspecbounout
    bufspecbc(ispec)=aboun_species(ispec)%name
  end do 
  ncstat=nf90_put_var(aoutbc_ncid,out_species_varid(2),bufspecbc)
  NCERR(__LINE__)

  ncstat=nf90_put_var(aoutini_ncid,out_species_varid(3),bufspecbc)
  NCERR(__LINE__)

 
contains

  !****************************************************
  subroutine create_dims_out_ea
    implicit none

    ncstat=nf90_def_dim(aoutea_ncid,'Time',        NF90_UNLIMITED,out_time_dimid(1))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutea_ncid,'DateStrLen',  dlen,          out_date_dimid(1))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutea_ncid,'west_east',   nzonal_domain,        out_zonal_dimid(1))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutea_ncid,'south_north', nmerid_domain,        out_merid_dimid(1))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutea_ncid,'bottom_top',  nlevemis,       out_layers_dimid(1))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutea_ncid,'vcoord_dim',  nverti,       out_verti_dimid(1))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutea_ncid,'SpStrLen',  spnamelen,  out_splen_dimid(1))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutea_ncid,'Species',    nemisa,  out_species_dimid(1))
NCERR(__LINE__)
  end subroutine create_dims_out_ea
  !***************************************************
    subroutine create_dims_out_eb
    implicit none

    ncstat=nf90_def_dim(aouteb_ncid,'Time',        NF90_UNLIMITED,out_time_dimid(4))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aouteb_ncid,'DateStrLen',  dlen,          out_date_dimid(4))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aouteb_ncid,'west_east',   nzonal_domain,        out_zonal_dimid(4))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aouteb_ncid,'south_north', nmerid_domain,        out_merid_dimid(4))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aouteb_ncid,'bottom_top',  1,       out_layers_dimid(4))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aouteb_ncid,'SpStrLen',  spnamelen,  out_splen_dimid(4))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aouteb_ncid,'Species',    nemisb,  out_species_dimid(4))
NCERR(__LINE__)
  end subroutine create_dims_out_eb
  !****************************************************
  
  subroutine create_dims_out_bc
    implicit none

    ncstat=nf90_def_dim(aoutbc_ncid,'Time',        NF90_UNLIMITED,out_time_dimid(2))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutbc_ncid,'DateStrLen',  dlen,          out_date_dimid(2))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutbc_ncid,'west_east',   nzonal_domain,        out_zonal_dimid(2))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutbc_ncid,'south_north', nmerid_domain,        out_merid_dimid(2))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutbc_ncid,'bottom_top',  nverti,       out_layers_dimid(2))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutbc_ncid,'SpStrLen',  spnamelen,  out_splen_dimid(2))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutbc_ncid,'Species',    nspecbounout,  out_species_dimid(2))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutbc_ncid,'h_boundary',  nhbound_domain,  out_hb_dimid)
NCERR(__LINE__)
  end subroutine create_dims_out_bc
  !  ****************************************************
  subroutine create_dims_out_ini
    implicit none

    ncstat=nf90_def_dim(aoutini_ncid,'Time',        1,out_time_dimid(3))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutini_ncid,'DateStrLen',  dlen,          out_date_dimid(3))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutini_ncid,'west_east',   nzonal_domain,        out_zonal_dimid(3))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutini_ncid,'south_north', nmerid_domain,        out_merid_dimid(3))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutini_ncid,'bottom_top',  nverti,       out_layers_dimid(3))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutini_ncid,'SpStrLen',  spnamelen,  out_splen_dimid(3))
NCERR(__LINE__)
    ncstat=nf90_def_dim(aoutini_ncid,'Species',    nspecbounout,  out_species_dimid(3))
NCERR(__LINE__)
  end subroutine create_dims_out_ini
  !****************************************************
  subroutine create_atts_out_ea
    implicit none

    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'Title',title)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'Sub-title',subtitle)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'Chimere_type','out')
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'Generating_process',generator)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'Conventions',conventions)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'Domain',domain)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'Chimere_version',version)
NCERR(__LINE__)

#if defined(IFORT)

    history = &
         'File '//fnout(1:len_trim(fnout))      // &
         ' was generated on '                   // &
         systime(1:len_trim(systime))           // &
         ' by '                                 // &
         usrname(1:len_trim(usrname))           // &
         ' on '                                 // &
         hname(1:len_trim(hname))               //'\n'C
    history = history(1:len_trim(history)-1)

#elif (defined(G95)+defined(PGI))

    history =                                       &
         'File '//fnout(1:len_trim(fnout))       // &
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

    ncstat=nf90_put_att(aoutea_ncid,NF90_GLOBAL,'history',history(1:len_trim(history)))
NCERR(__LINE__)

  end subroutine create_atts_out_ea
  !****************************************************
  subroutine create_atts_out_eb
    implicit none

    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'Title',title)
NCERR(__LINE__)
    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'Sub-title',subtitle)
NCERR(__LINE__)
    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'Chimere_type','out')
NCERR(__LINE__)
    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'Generating_process',generator)
NCERR(__LINE__)
    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'Conventions',conventions)
NCERR(__LINE__)
    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'Domain',domain)
NCERR(__LINE__)
    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'Chimere_version',version)
NCERR(__LINE__)

#if defined(IFORT)

    history = &
         'File '//fnout(1:len_trim(fnout))      // &
         ' was generated on '                   // &
         systime(1:len_trim(systime))           // &
         ' by '                                 // &
         usrname(1:len_trim(usrname))           // &
         ' on '                                 // &
         hname(1:len_trim(hname))               //'\n'C
    history = history(1:len_trim(history)-1)

#elif (defined(G95)+defined(PGI))

    history =                                       &
         'File '//fnout(1:len_trim(fnout))       // &
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

    ncstat=nf90_put_att(aouteb_ncid,NF90_GLOBAL,'history',history(1:len_trim(history)))
NCERR(__LINE__)

  end subroutine create_atts_out_eb
  !****************************************************
  subroutine create_atts_out_bc
    implicit none

    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'Title',title)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'Sub-title',subtitle)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'Chimere_type','out')
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'Generating_process',generator)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'Conventions',conventions)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'Domain',domain)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'Chimere_version',version)
NCERR(__LINE__)

#if defined(IFORT)

    history = &
         'File '//fnout(1:len_trim(fnout))      // &
         ' was generated on '                   // &
         systime(1:len_trim(systime))           // &
         ' by '                                 // &
         usrname(1:len_trim(usrname))           // &
         ' on '                                 // &
         hname(1:len_trim(hname))               //'\n'C
    history = history(1:len_trim(history)-1)

#elif (defined(G95)+defined(PGI))

    history =                                       &
         'File '//fnout(1:len_trim(fnout))       // &
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

    ncstat=nf90_put_att(aoutbc_ncid,NF90_GLOBAL,'history',history(1:len_trim(history)))
NCERR(__LINE__)

  end subroutine create_atts_out_bc
 !****************************************************
  subroutine create_atts_out_ini
    implicit none

    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'Title',title)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'Sub-title',subtitle)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'Chimere_type','out')
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'Generating_process',generator)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'Conventions',conventions)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'Domain',domain)
NCERR(__LINE__)
    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'Chimere_version',version)
NCERR(__LINE__)

#if defined(IFORT)

    history = &
         'File '//fnout(1:len_trim(fnout))      // &
         ' was generated on '                   // &
         systime(1:len_trim(systime))           // &
         ' by '                                 // &
         usrname(1:len_trim(usrname))           // &
         ' on '                                 // &
         hname(1:len_trim(hname))               //'\n'
    history = history(1:len_trim(history)-1)

#elif (defined(G95)+defined(PGI))

    history =                                       &
         'File '//fnout(1:len_trim(fnout))       // &
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

    ncstat=nf90_put_att(aoutini_ncid,NF90_GLOBAL,'history',history(1:len_trim(history)))
NCERR(__LINE__)

  end subroutine create_atts_out_ini
  !****************************************************
  subroutine create_vars_out_ea
    implicit none

    ! print*,' LON'
    ncstat=nf90_def_var(                                                         &
         aoutea_ncid,                                                               &
         'lon',                                                                  &
         NF90_DOUBLE,                                                             &
         (/out_zonal_dimid(1),out_merid_dimid(1)/),                               &
         out_lon_varid_ad(1)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_lon_varid_ad(1),                                                          &
         'units',                                                                &
         'degrees_east'                                                          &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_lon_varid_ad(1),                                                          &
         'long_name',                                                            &
         'Longitude'                                                             &
         )
NCERR(__LINE__)

    ! print*,' LAT'
    ncstat=nf90_def_var( &
         aoutea_ncid, &
         'lat', &
         NF90_DOUBLE, &
         (/out_zonal_dimid(1),out_merid_dimid(1)/),                                    &
         out_lat_varid_ad(1)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_lat_varid_ad(1),                                                          &
         'units',                                                                &
         'degrees_north'                                                         &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_lat_varid_ad(1),                                                          &
         'long_name',                                                            &
         'Latitude'                                                              &
         )
NCERR(__LINE__)

    ! print*,' AVCOORD'
    ncstat=nf90_def_var( &
         aoutea_ncid, &
         'a_vcoord', &
         NF90_FLOAT, &
         (/out_verti_dimid(1)/),                                    &
         out_avcoord_varid                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_avcoord_varid,                                                          &
         'units',                                                                &
         'no_units'                                                         &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_avcoord_varid,                                                          &
         'long_name',                                                            &
         'A_sigma_coefficient'                                                              &
         )
NCERR(__LINE__)

    ! print*,' BVCOORD'
    ncstat=nf90_def_var( &
         aoutea_ncid, &
         'b_vcoord', &
         NF90_FLOAT, &
         (/out_verti_dimid(1)/),                                    &
         out_bvcoord_varid                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_bvcoord_varid,                                                          &
         'units',                                                                &
         'no_units'                                                         &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutea_ncid,                                                               &
         out_bvcoord_varid,                                                          &
         'long_name',                                                            &
         'B_sigma_coefficient'                                                              &
         )
NCERR(__LINE__)
    
    ! SPECIES
    ncstat=nf90_def_var(aoutea_ncid,'species',NF90_CHAR,&
    (/out_splen_dimid(1),out_species_dimid(1)/),out_species_varid(1))
    
    ! TIMES'
    ncstat=nf90_def_var(                                                         &
         aoutea_ncid,                                                               &
         'Times',                                                                &
         NF90_CHAR,                                                              &
         (/out_date_dimid(1),out_time_dimid(1)/),                                      &
         out_times_varid_ad(1)                                                         &
         )
NCERR(__LINE__)
    
    
    ! EMISA SPECIES
    do ne=1,nemisa
      ! write names
      ncstat=nf90_def_var(                                              &
            aoutea_ncid,                                                     &
            emisa_species(ne)%name,                                         &
            NF90_DOUBLE,                                                          &
            (/out_zonal_dimid(1),out_merid_dimid(1),out_layers_dimid(1),out_time_dimid(1)/), &
            emisa_species(ne)%varid                                      &
            )
       NCERR(__LINE__)
       long_name=emisa_species(ne)%name
       ncstat=nf90_put_att(                                                      &
            aoutea_ncid,                                                            &
            emisa_species(ne)%varid,                                        &
            'long_name',                                                         &
            long_name(1:len_trim(long_name))//' aemissions'                   &
            )
       NCERR(__LINE__)
    end do
  end subroutine create_vars_out_ea
    !****************************************************
  subroutine create_vars_out_eb
    implicit none

!    print*,' LON'
    ncstat=nf90_def_var(                                                         &
         aouteb_ncid,                                                               &
         'lon',                                                                  &
         NF90_DOUBLE,                                                             &
         (/out_zonal_dimid(4),out_merid_dimid(4)/),                               &
         out_lon_varid_ad(4)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aouteb_ncid,                                                               &
         out_lon_varid_ad(4),                                                          &
         'units',                                                                &
         'degrees_east'                                                          &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aouteb_ncid,                                                               &
         out_lon_varid_ad(4),                                                          &
         'long_name',                                                            &
         'Longitude'                                                             &
         )
NCERR(__LINE__)

!    print*,' LAT'
    ncstat=nf90_def_var( &
         aouteb_ncid, &
         'lat', &
         NF90_DOUBLE, &
         (/out_zonal_dimid(4),out_merid_dimid(4)/),                                    &
         out_lat_varid_ad(4)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aouteb_ncid,                                                               &
         out_lat_varid_ad(4),                                                          &
         'units',                                                                &
         'degrees_north'                                                         &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aouteb_ncid,                                                               &
         out_lat_varid_ad(4),                                                          &
         'long_name',                                                            &
         'Latitude'                                                              &
         )
NCERR(__LINE__)

    
!    print*,'SPECIES'
    ncstat=nf90_def_var(aouteb_ncid,'species',NF90_CHAR,&
    (/out_splen_dimid(4),out_species_dimid(4)/),out_species_varid(4))
    
!    print*,'TIMES'
    ncstat=nf90_def_var(                                                         &
         aouteb_ncid,                                                               &
         'Times',                                                                &
         NF90_CHAR,                                                              &
         (/out_date_dimid(4),out_time_dimid(4)/),                                      &
         out_times_varid_ad(4)                                                         &
         )
NCERR(__LINE__)
    
    
!    print*,'EMISB SPECIES'
    do ns=1,nemisb
       ! write names
         ncstat=nf90_def_var(                                              &
            aouteb_ncid,                                                     &
            emisb_species(ns)%name,                                         &
            NF90_DOUBLE,                                                          &
            (/out_zonal_dimid(4),out_merid_dimid(4),out_layers_dimid(4),out_time_dimid(4)/), &
            emisb_species(ns)%varid                                      &
            )
       NCERR(__LINE__)
       long_name=emisb_species(ns)%name
       ncstat=nf90_put_att(                                                      &
            aouteb_ncid,                                                            &
            emisb_species(ns)%varid,                                        &
            'long_name',                                                         &
            long_name(1:len_trim(long_name))//' bemissions'                   &
            )
       NCERR(__LINE__)
    end do
  end subroutine create_vars_out_eb
  !****************************************************
  subroutine create_vars_out_bc
    implicit none

    !print*,' LON'
    ncstat=nf90_def_var(                                                         &
         aoutbc_ncid,                                                               &
         'lon',                                                                  &
         NF90_DOUBLE,                                                             &
         (/out_zonal_dimid(2),out_merid_dimid(2)/),                                    &
         out_lon_varid_ad(2)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutbc_ncid,                                                               &
         out_lon_varid_ad(2),                                                          &
         'units',                                                                &
         'degrees_east'                                                          &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutbc_ncid,                                                               &
         out_lon_varid_ad(2),                                                          &
         'long_name',                                                            &
         'Longitude'                                                             &
         )
NCERR(__LINE__)

    !print*,' LAT'
    ncstat=nf90_def_var( &
         aoutbc_ncid, &
         'lat', &
         NF90_DOUBLE, &
         (/out_zonal_dimid(2),out_merid_dimid(2)/),                                    &
         out_lat_varid_ad(2)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutbc_ncid,                                                               &
         out_lat_varid_ad(2),                                                          &
         'units',                                                                &
         'degrees_north'                                                         &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutbc_ncid,                                                               &
         out_lat_varid_ad(2),                                                          &
         'long_name',                                                            &
         'Latitude'                                                              &
         )
NCERR(__LINE__)

    !print*,' SPECIES'
    ncstat=nf90_def_var(aoutbc_ncid,'species',NF90_CHAR,&
    (/out_splen_dimid(2),out_species_dimid(2)/),out_species_varid(2))
    
    !print*,' TIMES'
    ncstat=nf90_def_var(                                                         &
         aoutbc_ncid,                                                               &
         'Times',                                                                &
         NF90_CHAR,                                                              &
         (/out_date_dimid(2),out_time_dimid(2)/),                                      &
         out_times_varid_ad(2)                                                         &
         )
NCERR(__LINE__)
    
    
    !print*,' TOP CONC'
    ncstat=nf90_def_var(                                              &
            aoutbc_ncid,                                                     &
            'top_conc',                                         &
            NF90_DOUBLE,                                                          &
            (/out_species_dimid(2),out_zonal_dimid(2),out_merid_dimid(2),out_time_dimid(2)/), &
            atopconc_conc_varid                                      &
            )
NCERR(__LINE__)
    
    !print*,' LAT CONC'
    ncstat=nf90_def_var(                                              &
            aoutbc_ncid,                                                     &
            'lat_conc',                                         &
            NF90_DOUBLE,                                                          &
            (/out_species_dimid(2),out_hb_dimid,out_layers_dimid(2),out_time_dimid(2)/), &
            alatconc_conc_varid                                      &
            )
NCERR(__LINE__)

  end subroutine create_vars_out_bc
  !****************************************************
  subroutine create_vars_out_ini
    implicit none

    !print*,' LON'
    ncstat=nf90_def_var(                                                         &
         aoutini_ncid,                                                               &
         'lon',                                                                  &
         NF90_DOUBLE,                                                             &
         (/out_zonal_dimid(3),out_merid_dimid(3)/),                                    &
         out_lon_varid_ad(3)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutini_ncid,                                                               &
         out_lon_varid_ad(3),                                                          &
         'units',                                                                &
         'degrees_east'                                                          &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutini_ncid,                                                               &
         out_lon_varid_ad(3),                                                          &
         'long_name',                                                            &
         'Longitude'                                                             &
         )
NCERR(__LINE__)

    !print*,' LAT'
    ncstat=nf90_def_var( &
         aoutini_ncid, &
         'lat', &
         NF90_DOUBLE, &
         (/out_zonal_dimid(3),out_merid_dimid(3)/),                                    &
         out_lat_varid_ad(3)                                                           &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutini_ncid,                                                               &
         out_lat_varid_ad(3),                                                          &
         'units',                                                                &
         'degrees_north'                                                         &
         )
NCERR(__LINE__)
    ncstat=nf90_put_att(                                                         &
         aoutini_ncid,                                                               &
         out_lat_varid_ad(3),                                                          &
         'long_name',                                                            &
         'Latitude'                                                              &
         )
NCERR(__LINE__)

    !print*,' SPECIES'
    ncstat=nf90_def_var(aoutini_ncid,'species',NF90_CHAR,&
    (/out_splen_dimid(3),out_species_dimid(3)/),out_species_varid(3))
    
    !print*,' TIMES'
    ncstat=nf90_def_var(                                                         &
         aoutini_ncid,                                                               &
         'Times',                                                                &
         NF90_CHAR,                                                              &
         (/out_date_dimid(3),out_time_dimid(3)/),                                      &
         out_times_varid_ad(3)                                                         &
         )
NCERR(__LINE__)
    
    !print*,' ACONCENTRATION FIELDS'
    do ns=1,nspecbounout
       ncstat=nf90_def_var(                                                      &
            aoutini_ncid,                                                            &
            aboun_species(ns)%name,                                        &
            NF90_DOUBLE,                                                          &
            (/out_zonal_dimid(3),out_merid_dimid(3),out_layers_dimid(3)/),  &
            aconc_species(ns)%varid                                               &
            )
       NCERR(__LINE__)
       ncstat=nf90_put_att(                                                  &
            aoutini_ncid,                                                        &
            aconc_species(ns)%varid,                                          &
            'units',                                                         &
            'molecules/cm3'                                                  &
            )
       NCERR(__LINE__)

    end do
  end subroutine create_vars_out_ini

end subroutine ainiout