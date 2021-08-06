subroutine ainiaconc 

  !  Initialization of species adjointized concentrations                             
  !  INPUT : ---                                                          
  !  OUTPUT: ACONC   Array of adjointized concentrations. Only active concs are initial
  !          ACONCO  Preceding array of adjointized concentrations (used by the TWOSTEP)

  use chimere_consts
  use chimere_common
  use wholedomain_common
  use netcdf

  
  implicit none 

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'ainiconc.f90')  

  integer :: iidate,ispec,ns
  integer :: times,dl,i_nmerid,i_nzonal,i_nverti,i_species,i_splen
  real(kind=8),allocatable,dimension(:,:,:) :: buf3d
  real(kind=8),allocatable,dimension(:,:,:) :: buf3do
  character(len=dlen) :: cur_date              ! current time
  character(len=23),allocatable, dimension(:) :: buf1d_species 
  integer,parameter :: spnamelen = 23
  character(len=spnamelen) :: nametl
  
! netCDF stuff
  integer :: ncid               ! ID of init file
  integer :: ncstat             ! return code for netCDF functions
  integer :: time_dimid         ! Time dimension IDs
  integer :: date_dimid         ! Date string dimension IDs
  integer :: spstr_dimid        ! Species string dimension ID
  integer :: species_dimid
  integer :: zonal_dimid
  integer :: merid_dimid 
  integer :: bt_dimid 
  integer :: start_time_idx                      ! netcdf index of the start time slot
  integer :: times_varid        ! Times variable IDs
  integer :: species_varid      ! 
  integer :: varid              ! current variable ID
  ! functions
  character(len=dlen),external :: numeric2mm5date
  integer,external             ::  mm5date2numeric
  
  
  aconc=dzero
  aconco=dzero	
  aconcini=dzero
  aconcsave=dzero
  if (iopainit.ne.0) then
    ! open init file
    fnainit='aini.nc'
    print*,'open file', fnainit
    ncstat=nf90_open(fnainit,NF90_NOWRITE,ncid)
NCERR(__LINE__)
    
    ! read dimensions from init file
    call read_dims
    
    ! allocate storage for a temporary buffer
    allocate(buf3d(nzonal_domain+(6*nzdoms),nmerid_domain+(6*nmdoms),nverti+1))
    allocate(buf3do(nzonal_domain,nmerid_domain,nverti))
    
    ! Verify the presence of the required time slot in init file
    call check_time(start_time_idx)
    
    ! read species in the order CHIMERE waits for
    allocate(buf1d_species(i_species))
    ncstat=nf90_inq_varid(ncid,"species",species_varid)
NCERR(__LINE__)
    ncstat=nf90_get_var(            &
            ncid,                      &
            species_varid,                     &
            buf1d_species,                     &
            (/1,1/),  &
            (/i_splen,i_species/) &
            )
NCERR(__LINE__)
    do ispec=1,i_species
      call findspec(buf1d_species(ispec),ns)
      if(ns.gt.0) then 
        nametl=species(ns)%name(1:len_trim(species(ns)%name))//'_ad'
        ncstat=nf90_inq_varid(ncid,nametl,varid)
        NCERR(__LINE__)
        ! read values for this species and this time slot
        ncstat=nf90_get_var(            &
             ncid,                      &
             varid,                     &
             buf3d,                     &
             (/1,1,1,start_time_idx/),  &
             (/nzonal_domain+(6*nzdoms),nmerid_domain+(6*nmdoms),nverti+1,1/) &
             )
        NCERR(__LINE__)
        ! copy to "conc" array
        aconcsave(ns,:,:,:)=buf3d(:,:,:)
        nametl=species(ns)%name(1:len_trim(species(ns)%name))//'_ad_o'
        ncstat=nf90_inq_varid(ncid,nametl,varid)
        NCERR(__LINE__)
        ! read values for this species and this time slot
        ncstat=nf90_get_var(            &
             ncid,                      &
             varid,                     &
             buf3do,                    &
             (/1,1,1,start_time_idx/),  &
             (/nzonal_domain,nmerid_domain,nverti,1/) &
             )
        NCERR(__LINE__)
       ! copy to "conc" array
      aconco(ns,:,:,:)=buf3do
      endif
    end do
    deallocate(buf3d)
    deallocate(buf3do)
    deallocate(buf1d_species)
    ncstat=nf90_close(ncid)
NCERR(__LINE__)

    
  endif

contains
  !********************************************************************

  subroutine read_dims
    ! dimensions IDs
    ncstat=nf90_inq_dimid(ncid,'Time',         time_dimid)
NCERR(__LINE__)

    ncstat=nf90_inq_dimid(ncid,'DateStrLen',   date_dimid)
NCERR(__LINE__)
    
    ncstat=nf90_inq_dimid(ncid,'SpStrLen',    spstr_dimid)
NCERR(__LINE__)

    ncstat=nf90_inq_dimid(ncid,'Species',   species_dimid)
NCERR(__LINE__)

    ncstat=nf90_inq_dimid(ncid,'west_east',   zonal_dimid)
NCERR(__LINE__)

    ncstat=nf90_inq_dimid(ncid,'south_north', merid_dimid)
NCERR(__LINE__)

    ncstat=nf90_inq_dimid(ncid,'bottom_top',     bt_dimid)
NCERR(__LINE__)

    ! lengths
    ncstat=nf90_inquire_dimension(ncid,time_dimid,  len=times)
NCERR(__LINE__)

    ncstat=nf90_inquire_dimension(ncid,date_dimid,  len=dl)
NCERR(__LINE__)
    if (dl/=dlen) stop 'iniconc : date format error in init file'

    ncstat=nf90_inquire_dimension(ncid,spstr_dimid,  len=i_splen)
NCERR(__LINE__)

    ncstat=nf90_inquire_dimension(ncid,species_dimid,  len=i_species)
NCERR(__LINE__)

    ncstat=nf90_inquire_dimension(ncid,zonal_dimid, len=i_nzonal)
NCERR(__LINE__)
    if (i_nzonal/=nzonal_domain+(6*nzdoms)) stop 'ainiaconc : X dimension error in init file'    

    ncstat=nf90_inquire_dimension(ncid,merid_dimid, len=i_nmerid)
NCERR(__LINE__)
    if (i_nmerid/=nmerid_domain+(6*nmdoms)) stop 'ainiaconc : Y dimension error in init file'    

    ncstat=nf90_inquire_dimension(ncid,bt_dimid, len=i_nverti)
NCERR(__LINE__)
    if (i_nverti/=nverti+1) stop 'ainiaconc : Z dimension error in init file'    

  end subroutine read_dims  

!***********
  subroutine check_time(idx)
    ! checks if start date/time for continued run is included in init file
    ! sets idx as the netcdf index of this start date/time

    integer,intent(out) :: idx
    logical :: found

    ncstat=nf90_inq_varid(ncid,'Times',times_varid)
NCERR(__LINE__)
    do idx=1,times
       ncstat=nf90_get_var(ncid,times_varid,cur_date,(/1,idx/),(/dlen,1/))
       NCERR(__LINE__)
       iidate=mm5date2numeric(cur_date)
       found=(iidate==idate(nhourrun))
       if (found) exit
    end do
    if (.not. found) stop 'ainiaconc : ainit file does not contain the requested start time'
  end subroutine check_time


end subroutine ainiaconc
