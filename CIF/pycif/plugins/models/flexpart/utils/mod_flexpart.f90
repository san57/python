!====================================================================
!
! Fortran subroutines used to create python module 'mod_flexpart'
! containing functions for reading flexpart binary output
!
! Taken from flexinvert+ (mod_flexpart.f90)
!
! Compile this by
! f2py -c mod_flexpart.f90 -m mod_flexpart
!
! TODO: add makefile/automatic compile on setup
!
!====================================================================


! --------------------------------------------------
! read_header
! --------------------------------------------------
!> read_header
!! 
!! Purpose:  Reads the flexpart header file and
!!           initiates general flexpart variables.
!!
! --------------------------------------------------

subroutine read_header(filename, numpoint, ibdate, ibtime, releases, &
     numx, numy, bndx, bndy, delx, dely, xshift, ageclass, trajdays, ndgrid, outheight)


  implicit none

  integer, parameter :: maxpoint=1000, maxspec=10, maxlev=50

  character(len=256),  intent(in)  :: filename
  integer,                      intent(out) :: numpoint
  integer,                      intent(out) :: ibdate, ibtime   
  integer, dimension(maxpoint), intent(out) :: releases           
  integer,                      intent(out) :: numx, numy, xshift
  integer,                      intent(out) :: ageclass ! length of backtrajectories (s)
  integer,                      intent(out) :: trajdays ! length of backtrajectories (days)
  integer,                      intent(out) :: ndgrid   ! number of footprints per day
  real,                         intent(out) :: bndx, bndy, delx, dely
  real, dimension(maxlev),      intent(out) :: outheight ! height of model layers for grid_init files (m)



  character(len=256)                           :: flexversion
  logical                                      :: lexist
  integer                                      :: nageclass, imethod, nspec
  integer                                      :: loutstep, loutaver, loutsample
  integer                                      :: jjjjmmdd, ihmmss
  integer            :: nzgrid                   ! number of vertical layers in grid_init files
  real                                         :: xp2, yp2
  character(len=10), dimension(maxspec)        :: species
  integer, dimension(:), allocatable           :: ireleasestart, ireleaseend
  integer, dimension(:), allocatable           :: npart, nkind
  character(len=45), dimension(:), allocatable :: compoint
  real, dimension(:), allocatable              :: xpoint, ypoint, zpoint1, zpoint2
  real, dimension(:,:), allocatable            :: xmass
  integer, dimension(10)                       :: lage
  integer                                      :: n, i, j, ierr 

  inquire ( file=trim(filename),exist=lexist )
  if ( .not.lexist ) then
    write(*,*) 'ERROR: cannot find: '//trim(filename)
    stop
  endif
  open(100,file=filename,form='unformatted',action='read',status='old',iostat=ierr)
  if(ierr.gt.0) then
    write(*,*) 'ERROR: cannot read: '//trim(filename)
    stop
  endif

  nzgrid=1

  read(100) ibdate, ibtime, flexversion
  read(100) loutstep, loutaver, loutsample
  read(100) bndx, bndy, numx, numy, delx, dely
  read(100) nzgrid, (outheight(i), i = 1, nzgrid)
  read(100) jjjjmmdd, ihmmss
  read(100) nspec, nzgrid
  nspec=nspec/3
  do n=1,nspec
    read(100) nzgrid, species(n)
    read(100) nzgrid, species(n)
    read(100) nzgrid, species(n)
  end do
  read(100) numpoint

  allocate( ireleasestart(numpoint) )
  allocate( ireleaseend(numpoint) )
  allocate( xpoint(numpoint) )
  allocate( ypoint(numpoint) )
  allocate( zpoint1(numpoint) )
  allocate( zpoint2(numpoint) )
  allocate( npart(numpoint) )
  allocate( nkind(numpoint) )
  allocate( xmass(numpoint,nspec) )
  allocate( compoint(numpoint) )

  do i = 1, numpoint
    read(100) ireleasestart(i), ireleaseend(i)
    read(100) xpoint(i), ypoint(i), xp2, yp2, zpoint1(i), zpoint2(i)
    read(100) npart(i), nkind(i)
    read(100) compoint(i)
    do j=1,nspec
      read(100)
      read(100)
      read(100) xmass(i,j)
    end do
  end do
  read(100) imethod
  read(100) nageclass,(lage(i), i = 1, nageclass)
  close(100)

  releases(1:numpoint) = ireleasestart(:)
  ageclass = lage(1)
  trajdays = ageclass/3600/24
  ndgrid = abs(24*3600/loutstep)

  if( real(numx)*delx.eq.360. ) then
! global domain - check offset relative to 180W
    xshift = int((bndx + 180.)/delx) 
    bndx = -180.
  else 
! not global domain
    xshift = 0
  endif


end subroutine read_header


! --------------------------------------------------
! read_grid
! --------------------------------------------------
!> read_grid
!!
!! Purpose:  Reads the grid_time files where there
!!           is one grid_time file per release.
!!
!! ESO: Remove interaction with 'filedates' here, not needed
! --------------------------------------------------

subroutine read_grid(grid, filename, filedates, numx, numy, maxngrid, &
     & xshift, ngrid,  gtime, ndgrid, trajdays)


  implicit none


  character(len=256),                  intent(in)     :: filename, filedates
  integer,                             intent(in)     :: numx, numy, xshift
  integer,                             intent(in)     :: ndgrid, trajdays
  integer,                             intent(out)    :: ngrid
  real, dimension(numx,numy,maxngrid), intent(out)    :: grid
  real(kind=8), dimension(maxngrid),   intent(out)    :: gtime

  integer, intent(in)                                 :: maxngrid
  real, parameter                                     :: scaleconc=1.e12
  real, parameter                                     :: smallnum=1.e-38
  logical                                             :: lexist
  integer                                             :: ierr
  integer                                             :: nread
  integer                                             :: sp_count_i, sp_count_r
  integer                                             :: jjjjmmdd, hhmiss
  integer                                             :: ii, ir, n, nt, jy, ix 
  real                                                :: fact
  integer, allocatable, dimension(:)                       :: sparse_dump_i
  real, allocatable, dimension(:)                          :: sparse_dump_r

  real, dimension(:,:,:), allocatable                 :: work
  integer(kind=4), dimension(:), allocatable          :: times
  real(kind=8), dimension(:), allocatable             :: dates
  real(kind=8)                                        :: jdrel,jdate
  real(kind=8), dimension(:), allocatable             :: jdtime


! read footprint dates
  inquire ( file=trim(filedates),exist=lexist )
  if ( .not.lexist ) then
    write(*,*) 'ERROR: cannot find '//trim(filedates)
    stop
  endif
  call get_nread(filedates, nread)
  allocate( dates(nread) )
  call read_dates(filedates, nread, dates)

! initialize
!  maxngrid=ndgrid*trajdays+2

  allocate( work(numx,numy,maxngrid) )
  allocate( times(maxngrid) )
  allocate( jdtime(maxngrid) )
!  allocate( grid(numx,numy,maxngrid))
!  allocate( gtime(maxngrid) )
  allocate( sparse_dump_i(numx*numy), sparse_dump_r(numx*numy) )

  work(:,:,:) = 0.
  grid(:,:,:) = 0.
  times(:) = 0
  gtime(:) = 0.

! open grid_time file
  inquire ( file=trim(filename),exist=lexist )
  if ( .not.lexist ) then
    write(*,*) 'WARNING: cannot find '//trim(filename)
    go to 10
  endif
  open(100,file=trim(filename),form='unformatted',action='read',status='old',iostat=ierr)
  if( ierr.ne.0 ) then
    write(*,*) 'WARNING: cannot read '//trim(filename)
    go to 10
  endif

! read footprints
  do nt = 1, maxngrid
    read(100,iostat=ierr,end=20) jjjjmmdd
    read(100,iostat=ierr,end=20) hhmiss
! note: assumes loutaverage of at least 1 hour
    times(nt)=jjjjmmdd*100+hhmiss/10000
    fact = 1.
    read(100,iostat=ierr,end=20) sp_count_i
    read(100,iostat=ierr,end=20) (sparse_dump_i(ix), ix=1,sp_count_i)
    read(100,iostat=ierr,end=20) sp_count_r
    read(100,iostat=ierr,end=20) (sparse_dump_r(ix), ix=1,sp_count_r)
    ii = 0
    do ir = 1, sp_count_r
      if ((sparse_dump_r(ir)*fact).gt.smallnum) then
        ii = ii + 1
        n = sparse_dump_i(ii)
        fact = fact*(-1.)
      else
        n = n + 1
      endif
      jy = (n - numx * numy)/numx
      ix = n - numx * numy - numx * jy
      work(ix+1,jy+1,nt) = abs(sparse_dump_r(ir)) * scaleconc
    end do
  end do
20 continue 
  !if ( ierr.ne.0 ) write(*,fmt='(A33,1X,I3)') 'WARNING: read_grid stopped at nt ',nt
  ngrid = nt - 1

  close(100)

! shift grid by xshift along longitudinal dimension
  work = cshift(work, shift=-1*xshift, dim=1)

! calculate footprint time in julian days
  do nt = 1, ngrid
    hhmiss = (times(nt) - (times(nt)/100)*100)*10000
    call juldate(jdate, times(nt)/100,hhmiss)
    jdtime(nt) = jdate
  end do

! index to footprints for current release
! jd is julian day of release
! ndgrid is number of grids per day
! eso: we don't need for now
!  jdrel = jd - 1./real(ndgrid)

! reverse order in time dimension
  do nt = 1, ngrid
    n = ngrid - nt + 1
    grid(:,:,n) = work(:,:,nt)
    gtime(n) = jdtime(nt)
  end do

10 continue


end subroutine read_grid

! --------------------------------------------------
! get_nread
! --------------------------------------------------
!> get_nread
!!
!! Purpose:  Gets the number of time steps in the
!!           grid_time files.
!!
! --------------------------------------------------

subroutine get_nread(filename, nread)

  implicit none

  character(len=256),          intent(in)     :: filename
  integer,                              intent(out) :: nread

  character(len=14), dimension(1000)                   :: tmp
  logical                                              :: lexist
  integer                                              :: ierr, n

  inquire ( file=trim(filename),exist=lexist )
  if ( .not.lexist ) then
    write(*,*) 'ERROR: cannot find '//trim(filename)
    stop
  endif
  open(100,file=trim(filename),action='read',status='old',iostat=ierr)
  if( ierr.ne.0 ) then
    write(*,*) 'ERROR: cannot read '//trim(filename)
    stop
  endif

  n = 1
  do while (ierr.eq.0)
    read(100,fmt='(A14)',iostat=ierr,end=10) tmp(n)
    n = n + 1
  enddo
10 continue

  close(100)

  nread = n - 1

end subroutine get_nread

! --------------------------------------------------
! read_dates
! --------------------------------------------------
!> read_dates
!! 
!! Purpose:  Reads the dates for each time step in
!!           the grid_time files.
!!
! --------------------------------------------------

subroutine read_dates(filename, nread, dates)

  implicit none

  character(len=256),          intent(in)     :: filename
  integer,                              intent(in)     :: nread
  real(kind=8), dimension(nread),       intent(in out) :: dates

  character(len=14)                                    :: tmp
  integer                                              :: ierr, n

  open(100,file=filename,action='read',status='old',iostat=ierr)
  if( ierr.ne.0 ) then
    write(*,*) 'ERROR: cannot read: '//trim(filename)
    stop
  endif

  do n = 1, nread
    read(100,fmt='(A14)',iostat=ierr,end=10) tmp
    read(tmp,*) dates(n)
  enddo
10 continue

  close(100)

end subroutine read_dates


subroutine caldate(jdate, yyyymmdd, hhmiss)

  real(kind=8), intent(in) :: jdate
  integer, intent(out)     :: yyyymmdd, hhmiss
  integer                  :: yyyy, mm, dd, hh, mi, ss
  integer                  :: julday, ja, jb, jc, jd, je, jalpha
  integer, parameter       :: igreg = 2299161

  julday=int(jdate)
  if(julday.ge.igreg) then
    jalpha=int(((julday-1867216)-0.25)/36524.25)
    ja=julday+1+jalpha-int(0.25*jalpha)
  else
    ja=julday
  endif
  jb=ja+1524
  jc=int(6680.+((jb-2439870)-122.1)/365.25)
  jd=365*jc+int(0.25*jc)
  je=int((jb-jd)/30.6001)
  dd=jb-jd-int(30.6001*je)
  mm=je-1
  if(mm.gt.12) mm=mm-12
  yyyy=jc-4715
  if(mm.gt.2) yyyy=yyyy-1
  if(yyyy.le.0) yyyy=yyyy-1

  yyyymmdd=10000*yyyy+100*mm+dd
  hh=int(24.*(jdate-float(julday)))
  mi=int(1440.*(jdate-float(julday))-60.*float(hh))
  ss=nint(86400.*(jdate-float(julday))-3600.*float(hh))-60.*float(mi)
  if(ss.eq.60) then  ! 60 seconds = 1 minute
    ss=0
    mi=mi+1
  endif
  if(mi.eq.60) then
    mi=0
    hh=hh+1
  endif
  hhmiss=10000*hh+100*mi+ss

end subroutine caldate

! --------------------------------------------------
! julday
! --------------------------------------------------
!> julday
!!
!! Purpose:  Converts the date and time format 
!!           YYYYMMDD and HHMMSS to a julian date 
!!           number.
!!
! --------------------------------------------------

subroutine juldate(outdate,yyyymmdd,hhmiss)

  real(kind=8), intent(out) :: outdate
  integer, intent(in) :: yyyymmdd,hhmiss
  integer :: yyyy,mm,hh,dd,mi,ss
  integer :: julday,jy,jm,ja
  integer, parameter :: igreg=15+31*(10+12*1582)

  yyyy=yyyymmdd/10000
  mm=(yyyymmdd-10000*yyyy)/100
  dd=yyyymmdd-10000*yyyy-100*mm
  hh=hhmiss/10000
  mi=(hhmiss-10000*hh)/100
  ss=hhmiss-10000*hh-100*mi

  if(yyyy.eq.0) print*, 'ERROR juldate: there is no year zero'
  if(yyyy.lt.0) yyyy=yyyy+1
  if(mm.gt.2) then
    jy=yyyy
    jm=mm+1
  else
    jy=yyyy-1
    jm=mm+13
  endif
  julday=int(365.25*jy)+int(30.6001*jm)+dd+1720995
  if (dd+31*(mm+12*yyyy).ge.igreg) then
    ja=int(0.01*jy)
    julday=julday+2-ja+int(0.25*ja)
  endif

  outdate=dble(float(julday))+dble(float(hh)/24.)&
       &+dble(float(mi)/1440.)+dble(float(ss)/86400.)

end subroutine juldate

! --------------------------------------------------
! calceomday
! --------------------------------------------------
!> calceomday
!! 
!! Purpose:  Calculates number of days in a given 
!!           year, month. Currently only considers 
!!           years from 1900
! --------------------------------------------------

integer function calceomday(yyyymm)

  integer, intent(in) :: yyyymm
  integer :: yyyy,mm
  integer, dimension(12) :: leapdays,days
  integer :: eomday

  leapdays=(/31,29,31,30,31,30,31,31,30,31,30,31/)
  days=(/31,28,31,30,31,30,31,31,30,31,30,31/)

  yyyy=floor(yyyymm/100.)
  mm=yyyymm-yyyy*100

  if((float(yyyy)/100.).eq.float(yyyy/100)) then
    if((float(yyyy)/400.).eq.float(yyyy/400)) then
      eomday=leapdays(mm)
    else
      eomday=days(mm)
    endif
  else
    if((float(yyyy)/4.).eq.float(yyyy/4)) then
      eomday=leapdays(mm)
    else
      eomday=days(mm)
    endif
  endif

  calceomday=eomday

end function calceomday


! --------------------------------------------------
! read_factor
! --------------------------------------------------
!> read_factor
!!
!! Purpose:  Reads the correction factor (i.e. ratio
!!           of the wet to dry air density for 
!!           correcting mixing ratios for dry air).
!!
! --------------------------------------------------

subroutine read_factor(filename, nread, numx, numy, xshift, factor)

!  use mod_var

  implicit none

  character(len=256),          intent(in)     :: filename
  integer,                              intent(in)     :: nread
  integer,                              intent(in)     :: numx, numy, xshift
  real, dimension(numx,numy,nread),     intent(in out) :: factor

  real, parameter                                      :: smallnum=1.e-38
  logical                                              :: lexist
  integer                                              :: ierr
  integer                                              :: sp_count_i, sp_count_r
  integer                                              :: ii, ir, n, nt, jy, ix
  real                                                 :: fact
  integer, dimension(numx*numy)                        :: sparse_dump_i
  real, dimension(numx*numy)                           :: sparse_dump_r

! open factor file
  inquire ( file=trim(filename),exist=lexist )
  if ( .not.lexist ) then
    write(*,*) 'WARNING: cannot find '//trim(filename)
    go to 10
  endif
  open(100,file=trim(filename),form='unformatted',action='read',status='old',iostat=ierr)
  if( ierr.ne.0 ) then
    write(*,*) 'WARNING: cannot read '//trim(filename)
    go to 10
  endif

! read footprints
  do nt = 1, nread
    fact = 1.
    read(100) sp_count_i
    read(100) (sparse_dump_i(ix), ix=1,sp_count_i)
    read(100) sp_count_r
    read(100) (sparse_dump_r(ix), ix=1,sp_count_r)
    ii = 0
    do ir = 1, sp_count_r
      if ((sparse_dump_r(ir)*fact).gt.smallnum) then
        ii = ii + 1
        n = sparse_dump_i(ii)
        fact = fact*(-1.)
      else
        n = n + 1
      endif
      jy = (n - numx * numy)/numx
      ix = n - numx * numy - numx * jy
      factor(ix+1,jy+1,nt) = abs(sparse_dump_r(ir))
    end do
  end do

  close(100)

! shift grid by xshift along longitudinal dimension
  factor = cshift(factor, shift=-1*xshift, dim=1)

10 continue

end subroutine read_factor


! --------------------------------------------------
! read_init
! --------------------------------------------------
!> read_init
!! 
!! Purpose:  Reads the grid_initial files where 
!!           there is one file per release.
!!
! --------------------------------------------------

subroutine read_init(filename, gridinit, nxgrid,nygrid,nzgrid)

!  use mod_var

  implicit none

  character(len=256),                    intent(in)     :: filename
  integer, intent(in)                                   :: nxgrid, nygrid, nzgrid
  real, dimension(nxgrid,nygrid,nzgrid), intent(in out) :: gridinit

  real, parameter                                       :: scaleconc=1.e12
  real, parameter                                       :: smallnum=1.e-38
  logical                                               :: lexist
  integer                                               :: ierr
  integer                                               :: sp_count_i, sp_count_r
  integer                                               :: jjjjmmdd, hhmiss
  integer                                               :: ii, ir, n, jy, ix, kz
  real                                                  :: fact
  integer, dimension(nxgrid*nygrid*nzgrid)              :: sparse_dump_i
  real, dimension(nxgrid*nygrid*nzgrid)                 :: sparse_dump_r

  gridinit(:,:,:)=0.

  inquire ( file=trim(filename),exist=lexist )
  if ( .not.lexist ) then
    write(*,*) 'WARNING: cannot find '//trim(filename)
    go to 10
  endif
  open(100, file=trim(filename), form='unformatted', action='read', status='old', iostat=ierr)
  if ( ierr.ne.0 ) then
    write(*,*) 'WARNING: cannot read '//trim(filename)
    go to 10
  endif

  read(100) jjjjmmdd
  read(100) hhmiss

  fact = 1.
  read(100) sp_count_i
  read(100) (sparse_dump_i(ix), ix=1, sp_count_i)
  read(100) sp_count_r
  read(100) (sparse_dump_r(ix), ix=1, sp_count_r)
  ii = 0
  do ir = 1, sp_count_r
    if ((sparse_dump_r(ir)*fact).gt.smallnum) then
      ii = ii + 1
      n = sparse_dump_i(ii)
      fact = fact*(-1.)
    else
      n = n + 1
    endif
    kz = n/(nxgrid*nygrid)
    jy = (n - kz*nxgrid*nygrid)/nxgrid
    ix = n - nxgrid*nygrid*kz-nxgrid*jy
    gridinit(ix+1,jy+1,kz) = abs(sparse_dump_r(ir))*scaleconc
  end do

  close(100)

10 continue

end subroutine read_init
