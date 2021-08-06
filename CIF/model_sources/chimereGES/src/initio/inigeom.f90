subroutine inigeom

  !  Initialization of the model geometry
  !  INPUT : XLATI    Latitudes  of the cell centers
  !          XLONG    Longitudes of the cell centers
  !          XLATIC   Latitudes  of the cell corners
  !          XLONGC   Longitudes of the cell corners
  !  OUTPUT: CLATI    Cosines of latitudes
  !          SLATI    Sines of latitudes
  !          XSIZE    Zonal extent of cells
  !          YSIZE    Meridional extent of cells
  !          XBASX    Normalized vector in the dir. of cell west  side: x-coord
  !          XBASY    Normalized vector in the dir. of cell west  side: y-coord
  !          YBASX    Normalized vector in the dir. of cell south side: x-coord
  !          YBASY    Normalized vector in the dir. of cell south side: y-coord

  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'inigeom.f90')
  
 !*****************************************************************************************
  integer :: izo,ime
  integer :: ncstat,ifncoorn,ifnpol
  real xlongc(nzonal_domain+1,nmerid_domain+1),xlatic(nzonal_domain+1,nmerid_domain+1),pol(nzonal_domain,nmerid_domain)
  real area

  real(kind=8) :: dxx,dxy,dyx,dyy
  real(kind=8) :: dxr1,dxr2,dxx1,dxx2,dxy1,dxy2
  real(kind=8) :: dyr1,dyr2,dyx1,dyx2,dyy1,dyy2
  !*****************************************************************************************
  !print*,'  Reading horizontal coordinates of corners'
  ifncoorn=11
  open(ifncoorn,file=fniCOOcorn)
  do ime=1,nmerid_domain+1
    do izo=1,nzonal_domain+1
      read(ifncoorn,*)xlongc(izo,ime),xlatic(izo,ime),area
    enddo
  enddo
  close(ifncoorn)

  !print*,'polar?',polar
  if (polar.eq.1) then
  !print*,'  Reading polar indicator'
  ifnpol=12
  open(ifnpol,file=fnipol)
  endif
  do ime=1,nmerid_domain
    do izo=1,nzonal_domain
      if (polar.eq.1) then
        read(ifnpol,*)pol(izo,ime)
      else
        pol(izo,ime)=0.
      endif
    enddo
  enddo
  if (polar.eq.1) close(ifnpol)

  
  !print*,'  Calculation of the size of the cells (cm) for transport'
  !  and coordinates of grid mesh to be used 
  !  for the calculation of zenithal angles
  do ime=1,nmerid_domain
    do izo=1,nzonal_domain
        clati(izo,ime) = cos(pi*xlati(izo,ime)/180d0)
        slati(izo,ime) = sin(pi*xlati(izo,ime)/180d0)
    enddo
  enddo


  do ime=1,nmerid_domain
    do izo=1,nzonal_domain
        if(pol(izo,ime).eq.0.) then
            dxx = (modulo(xlongc(izo+1,ime)           &
                   -xlongc(izo,ime)+180d0,360d0)-180d0)  &
                    * pi*earthr*clati(izo,ime)/180d0
            dxy = (xlatic(izo+1,ime)-xlatic(izo,ime)) &
                    * pi*earthr/180d0
            dyx = (modulo((xlongc(izo,ime+1)          &
                   -xlongc(izo,ime))+180d0,360d0)-180d0) &
                    * pi*earthr*clati(izo,ime)/180d0
            dyy = (xlatic(izo,ime+1)-xlatic(izo,ime)) &
                    * pi*earthr/180d0
        else        
            
            ! If polar grid cell, specific computations   
            dxr1 = (90-xlatic(izo,ime)) &
                    * pi*earthr/180d0
            dxx1 = dxr1 * sin(xlongc(izo,ime)*pi/180d0)
            dxy1 = - dxr1 * cos(xlongc(izo,ime)*pi/180d0)
            dxr2 = (90-xlatic(izo+1,ime)) &
                    * pi*earthr/180d0
            dxx2 = dxr2 * sin(xlongc(izo+1,ime)*pi/180d0)
            dxy2 = - dxr2 * cos(xlongc(izo+1,ime)*pi/180d0)
            dxx = dxx2 - dxx1
            dxy = dxy2 - dxy1   
            
            dyr1 = (90-xlatic(izo,ime)) &
                    * pi*earthr/180d0
            dyx1 = dyr1 * sin(xlongc(izo,ime)*pi/180d0)
            dyy1 = - dyr1 * cos(xlongc(izo,ime)*pi/180d0)
            dyr2 = (90-xlatic(izo,ime+1)) &
                    * pi*earthr/180d0
            dyx2 = dyr2 * sin(xlongc(izo,ime+1)*pi/180d0)
            dyy2 = - dyr2 * cos(xlongc(izo,ime+1)*pi/180d0)
            dyx = dyx2 - dyx1
            dyy = dyy2 - dyy1
               
        endif 
        xsize(izo,ime) = sqrt(dxx*dxx+dxy*dxy)
        xbasx(izo,ime) = dxx/xsize(izo,ime)
        xbasy(izo,ime) = dxy/xsize(izo,ime)
        
        ysize(izo,ime) = sqrt(dyx*dyx+dyy*dyy)
        ybasx(izo,ime) = dyx/ysize(izo,ime)
        ybasy(izo,ime) = dyy/ysize(izo,ime)
        
    enddo
  enddo

  !print*,' Boundaries',nzonal_domain,nmerid_domain
  xsize(0,:) = xsize(1,:)
  xsize(-1,:) = xsize(1,:)
  xbasx(0,:) = xbasx(1,:)
  xbasy(0,:) = xbasy(1,:)
  
  xsize(nzonal_domain+1,:) = xsize(nzonal_domain,:) 
  xsize(nzonal_domain+2,:) = xsize(nzonal_domain,:)
  xbasx(nzonal_domain+1,:) = xbasx(nzonal_domain,:)
  xbasy(nzonal_domain+1,:) = xbasy(nzonal_domain,:)

  xsize(:,0) = xsize(:,1) 
  xsize(:,-1) = xsize(:,1)
  xbasx(:,0) = xbasx(:,1)
  xbasy(:,0) = xbasy(:,1)
  
  xsize(:,nmerid_domain+1) = xsize(:,nmerid_domain) 
  xsize(:,nmerid_domain+2) = xsize(:,nmerid_domain)
  xbasx(:,nmerid_domain+1) = xbasx(:,nmerid_domain)
  xbasy(:,nmerid_domain+1) = xbasy(:,nmerid_domain)


  ysize(0,:) = ysize(1,:) 
  ysize(-1,:) = ysize(1,:)
  ybasx(0,:) = ybasx(1,:)
  ybasy(0,:) = ybasy(1,:)
  
  ysize(nzonal_domain+1,:) = ysize(nzonal_domain,:) 
  ysize(nzonal_domain+2,:) = ysize(nzonal_domain,:)
  ybasx(nzonal_domain+1,:) = ybasx(nzonal_domain,:)
  ybasy(nzonal_domain+1,:) = ybasy(nzonal_domain,:)

  ysize(:,0) = ysize(:,1) 
  ysize(:,-1) = ysize(:,1)
  ybasx(:,0) = ybasx(:,1)
  ybasy(:,0) = ybasy(:,1)
  
  ysize(:,nmerid_domain+1) = ysize(:,nmerid_domain) 
  ysize(:,nmerid_domain+2) = ysize(:,nmerid_domain)
  ybasx(:,nmerid_domain+1) = ybasx(:,nmerid_domain)
  ybasy(:,nmerid_domain+1) = ybasy(:,nmerid_domain)


end subroutine inigeom
