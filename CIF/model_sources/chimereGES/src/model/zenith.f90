subroutine zenith

!--------------------------------------------------------------------------------------
! Calculates current cosine of zenithal angle 		
! Calculation based on the M.Jacobson book 'Fundamentals of atmospheric modeling' (p.280)
!	
! INPUT:		
!    soltim: number of entire days between the winter solstice and the current day. In hours.
!    tim0: the current hour of the day (from 0 to 23)
!    tim1: the number of hour between the winter solstice and the current hour
!    phi: the latitude of the studied location (converted in radians)
!    ts: the number of seconds past local noon.
!    Njd: the number of days between the winter solstice and the current day. In days.
!    epsob: the obliquity of the ecliptic
!    Lm: mean longitude of the sun (in degrees)
!    gm: mean anomaly of the sun (in degrees)
!    lambdaec: ecliptic longitude of the sun
!    delta: solar declination angle (in radians)
!    Ha: local hour angle of the sun (converted in radians)
! OUTPUT:   
!    ZENILOC  Currrent cosines of zenithal angles
!--------------------------------------------------------------------------------------
  use worker_common
  implicit none
  integer :: izo,ime
  real(kind=8) :: tim0, tim1, zenitmp
  real(kind=8) :: phi, Njd, epsob, Lm, gm, lambdaec, ts, Ha, delta
  real(kind=8), parameter :: pirad=pi/180d0
!--------------------------------------------------------------------------------------
!  Declination at 12Z                                                   

   tim0 = ihour(ihourrun) + thour 
   tim1 = soltim + ihourrun 

   do ime=1,nmerid
   do izo=1,nzonal
     
      phi=xlati(izo,ime)*pirad
      Njd=tim1/24.
      epsob=23.439-4.d-7*Njd
      Lm=280.46+0.9856474*Njd
      gm=357.528+0.9856003*Njd
      lambdaec=Lm+1.915*sin(gm*pirad)+0.02*sin(2*gm*pirad)
      delta=dasin(sin(epsob*pirad)*sin(lambdaec*pirad))
      ts=tim0-12.0+xlong(izo,ime)*24d0/360d0
      if(ts.lt.0)ts=ts+24.
      Ha=2.0*pi*ts/24.
      zenitmp=sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(Ha)
      zeniloc(izo,ime)=zenitmp

   enddo
   enddo

end subroutine zenith
