function vanleer_tl(cupup,cup,cdo,dx,u,dt,cupup_tl,cup_tl,cdo_tl) 

  implicit none 

  !****************************************************************
  ! function arguments
  real(kind=8) :: cupup,cupup_tl
  real(kind=8) :: cup,cup_tl
  real(kind=8) :: cdo,cdo_tl
  real(kind=8) :: dx
  real(kind=8) :: u
  real(kind=8) :: dt

  ! local variables
  real(kind=8) :: cfl,vanleer_tl,minmod

  !*****************************************************************
  cfl = u*dt/dx 
  if(cfl.gt.1.) then 
    print *,'*** WARNING: CFL > 1 IN VANLEER scheme :',cfl 
  endif

  !  Nonmonotonic cases                         

  if(cup.lt.min(cupup,cdo)) then 
    vanleer_tl = cup_tl 
  elseif(cup.gt.max(cupup,cdo)) then 
    vanleer_tl = cup_tl 
  else 

    !  Monotonic case                             

    if(abs(cdo-cup).lt.abs(cup-cupup)) then 
      minmod = cdo_tl - cup_tl 
    else 
      minmod = cup_tl - cupup_tl 
    endif
    vanleer_tl = cup_tl + 0.5*max(1d0-cfl,0d0)*minmod 
  endif

END function vanleer_tl


function vanleer_nonunif_tl(cupup,cup,cdo,dxupup,dx,dxdo,u,dt,cupup_tl,cup_tl,cdo_tl) 


  implicit none 

  !****************************************************************
  ! function arguments
  real(kind=8) :: cupup,cupup_tl
  real(kind=8) :: cup,cup_tl
  real(kind=8) :: cdo,cdo_tl
  real(kind=8) :: dx,dxupup,dxdo
  real(kind=8) :: u
  real(kind=8) :: dt

  ! local variables
  real(kind=8) :: cfl,vanleer_nonunif_tl,minmod


  !*****************************************************************

  !  Nonmonotonic cases                         

  if(cup.lt.min(cupup,cdo)) then 
    vanleer_nonunif_tl = cup_tl 
  elseif(cup.gt.max(cupup,cdo)) then 
    vanleer_nonunif_tl = cup_tl 
  else 

    !  Monotonic case                             

    if(abs((cdo-cup)/(dx+dxdo)).lt.abs((cup-cupup)/(dx+dxupup))) then 
      cfl = 2.*u*dt/(dxdo+dx) 
      if(cfl.gt.1.) then 
!       print *,'*** WARNING: CFL > 1 IN VANLEER scheme :',cfl
        vanleer_nonunif_tl=cup_tl
      else
        vanleer_nonunif_tl=(cup_tl*dxdo+cdo_tl*dx)/(dx+dxdo)-cfl/2.*(cdo_tl-cup_tl)
      endif
    else 
      cfl=2.*u*dt/(dxupup+dx)
      if(cfl.gt.1.) then 
!       print *,'*** WARNING: CFL > 1 IN VANLEER scheme :',cfl
        vanleer_nonunif_tl=cup_tl
      else
        vanleer_nonunif_tl=cup_tl+dx*(cup_tl-cupup_tl)/(dx+dxupup)-cfl/2.*(cup_tl-cupup_tl) 
      endif
    endif
  endif

END function vanleer_nonunif_tl
