subroutine inidepo

  !  This routine reads parameters for deposition and prepares intermediary
  !  parameters used in subroutine DEPVEL

  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

  !*****************************************************************************************
  ! constants
  real(kind=8),parameter :: rNu = 0.15d0      ! kinetmatic viscosity of air
  real(kind=8),parameter :: DH2O = 0.25d0     ! molecular diffusivity of pollutant
  ! local variables
  real(kind=8),dimension(nspec) :: dMx
  real,dimension(1000) :: cread
  real(kind=8) :: rMx,rHx,rf0,rRwat,DH2O_Dx
  real(kind=8) :: rKoa0,rAoa,rToa,m1,m2,m3,n1,n2,n3,k1,k2
  real(kind=8) :: rKh0,rTh,rAh,rksea
  real(kind=8) :: rkoc,rkdoc,rdiffair,rdiffwat,rksoil
  character(len=15) :: charspec,dumchar
  integer :: i,j,ns,nd,nl,nv
  integer :: izo,ime
  integer :: itywet
  integer :: ifndepoespe,ifnlanduse,ifndepopars,ifnwetd
  
   character*15 cspec
   character*100000 usedval,notusedval

  !*****************************************************************************************

  !  Initialisations

  indepo = 0


  !  Reading deposited species and related parameters

  notusedval=''
  call opfi(ifndepoespe,fndepoespe,'f','o',0)
  ndepo = 0
  do nd=1,100000
     read(ifndepoespe,*,end=1001)charspec,rMx,rHx,rf0,rRwat
     call findspec(charspec,ns)
     if(ns.gt.0) then
        ndepo = ndepo + 1
        indepo(ns) = ndepo
        dMx(ndepo) = real(rMx)
        dHx(ndepo) = rHx
        df0(ndepo) = rf0
     else
	notusedval=notusedval(1:len_trim(notusedval))  &
	&      //charspec(1:len_trim(charspec))//'; '
     endif
  enddo
1001 continue
  close(ifndepoespe)
  print*,'-------------------------------------------'
  print *,' o DEPO_SPEC ignored species: '
  if(len_trim(notusedval).gt.1) print*,'   ',notusedval(1:len_trim(notusedval)-1)

  
  !  Quasi-laminary boundary layer resistance Rb factor
  do nd=1,ndepo
     DH2O_Dx=sqrt(dMx(nd)/18d0)
     factRb(nd) = ((rNu/(DH2O*prandtl))              &
          &              * DH2O_Dx)**(2d0/3d0)       &
          &              * 2d0/vkarm
  enddo

  !  Mesophylle resistance Rm factor

  do nd=1,ndepo
     Rm(nd)=1d-2/((dHx(nd)/3d3)+(1d2*df0(nd)))
  enddo

  !  Stomatal resistance Rstom factors

  do nd=1,ndepo
     factD(nd)=sqrt(dMx(nd)/48d0)
  enddo

  !  EMEP uses a different landuse classification than the CHIMERE 9-categ
  !  one. The next loops provide the
  !  correspondence factors between 9-class LU and EMEP deposition LU
  !  This should be improved later by the use of a LU inventory adapted to
  !  deposition module


  fveg=dzero
  do ime=1,nmerid_domain
    do izo=1,nzonal_domain
        fveg(izo,ime,5,1) = 7d-1
        fveg(izo,ime,6,1) = 3d-1
        fveg(izo,ime,9,2) = dun
        fveg(izo,ime,13,3)= dun
        fveg(izo,ime,14,4)= dun
        fveg(izo,ime,16,5)= dun
        fveg(izo,ime,10,6)= dun
        fveg(izo,ime,1,7) = dun
        fveg(izo,ime,2,8) = dun
        fveg(izo,ime,14,9)= dun
    end do
  end do

  call opfi(ifndepopars,fndepopars,'f','o',0)
  read(ifndepopars,*,END=1010)dumchar
  read(ifndepopars,*)(gmax(nv),nv=1,nvegtype)
  read(ifndepopars,*)(fmin(nv),nv=1,nvegtype)
  read(ifndepopars,*)(deptmin(nv),nv=1,nvegtype)
  read(ifndepopars,*)(deptopt(nv),nv=1,nvegtype)
  read(ifndepopars,*)(deptmax(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depalph(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depvpd1(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depvpd2(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depsgs(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depegs(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depsgl(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depegl(nv),nv=1,nvegtype)
  read(ifndepopars,*)(deplai1(nv),nv=1,nvegtype)
  read(ifndepopars,*)(deplai2(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depphe0(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depphe1(nv),nv=1,nvegtype)
  read(ifndepopars,*)(depphe2(nv),nv=1,nvegtype)
  read(ifndepopars,*)(zcanopy(nv),nv=1,nvegtype)
  read(ifndepopars,*)(RGSO3(nv),nv=1,nvegtype)
  read(ifndepopars,*)(RGSSO2(nv),nv=1,nvegtype)
  read(ifndepopars,*)(so2rh(nv),nv=1,nvegtype)
1010 continue
  close(ifndepopars)

!  print*,'  Reading landuse fractions'

  call opfi(ifnlanduse,fnlanduse,'f','o',0)
  do ime=1,nmerid_domain
    do izo=1,nzonal_domain
        read(ifnlanduse,*) (dland(izo,ime,nl),nl=1,nlduse)
    end do
  end do
  close(ifnlanduse)
  !ccccccccccccccccccccccccc
!  print*,'  WET SCAVENGING OF GAS c'
  !ccccccccccccccccccccccccc

  !  Initialisations

  do ns=1,nspec
     inwetd(ns,1) = 0
     inwetd(ns,2) = 0
  enddo
  nwetd(1) = 0
  nwetd(2) = 0

  call opfi(ifnwetd,fnwetd,'f','o',0)

  notusedval=''
  do i=1,1000000
     read(ifnwetd,*,end=5765)charspec,itywet,(cread(j),j=1,3)
     call findspec(charspec,ns)
     if(ns.gt.0) then
        nwetd(itywet) = nwetd(itywet) + 1
        inwetd(ns,itywet) = nwetd(itywet)
     else
	notusedval=notusedval(1:len_trim(notusedval))  &
	&      //charspec(1:len_trim(charspec))//'; '
     endif
  enddo
5765 continue
  close(ifnwetd)
  print *,' o WET DEPOSITION ignored species: ',notusedval(1:len_trim(notusedval)-1)


end subroutine inidepo
