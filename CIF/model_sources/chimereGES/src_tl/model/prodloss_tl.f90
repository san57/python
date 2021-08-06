
subroutine prodloss_tl(ns,izo,ime,ivert,topr,tolo) 

  !  Calculates the sum of production and loss terms for the current cell of (global)
  !  coordinates izo,ime,ivert

  use worker_common

  implicit none
  

  !*************************************************************************
  ! subroutine arguments
  integer      :: ns
  integer      :: izo,ime,ivert
  real(kind=8) :: topr
  real(kind=8) :: tolo

  ! local variables
  real(kind=8) :: chlo,delo,trlo,wdlo
  real(kind=8) :: chpr,empr,trpr


  !*************************************************************************
  chpr=0d0
  chlo=0d0
  delo=0d0
  empr=0d0
  trpr=0d0
  trlo=0d0
  wdlo=0d0
  if(usechemistry.ne.0) call chemistry_tl  (ns,izo,ime,ivert,chpr,chlo) 
  if(usedepos.ne.0)     call deposition_tl (ns,izo,ime,ivert,delo) 
  if(useemissions.ne.0) call emissions_tl  (ns,izo,ime,ivert,empr) 
  if(usetransmix.ne.0)  call transmix_tl   (ns,izo,ime,ivert,trpr,trlo) 
  if(usewetdepos.ne.0)  call wdeposition_tl(ns,izo,ime,ivert,wdlo)
  topr = empr + trpr + chpr 
  tolo = trlo + chlo + delo + wdlo
 
END subroutine prodloss_tl
