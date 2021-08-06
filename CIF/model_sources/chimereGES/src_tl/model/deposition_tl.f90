subroutine deposition_tl(ns,izo,ime,ivert,delo) 

  !  Loss terms due to deposition                                         

  use chimere_consts
  use worker_common

  implicit none

  !*****************************************************************************************
  !subroutine arguments
  integer      :: ns            ! Species number
  integer      :: izo,ime,ivert ! coordinates of current box
  real(kind=8) :: delo          ! Deposition loss of species NS in box NB
  
  
  !*****************************************************************************************
  delo = dzero 

  if(ivert.eq.1.and.indepo(ns).ne.0) then
      delo = conc_tl(ns,izo,ime,ivert)*depoloc(indepo(ns),izo,ime)
  endif
  drydepi_tl(ns,izo,ime,ivert) = delo*thlayloc(izo,ime,ivert)

END subroutine deposition_tl
