subroutine wdeposition_tl(ns,izo,ime,ivert,wdlo)

  !  Loss terms due to wet deposition

  use chimere_consts
  use worker_common

  implicit none


  ! subroutine arguments
  integer,intent(in) :: ns
  integer,intent(in) :: izo,ime,ivert
  real(kind=8) :: wdlo,tho


  tho = dun / dtr2
  wdlo = dzero

  !  For gases

  if(inwetd(ns,1).ne.0) then
     wdlo = conc_tl(ns,izo,ime,ivert)*min(wetdr1(inwetd(ns,1),izo,ime,ivert),tho)
  endif

  if(inwetd(ns,2).ne.0) then
     wdlo = wdlo + conc_tl(ns,izo,ime,ivert)*min(wetdr2(inwetd(ns,2),izo,ime,ivert),tho)
  endif

  wetdepi_tl(ns,izo,ime,ivert) = wdlo*thlayloc(izo,ime,ivert)

END subroutine wdeposition_tl
