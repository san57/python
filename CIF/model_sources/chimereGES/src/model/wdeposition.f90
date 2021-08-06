subroutine wdeposition(ns, izo, ime, ivert, wdlo)

    !  Loss terms due to wet deposition

    use chimere_consts
    use worker_common

    implicit none


    ! subroutine arguments
    integer, intent(in) :: ns
    integer, intent(in) :: izo, ime, ivert
    real(kind = 8) :: wdlo,tho
    
    tho=dun/dtr2
    wdlo = dzero

    !  For gases

    if(inwetd(ns, 1).ne.0) then
!        wdlo = conc(ns, izo, ime, ivert) * wetdr1(inwetd(ns, 1), izo, ime, ivert)
! DEBUG version 2017?
      wdlo = conc(ns, izo, ime, ivert) * min(wetdr1(inwetd(ns, 1), izo, ime, ivert), tho)
    endif

    if(inwetd(ns, 2).ne.0) then
!        wdlo = wdlo + conc(ns, izo, ime, ivert) * wetdr2(inwetd(ns, 2), izo, ime, ivert)
! DEBUG version 2017?
      wdlo = wdlo + conc(ns, izo, ime, ivert) * min(wetdr2(inwetd(ns, 2), izo, ime, ivert), tho)  
    endif

    wetdepi(ns, izo, ime, ivert) = wdlo * thlayloc(izo, ime, ivert)

END subroutine wdeposition
