subroutine awdeposition(ns, izo, ime, ivert, awdlo)

    !!!!!!!!!!!!!! ADJOINT of wdeposition - ONLY against conc() !!!!!!
    !  Loss terms due to wet deposition

    use chimere_consts        !GD 20-04-10  to be checked
    use worker_common        !GD 20-04-10  to be checked

    implicit none


    ! subroutine arguments
    integer, intent(in) :: ns
    integer, intent(in) :: izo, ime, ivert
    real(kind = 8) :: awdlo,tho
    
    tho=dun/dtr2
    !  For gases

    if(inwetd(ns, 1).ne.0) then
        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                &      min(wetdr1(inwetd(ns, 1), izo, ime, ivert),tho) * awdlo
    endif

    if(inwetd(ns, 2).ne.0) then
        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                &      min(wetdr2(inwetd(ns, 2), izo, ime, ivert),tho) * awdlo
    endif

END subroutine awdeposition
