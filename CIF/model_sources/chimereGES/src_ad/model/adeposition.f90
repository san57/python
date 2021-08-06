subroutine adeposition(ns, izo, ime, ivert, adelo)

    !  Loss terms due to deposition

    use chimere_consts
    use worker_common

    implicit none

    !*****************************************************************************************
    !subroutine arguments
    integer :: ns            ! Species number
    integer :: izo, ime, ivert ! coordinates of current box
    real(kind = 8) :: adelo          ! Deposition loss of species NS in box NB


    !*****************************************************************************************

    if(ivert.eq.1.and.indepo(ns).ne.0) then
        !   adjoint of delo = conc(ns,izo,ime,ivert)*depoloc(indepo(ns),izo,ime)
        aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) &
                & + adelo * depoloc(indepo(ns), izo, ime)
    endif

END subroutine adeposition
