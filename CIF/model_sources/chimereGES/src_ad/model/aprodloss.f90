subroutine aprodloss(ns, izo, ime, ivert, atopr, atolo)

    !!!!!!!!!!!  ADJOINT of prodloss subroutine - ONLY for gases !!!!!!!!!!!!!!!!!!!!!
    !  Calculates the sum of production and loss terms for the current cell of (global)
    !  coordinates izo,ime,ivert

    use worker_common                        !GD 20-04-10  to be checked

    implicit none


    !*************************************************************************
    ! subroutine arguments
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: atopr
    real(kind = 8) :: atolo

    !*************************************************************************
    if(usewetdepos.ne.0)  call awdeposition(ns, izo, ime, ivert, atolo)
    if(usetransmix.ne.0)  call atransmix   (ns, izo, ime, ivert, atopr, atolo)
    if(useemissions.ne.0) call aemissions  (ns, izo, ime, ivert, atopr)
    if(usedepos.ne.0)     call adeposition (ns, izo, ime, ivert, atolo)
    if(usechemistry.ne.0) call achemistry  (ns, izo, ime, ivert, atopr, atolo)

END subroutine aprodloss
