subroutine alocvalues

    use chimere_consts
    use worker_common

    implicit none

    include 'mpif.h'

    !**********************************************************************************
    real(kind = 8) :: w1, w2, rml, psat, clwckgkg, kresult, kresult_max
    integer :: i, ih, idt, ip, ik, ins
    integer :: ne, ns
    integer :: izo, ime, ivert, iemisb
    integer :: izo_max, ime_max, ivert_max, ns_max
    integer :: ierr

    !**********************************************************************************
    ! airmloc has been evaluated in master and sent (hopefully) to worker
    ! before call to this routine.

    !  Interpolation weights

    w1 = thour
    w2 = dun - w1
    ih = ihour(ihourrun)
    idt = idtyp(ihourrun)

    !  Interpolation of biogenic emissions
    !  Anthropic emissions are NOT interpolated
    !***************************
    ! adjoint of emisb
    ! emisbloc(1:nemisb,1:nzonal,1:nmerid) = w2*emisb(1:nemisb,1:nzonal,1:nmerid,1) + w1*emisb(1:nemisb,1:nzonal,1:nmerid,2)
    aemisb(ihourrun, 1:nemisb, 1:nzonal, 1:nmerid) = aemisb(ihourrun, 1:nemisb, 1:nzonal, 1:nmerid) + w2 * aemisbloc(1:nemisb, 1:nzonal, 1:nmerid)
    aemisb(ihourrun + 1, 1:nemisb, 1:nzonal, 1:nmerid) = aemisb(ihourrun + 1, 1:nemisb, 1:nzonal, 1:nmerid) + w1 * aemisbloc(1:nemisb, 1:nzonal, 1:nmerid)

end subroutine alocvalues
