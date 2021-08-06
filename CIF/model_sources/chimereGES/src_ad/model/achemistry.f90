subroutine achemistry(ns, izo, ime, ivert, achpr, achlo)
    !!!!!!!!!!!  ADJOINT of chemistry subroutine !!!!!!!!!!!!!!!!!!!!!

    use chimere_consts
    use worker_common

    implicit none

    include 'mpif.h'


    !*****************************************************************************************
    ! subroutine arguments
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: achpr, achlo
    ! local variables
    integer :: kr, it, ir
    integer :: nr
    real(kind = 8) :: trat, atrat

    !*****************************************************************************************

    !  Production Terms

    do kr = 1, kreacp(ns)
        ir = ireacp(ns, kr)
        trat = rate(ir, izo, ime, ivert)
        atrat = achpr * (wgstl(izo, ime, ivert) * stoi(ns, ir, istoit(izo, ime, ivert))           &
                & + wgsth(izo, ime, ivert) * stoi(ns, ir, istoit(izo, ime, ivert) + 1))
        do nr = 1, nreactants(ir)
            trat = trat * conc(irctt(ir, nr), izo, ime, ivert)
        enddo
        do nr = 1, nreactants(ir)
            if(conc(irctt(ir, nr), izo, ime, ivert).ne.dzero)then
                aconc(irctt(ir, nr), izo, ime, ivert) = aconc(irctt(ir, nr), izo, ime, ivert) + &
                        atrat * trat / conc(irctt(ir, nr), izo, ime, ivert)
            endif
        enddo
    enddo

    !  Loss Terms

    do kr = 1, kreacl(ns)
        ir = ireacl(ns, kr)
        trat = rate(ir, izo, ime, ivert)
        atrat = achlo
        do it = 1, nreactants(ir)
            trat = trat * conc(irctt(ir, it), izo, ime, ivert)
        enddo
        do it = 1, nreactants(ir)
            if(conc(irctt(ir, it), izo, ime, ivert) .ne.dzero)then
                aconc(irctt(ir, it), izo, ime, ivert) = aconc(irctt(ir, it), izo, ime, ivert) + &
                        atrat * trat / conc(irctt(ir, it), izo, ime, ivert)
            endif
        enddo
    enddo

    return
END subroutine achemistry
