subroutine aemissions(ns, izo, ime, ivert, aempr)

    !  Production terms due to emissions
    !  INPUT:  NS       Species number
    !  OUTPUT: EMPR     Emission production of species NS in current cell izo,ime,ivert

    use worker_common

    implicit none


    !******************************************************************************
    ! subroutine arguments
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: aempr

    !******************************************************************************

    !  Surface emissions

    if(ivert.eq.1) then

        !  Anthropic emissions

        if(inemisa(ns).ne.0) then
            !      adjoint of empr = empr + emisaloc(inemisa(ns),izo,ime,ivert)/thlayloc(izo,ime,ivert)
            aemisaloc(ihourrun, inemisa(ns), izo, ime, ivert) = aemisaloc(ihourrun, inemisa(ns), izo, ime, ivert) &
                    & + aempr / thlayloc(izo, ime, ivert)
        endif

        !  Biogenic emissions

        if(inemisb(ns).ne.0) then
            !      adjoint of empr = empr + emisbloc(inemisb(ns),izo,ime)/thlayloc(izo,ime,ivert)
            aemisbloc(inemisb(ns), izo, ime) = aemisbloc(inemisb(ns), izo, ime) &
                    & + aempr / thlayloc(izo, ime, ivert)
        endif
    endif

    !  Altitude emissions

    if(ivert.gt.1.and.ivert.lt.nlevemis + 1) then
        if(inemisa(ns).ne.0) then
            !      adjoint of empr = empr + emisaloc(inemisa(ns),izo,ime,ivert)/thlayloc(izo,ime,ivert)
            aemisaloc(ihourrun, inemisa(ns), izo, ime, ivert) = aemisaloc(ihourrun, inemisa(ns), izo, ime, ivert) &
                    & + aempr / thlayloc(izo, ime, ivert)
        endif
    endif

    return
END subroutine aemissions
