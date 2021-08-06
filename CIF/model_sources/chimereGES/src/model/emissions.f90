subroutine emissions(ns, izo, ime, ivert, empr)

    !  Production terms due to emissions
    !  INPUT:  NS       Species number
    !  OUTPUT: EMPR     Emission production of species NS in current cell izo,ime,ivert

    use worker_common

    implicit none


    !******************************************************************************
    ! subroutine arguments
    integer :: ns
    integer :: izo, ime, ivert
    real(kind = 8) :: empr

    !******************************************************************************
    !  Initialization

    empr = 1d-17

    !  Surface emissions

    if(ivert.eq.1) then

        !  Anthropic emissions
        if(inemisa(ns).ne.0) then
            empr = empr + emisaloc(inemisa(ns), izo, ime, ivert)                     &
                    & / thlayloc(izo, ime, ivert)
        endif

        !  Biogenic emissions

        if(inemisb(ns).ne.0) then
            empr = empr + emisbloc(inemisb(ns), izo, ime)                         &
                    & / thlayloc(izo, ime, ivert)
        endif 

    endif


    !  Altitude emissions

    if(ivert.gt.1.and.ivert.lt.nlevemis + 1) then
        if(inemisa(ns).ne.0) then
            empr = empr + emisaloc(inemisa(ns), izo, ime, ivert)             &
                    & / thlayloc(izo, ime, ivert)
        endif
    endif

    return
END subroutine emissions
