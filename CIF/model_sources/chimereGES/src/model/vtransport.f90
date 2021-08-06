subroutine vtransport

    !  Calculation of VERTICAL first-order UPWIND flux rates for every model
    !  cell. Using the first-order scheme, the transport terms, in the gover
    !  equation, are of the form:
    !  dCi/dt = SUM ( AijCj ) - BiCi
    !  The AijCj terms represent incoming fluxes, and the sum runs over some
    !  Of the adjacent boxes j of box i. The factors Aij are called the Flux
    !  rates, and are calculated in this subroutine, as well as the outgoing
    !  flux rates Bi.
    !
    !  Vertical fluxes contain vertical TRANSPORT fluxes and vertical MIXING
    !  fluxes. All are processed using this first-order scheme. Vertical TRA
    !  fluxes are DIAGNOSED in order to achieve instantaneous mass balance.
    !  this purpose, vertical winds are estimated cell by cell in an air col
    !  from bottom to top.
    !
    !  Vertical transport+mixing is written in such a way that transilient m
    !  is possible. Hence the matrix VFLUXI(NBOX,NVERTI+1) describes the inc
    !  fluxes along the vertical. For instance VFLUXI(NB,NV) gives the flux
    !  coming from layer NV into box NV (above or below). When NV=NVERTI+1,
    !  is the flux coming from the atmosphere above the model column. In thi
    !  version, however, mixing has a Kz local formulation.
    !
    !  This routine is called every physical time step. It MUST be called AF
    !  a call to routine HTRANSPORT since vertical winds, which are calculat
    !  here, depend on the horizontal fluxes.
    !
    !  INPUT : WINXLOC    Mixing rate between current and above layer (a vel
    !          THLAYLOC   Current thicknesses of the model layers
    !          AIRMLOC    Current density field
    !          AIRM       Hourly density field
    !          FLUXW      Horizontal flux rate at the western boundary
    !          FLUXE      Horizontal flux rate at the eastern boundary
    !          FLUXS      Horizontal flux rate at the southern boundary
    !          FLUXN      Horizontal flux rate at the northern boundary
    !  OUTPUT: WINVLOC    Vertical wind field
    !          VFLUXI     Vertical incoming flux rates
    !          VFLUXO     Vertical outgoing flux rates

    use chimere_consts
    use worker_common

    implicit none

    !*************************************************************************
    integer :: izo, ime, ivert
    real(kind = 8) :: hlow, hupp
    real(kind = 8) :: dratio, dtende
    real(kind = 8) :: belowflux, fluxbal, fluxtmp, add


    !*************************************************************************
    !  Initialization of the number of vertically importing cells (NVIMPORT)
    !  The first importing cell is the upper one. The second is the lower on
    !  One also initializes the outgoing flux rate VFLUXO

    vfluxo = dzero
    vfluxi = dzero

    !  Loop on model columns

    do ime = 1, nmerid
        do izo = 1, nzonal

            !  Loop from bottom to top layers

            do ivert = 1, nverti
                hlow = thlayloc(izo, ime, ivert)
                if(ivert.lt.nverti) hupp = thlayloc(izo, ime, ivert + 1)
                dratio = airmloc(izo, ime, ivert) / airmloc(izo, ime, ivert + 1)
                dtende = dtenloc(izo, ime, ivert)

                !  Calculation of mass flux imbalance due to horizontal advection
                !  and vertical advection from below.

                if(ivert.gt.1) then
                    if(winvloc(izo, ime, ivert - 1).gt.0.)then
                        belowflux = airmloc(izo, ime, ivert - 1) * winvloc(izo, ime, ivert - 1)   &
                                / hlow
                    else
                        belowflux = airmloc(izo, ime, ivert) * winvloc(izo, ime, ivert - 1)   &
                                / hlow
                    endif

                else
                    belowflux = dzero
                endif

                ! fluxbal in [s^-1]
                fluxbal = belowflux       
                add = max(fluxw(izo, ime, ivert), dzero) * airmloc(izo - 1 ,ime, ivert)
                fluxbal = fluxbal + add
                add = max(-fluxe(izo, ime, ivert), dzero) * airmloc(izo + 1 ,ime, ivert) 
                fluxbal = fluxbal + add                             
                add = max(fluxs(izo, ime, ivert), dzero) * airmloc(izo, ime -1, ivert)
                fluxbal = fluxbal + add 
                add = max(-fluxn(izo, ime, ivert), dzero) * airmloc(izo, ime + 1, ivert) 
                fluxbal = fluxbal + add 
                add = max(-fluxw(izo, ime, ivert), dzero) * airmloc(izo, ime, ivert) 
                fluxbal = fluxbal - add 
                add = max(fluxe(izo, ime, ivert), dzero) * airmloc(izo, ime, ivert)
                fluxbal = fluxbal - add 
                add = max(-fluxs(izo, ime, ivert), dzero) * airmloc(izo, ime, ivert)
                fluxbal = fluxbal - add 
                add = max(fluxn(izo, ime, ivert), dzero) * airmloc(izo, ime, ivert)
                fluxbal = fluxbal - add 
                add = dtende
                fluxbal = fluxbal - add    

                !  Calculation of vertical wind and transport flux rates
                if(fluxbal.gt.dzero) then
                    add = hlow / airmloc(izo, ime, ivert)
                    winvloc(izo, ime, ivert) = fluxbal * add
                else
                    add =  hlow / airmloc(izo, ime, ivert + 1)
                    winvloc(izo, ime, ivert) = fluxbal * add
                endif
    
                !  Calculation of mixing fluxes between layer IVERT and IVERT+1

                fluxtmp = winxloc(izo, ime, ivert) / hlow
                vfluxi(izo, ime, ivert, ivert + 1) = vfluxi(izo, ime, ivert, ivert + 1) &
                        + fluxtmp * dratio
                vfluxo(izo, ime, ivert) = vfluxo(izo, ime, ivert) + fluxtmp
                if(ivert.lt.nverti) then
                    fluxtmp = winxloc(izo, ime, ivert) / hupp
                    vfluxi(izo, ime, ivert + 1, ivert) = vfluxi(izo, ime, ivert + 1, ivert) &
                            + fluxtmp
                    vfluxo(izo, ime, ivert + 1) = vfluxo(izo, ime, ivert + 1) + &
                            fluxtmp * dratio
                endif
            enddo
        enddo
    enddo
end subroutine vtransport

