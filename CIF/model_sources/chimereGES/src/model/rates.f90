subroutine rates

    !  Calculation of reaction rate coefficients (the K's and the J's)
    !  Reaction rates are classified into 20 groups depending on their
    !  analytical dependence to temperature, density, etc...
    !  Some groups are free, the user can define other groups.
    !  *** There is no particular logic in the classification! It is
    !      just as logic as ... chemistry itself.
    !  This routine also calculates addresses for T-dependent stoichiometry
    !
    !  INPUT :  ITYPERATE   Type of rate for each reaction
    !           TEMPLOC     Current temperature
    !           AIRMLOC     Current density
    !           SPHULOC     Current H2O
    !           ATTELOC     Current radiation attenuation
    !           IHORIZ      Horizontal address
    !           TABRATE     Array containing rate constants
    !           DTABRATE    Slopes for photolytic rates
    !           ZENILOC     Current cosines of zenithal angles
    !
    !  OUTPUT:  RATE        Array containing the K's and the J's

    use chimere_consts
    use worker_common

    implicit none

    !*****************************************************************************************
    real(kind = 8), parameter :: Rg = 8.205d-2
    real(kind = 8), parameter :: conv = 1d+3
    real(kind = 8), parameter :: mm_HNO3 = 6.3d+1

    integer :: i, ns, nr, nt
    integer :: izo, ime, ivert
    integer :: isdepono2
    integer :: ity
    real(kind = 8) :: te, ai, hu
    real(kind = 8) :: c1, c2, c3, c4, f1, f2, f3, f4, ex, factor, t
    real(kind = 8) :: wl, kc, ve
    real(kind = 8) :: cw, gama
    real(kind = 8) :: rhoa, diff, He, ft
    real(kind = 8) :: tho, za, diam
    real(kind = 8) :: ph, ph1

    real(kind = 8), dimension(ntemps) :: dtemp

    !*****************************************************************************************

    tho = dun / dtr2
    diff = 0.25d-4

    !  First finds NO2 deposition index for HONO formation reaction
    !  after Aumont et al 2002

    do ns = 1, nspec
        if(species(ns)%name.eq.'NO2') then
            isdepono2 = indepo(ns)
        endif
    enddo

    !  Calculation of rates

    do nr = 1, nreac
        ity = ityperate(nr)

        do ivert = 1, nverti
            do ime = 1, nmerid
                do izo = 1, nzonal

                    rate(nr, izo, ime, ivert) = dzero

                    !  Physical parameters

                    te = temploc(izo, ime, ivert)
                    ai = airmloc(izo, ime, ivert)
                    hu = sphuloc(izo, ime, ivert)
                    cw = clwcloc(izo, ime, ivert)
                    ph = phloc(izo, ime, ivert)

                    !  Constant rates

                    if(ity.eq.1) then
                        rate(nr, izo, ime, ivert) = tabrate(1, nr)

                        !  Arrhenius simplified rates

                    else if(ity.eq.2) then
                        rate(nr, izo, ime, ivert) = tabrate(1, nr) * exp(-tabrate(2, nr) / te)

                        !  Arrhenius complete

                    else if(ity.eq.3) then
                        rate(nr, izo, ime, ivert) = tabrate(1, nr) * exp(-tabrate(2, nr) / te)        &
                                * (300d0 / te)**tabrate(3, nr)

                        !  Troe or Fall-off

                    else if(ity.eq.4) then
                        c1 = tabrate(1, nr) * exp(-tabrate(2, nr) / te)                 &
                                * (300d0 / te)**tabrate(3, nr)
                        c2 = tabrate(4, nr) * exp(-tabrate(5, nr) / te)                 &
                                * (300d0 / te)**tabrate(6, nr)
                        c3 = ai * c1
                        c4 = c3 / c2
                        ex = dun / (dun + log10(c4)**2)
                        rate(nr, izo, ime, ivert) = c1 * tabrate(7, nr)**ex / (dun + c4)

                        !  Special case of HONO formation from deposition of NO2

                    else if(ity.eq.15) then
                        if(ivert.eq.1) then
                            rate(nr, izo, ime, ivert) = tabrate(1, nr) * depoloc(isdepono2, izo, ime)
                        else
                            rate(nr, izo, ime, ivert) = dzero
                        endif

                        !  Modified Troe for NO2+OH

                    else if(ity.eq.14) then
                        c1 = tabrate(1, nr) * exp(-tabrate(2, nr) / te)                 &
                                * (300d0 / te)**tabrate(3, nr)
                        c2 = tabrate(4, nr) * exp(-tabrate(5, nr) / te)                 &
                                * (300d0 / te)**tabrate(6, nr)
                        c3 = ai * c1
                        c4 = c3 / c2
                        ex = dun / (dun + ((log10(c4) - 0.12d0) / 1.2d0)**2)
                        rate(nr, izo, ime, ivert) = c1 * tabrate(7, nr)**ex / (dun + c4)

                        !  Photolytic rates

                    else if(ity.eq.5) then
                        rate(nr, izo, ime, ivert) = phrate(iphoto(nr), izo, ime, ivert)

                    else if(ity.eq.25) then
                        rate(nr, izo, ime, ivert) = tabrate(1, nr) / te

                        !  Special case of ozone photolysis

                    else if(ity.eq.13) then
                        factor = hu / (hu + ai * (0.02909d0 * exp(70d0 / te) + 0.06545d0 * exp(110d0 / te)))
                        rate(nr, izo, ime, ivert) = phrate(iphoto(nr), izo, ime, ivert) * factor

                        !  Special types of rate forms

                    else if(ity.eq.6) then
                        f1 = tabrate(1, nr) * exp(-tabrate(2, nr) / te)
                        f2 = tabrate(3, nr) * exp(-tabrate(4, nr) / te)
                        rate(nr, izo, ime, ivert) = f1 * f2 / (dun + f2)
                    else if(ity.eq.7) then
                        f1 = tabrate(1, nr) * exp(-tabrate(2, nr) / te)
                        f2 = tabrate(3, nr) * exp(-tabrate(4, nr) / te)
                        rate(nr, izo, ime, ivert) = f1 / (dun + f2)
                    else if(ity.eq.8) then
                        f1 = tabrate(1, nr) * exp(-tabrate(2, nr) / te)
                        f2 = tabrate(3, nr) * exp(-tabrate(4, nr) / te)
                        f3 = tabrate(5, nr) * exp(-tabrate(6, nr) / te)
                        f4 = tabrate(7, nr) * exp(-tabrate(8, nr) / te)
                        rate(nr, izo, ime, ivert) = 2d0 * sqrt(f1 * f2 * f3 * f4 / ((dun + f3) * (dun + f4)))
                    else if(ity.eq.9) then
                        f1 = tabrate(1, nr) * exp(-tabrate(2, nr) / te)
                        f2 = tabrate(3, nr) * exp(-tabrate(4, nr) / te)
                        f3 = tabrate(5, nr) * exp(-tabrate(6, nr) / te)
                        f4 = tabrate(7, nr) * exp(-tabrate(8, nr) / te)
                        f3 = f3 / (dun + f3)
                        f4 = f4 / (dun + f4)
                        rate(nr, izo, ime, ivert) = 2d0 * sqrt(f1 * f2) * (dun - sqrt(f3 * f4)) * (dun - f4) / (2d0 - f3 - f4)
                    else if(ity.eq.16) then
                        ft = exp(tabrate(3, nr) * (1. / te - 1. / 298.))
                        if(tabrate(4, nr).eq.0) then
                            c1 = 10d0**(-ph)
                            c1 = c1**(tabrate(1, nr))
                        else
                            c1 = 1d0 + 13d0 * (10d0**(-ph))
                        endif

                        rate(nr, izo, ime, ivert) = 1.d-20
                        if(cw.gt.1.d-11) then
                            rate(nr, izo, ime, ivert) = (1.d3 * cw * ((Rg * te)**2) / an) * ft * tabrate(2, nr) &
                                    * tabrate(5, nr) * tabrate(6, nr) / c1
                        endif
                    else if(ity.eq.23) then
                        ft = exp(tabrate(4, nr) * (1.d0 / te - 1.d0 / 298.d0))
                        ph1 = min(ph, 5.0)
                        c1 = 10.d0**(-ph1)
                        rate(nr, izo, ime, ivert) = ft * cw * tabrate(1, nr) * tabrate(2, nr) * &
                                tabrate(3, nr) * 8.314d0 * te * 1.d6 / (1.013d5 * conv)  &
                                / (c1**(tabrate(5, nr)))
                    else if(ity.eq.24) then
                        ft = exp(tabrate(4, nr) * (1.d0 / te - 1.d0 / 298.d0))
                        ph1 = min(ph, 6.0)
                        c1 = 10.d0**(-ph1)
                        rate(nr, izo, ime, ivert) = ft * cw * tabrate(1, nr) * tabrate(2, nr) * &
                                tabrate(3, nr) * 8.314d0 * te * 1.d6 / (1.013d5 * conv)  &
                                / (c1**(tabrate(5, nr)))
                    else if(ity.eq.20) then
                        rate(nr, izo, ime, ivert) = 1.d-10
                        ve = sqrt(8d+3 * 8.314d0 * te / (pi * tabrate(2, nr)))
                        if(tabrate(4, nr).eq.0.) then
                            gama = tabrate(1, nr)
                        else
                            za = (log(1.d0) - log(tabrate(1, nr))) / (log(273.d0) - log(298.d0))
                            gama = za * (log(te) - log(298.d0)) + log(tabrate(1, nr))
                            gama = exp(gama)
                            if(gama.gt.1.) gama = 1.
                            if(gama.lt.0.001) gama = 0.001
                        endif
                        if(tabrate(5, nr).ne.0.) then ! particles
                            diam = 10.d-6 ! cloud droplet diameter
                            c1 = diam / diff / 2.d0
                            c2 = 4.d0 / (gama * ve)
                            rate(nr, izo, ime, ivert) = 1.d-6 * (1.d+8 * cw / (diam * 1.d2)) / (c1 + c2)
                        endif
                        if(rate(nr, izo, ime, ivert).gt.tho) rate(nr, izo, ime, ivert) = tho
                    else if(ity.eq.21) then
                        rhoa = ai * 29.d0 * 1000.d0 / an
                        if(tabrate(1, nr).eq.1) then
                            if(cw.gt.1.d-11) then
                                wl = an * cw / (ai * 29.d0)
                            else
                                wl = dzero
                            endif
                            ve = sqrt(8d+3 * 8.314d0 * te / (pi * tabrate(3, nr)))
                            kc = tabrate(2, nr) / 2.d0 / diff + 4.d0 / ve / tabrate(4, nr)
                            rate(nr, izo, ime, ivert) = 6.d0 * wl * rhoa / (1.d+3 * tabrate(2, nr) * kc)
                        else
                            print *, '*** RATES: ERROR, ITY=21 and C1 != 1'
                            stop
                        endif
                        if(rate(nr, izo, ime, ivert).gt.tho) rate(nr, izo, ime, ivert) = tho
                    else if(ity.eq.22) then
                        He = tabrate(2, nr) * exp(-tabrate(3, nr) * (1.d0 / te - 1.d0 / 298.d0))
                        ve = sqrt(8d+3 * 8.314d0 * te / (pi * tabrate(4, nr)))
                        kc = tabrate(1, nr) / 2.d0 / diff + 4.d0 / ve / tabrate(5, nr) !
                        rate(nr, izo, ime, ivert) = 6.d2 / (8.314d0 * He * te * tabrate(1, nr) * kc)
                        if(rate(nr, izo, ime, ivert).gt.tho) rate(nr, izo, ime, ivert) = tho
                    else
                        print *, '*** ERROR: Reaction ', nr, ': Undefined rate type:', ity
                        stop
                    end if
                end do
            end do
        end do
    end do


    !  T-dependent Stoichiometry

    do nt = 1, ntemps - 1
        dtemp(nt) = tabtemp(nt + 1) - tabtemp(nt)
    enddo

    do ivert = 1, nverti
        do ime = 1, nmerid
            do izo = 1, nzonal
                t = temploc(izo, ime, ivert)
                if(t.lt.tabtemp(1)) then
                    istoit(izo, ime, ivert) = 1
                    wgstl(izo, ime, ivert) = dun
                    wgsth(izo, ime, ivert) = dzero
                elseif(t.lt.tabtemp(ntemps)) then
                    do nt = 1, ntemps - 1
                        if(t.ge.tabtemp(nt).and.t.lt.tabtemp(nt + 1)) then
                            istoit(izo, ime, ivert) = nt
                            wgstl(izo, ime, ivert) = (tabtemp(nt + 1) - t) / dtemp(nt)
                            wgsth(izo, ime, ivert) = dun - wgstl(izo, ime, ivert)
                        end if
                    end do
                else
                    istoit(izo, ime, ivert) = ntemps - 1
                    wgstl(izo, ime, ivert) = dzero
                    wgsth(izo, ime, ivert) = dun
                end if
            end do
        end do
    end do

end subroutine rates
