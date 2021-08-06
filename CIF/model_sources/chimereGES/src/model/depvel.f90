subroutine depvel

    !  Dry deposition module for gases taken mostly from Emberson et al. (2000)
    !  Most details are on the EMEP web site www.emep.int
    !  Some validation and details can also be found on Tuovinen et al. (2003, Atm. Env.)

    use chimere_consts
    use worker_common

    implicit none

    !*****************************************************************************************
    integer :: imon, nd, nl
    integer :: nv, ns
    integer :: izo, ime

    real(kind = 8) :: djul_loc
    real(kind = 8) :: Rin, Rb, Rc, Ra
    real(kind = 8) :: ustar, surftemp
    real(kind = 8) :: ts, vapp, vpdkPa, Rlow, dsnow, fsn, frh
    real(kind = 8) :: cz, beamrad, dirrad, difrad, sunlai, shalai, shapar, sunpar
    real(kind = 8) :: Glisun, Glisha, Gsto, Gnst
    real(kind = 8) :: O3Gext, SO2Rdry, SO2Rwet
    real(kind = 8) :: Rl, Rg, Rs

    real(kind = 8), dimension(nvegtype) :: dlai
    real(kind = 8), dimension(nvegtype) :: dsai
    real(kind = 8), dimension(nvegtype) :: fphen
    real(kind = 8), dimension(nvegtype) :: ftem
    real(kind = 8), dimension(nvegtype) :: flig
    real(kind = 8), dimension(nvegtype) :: fvpd
    real(kind = 8), dimension(nvegtype) :: fswp
    real(kind = 8), dimension(nvegtype) :: O3Gsto
    real(kind = 8), dimension(nvegtype) :: O3Gnst
    real(kind = 8), dimension(nvegtype) :: SO2Gnst
    real(kind = 8), dimension(nvegtype) :: NH3Gnst


    !*************************************************************************************
    !  Month

    imon = imonth(ihourrun)

    !  Loop on Horizontal cells

    do ime = 1, nmerid
        do izo = 1, nzonal

            !  Initialization

            do nd = 1, ndepo
                depoloc(nd, izo, ime) = dzero
            enddo


            !  Meteorological parameters

            ts = tem2loc(izo, ime) - t0k
            vapp = 611d0 * exp(17.27d0 * ts / (tem2loc(izo, ime) - 35.86d0))
            vpdkPa = 1d-3 * vapp * (dun - srehloc(izo, ime))

            !  Cold temperature resistance Rlow and presence of snow (0/1)

            Rlow = 1d1 * exp(-ts - 4d0)
            dsnow = dzero

            !  Ustar and Aerodynamic resistance

            ustar = ustaloc(izo, ime)
            Ra = aerrloc(izo, ime)

            if(zeniloc(izo, ime).gt.1d-3) then
                cz = zeniloc(izo, ime)
                beamrad = 4.57d0 * atteloc(izo, ime) * 1069d0 * exp(-0.21d0 / zeniloc(izo, ime))
                dirrad = beamrad * zeniloc(izo, ime)
                difrad = 0.13d0 * dirrad !DEBUG from 2017 version? previously *beamrad
            else
                cz = 1d-3
                dirrad = dzero
                difrad = dzero
            endif

            !  Stomatal resistance as a function of vegetation cover

            do nv = 1, nvegtype
                
                !  Temperature factor

                surftemp = max(deptmin(nv) + 0.1d0, ts)
                surftemp = min(surftemp, 2d0 * deptopt(nv) - deptmin(nv) - 0.1d0)
                ftem(nv) = dun - ((surftemp - deptopt(nv)) / (deptopt(nv) - deptmin(nv)))**2

                !  Vapor Pressure Deficit factor

                if(vpdkPa.le.depvpd1(nv)) then
                    fvpd(nv) = dun
                elseif(vpdkPa.ge.depvpd2(nv)) then
                    fvpd(nv) = 1d-1
                else
                    fvpd(nv) = dun - (depvpd1(nv) - vpdkPa) * (dun - 1d-1) &
                            & / (depvpd1(nv) - depvpd2(nv))
                endif

                !  Soil water pressure

                fswp(nv) = dun


                !  LAI, SAI and phenology

                dlai(nv) = dzero
                dsai(nv) = dzero
                fphen(nv) = dzero
                ! DEBUG local from 2017
                djul_loc=djul
                if(xlati(izo,ime).lt.0.) then
                  djul_loc=mod((djul-int(365./2.)),365d0)    
                  if(djul_loc.eq.0.) djul_loc=1d0
                endif
                if(djul_loc.ge.depsgs(nv).and.djul_loc.le.depegs(nv)) then
                    if(djul.lt.depsgs(nv) + depsgl(nv)) then
                        dlai(nv) = deplai1(nv) + (djul_loc - depsgs(nv)) &
                                & * (deplai2(nv) - deplai1(nv)) / depsgl(nv)
                        dsai(nv) = 5d0 * dlai(nv) / 3.5d0
                    elseif(djul_loc.gt.depegs(nv) - depegl(nv)) then
                        dlai(nv) = deplai1(nv) + (depegs(nv) - djul_loc) &
                                & * (deplai2(nv) - deplai1(nv)) / depegl(nv)
                        dsai(nv) = dlai(nv) + 1.5d0
                    else
                        dlai(nv) = deplai2(nv)
                        dsai(nv) = dlai(nv) + 1.5d0
                    endif
                    if(djul_loc.lt.depsgs(nv) + depphe1(nv)) then
                        fphen(nv) = depphe0(nv) + (djul_loc - depsgs(nv)) * (dun - depphe0(nv)) / depphe1(nv)
                    elseif(djul_loc.gt.depegs(nv) - depphe2(nv)) then
                        fphen(nv) = depphe0(nv) + (depegs(nv) - djul_loc) * (dun - depphe0(nv)) / depphe2(nv)
                    else
                        fphen(nv) = dun
                    endif
                endif
                if(nv.le.4) then
                    dsai(nv) = dlai(nv) + 1.
                elseif(nv.ge.8) then
                    dsai(nv) = dlai(nv)
                endif


                !  Light factor

                sunlai = (dun - exp(-5d-1 * dlai(nv) / cz)) * 2d0 * cz ! here 2d0=1/cos(60deg),
                                              !the constant inclination of leaves in the canopy 
                shalai = dlai(nv) - sunlai
                shapar = difrad * exp(-5d-1 * dlai(nv)**0.7) &
                        & + 0.07d0 * dirrad * (1.1d0 - 0.1d0 * dlai(nv)) * exp(-cz)
                sunpar = 5d-1 * dirrad / cz + shapar! here 5d-1 = cos(60deg),
                                              !the constant inclination of leaves in the canopy 
                Glisun = dun - exp(-depalph(nv) * sunpar)
                Glisha = dun - exp(-depalph(nv) * shapar)
                flig(nv) = (Glisun * sunlai + Glisha * shalai) / max(dlai(nv), 1d-20)

                !  Final O3 stomatal conductance in cm/s

                O3Gsto(nv) = dlai(nv) * gmax(nv) * fphen(nv) * flig(nv) &
                        & * max(fmin(nv), ftem(nv) * fvpd(nv) * fswp(nv)) / 410d0

                !  External plant conductance

                O3Gext = dsai(nv) / 25d0

                !  In-canopy resistance for ozone

                Rin = 1d2 * 14d0 * dsai(nv) * zcanopy(nv) / ustar

                !  Total Ozone non-stomatal conductance

                O3Gnst(nv) = O3Gext + dun / (1d-2 * (Rin + RGSO3(nv)) + Rlow + 2d1 * dsnow)

              enddo ! loop on nv=1,nvegtype

              !----------------------------------------------------------------
              !  SO2 non-stomatal resistance
              do nv=1,nvegtype

                if(acraloc(izo, ime).ge.2d0) then
                    fsn = dun
                else
                    fsn = exp(acraloc(izo, ime) - 2d0)
                endif
                if(srehloc(izo, ime).gt.so2rh(nv)) then
                    frh = (srehloc(izo, ime) - so2rh(nv)) / (dun - so2rh(nv))
                else
                    frh = dzero
                endif
                SO2Rdry = 1.8d0 * fsn + Rlow + 2d1 * dsnow
                SO2Rwet = 1.0d0 * fsn + Rlow + 2d1 * dsnow

                if(nv.le.10) then
                    SO2Gnst(nv) = (dun - frh) / SO2Rdry + frh / SO2Rwet
                else
                    SO2Gnst(nv) = 1d2 / RGSSO2(nv)
                endif

             enddo ! loop on nv=1,nvegtype

              !----------------------------------------------------------------
                !  NH3 non-stomatal conductance
              do nv=1,nvegtype
                NH3Gnst(nv)=5.d-1
                if(nv.le.10) then
                    if(ts.le.-5d0) then
                        NH3Gnst(nv) = 1d-1
                    else if(ts.le.0) then
                        NH3Gnst(nv) = 5d-1
                    else
                        NH3Gnst(nv) = 1d2 / (0.0455d0 &
                                & * (ts + 2d0) * exp((dun - srehloc(izo, ime)) / 0.07) &
                                & * 10d0**(-1.1099 * acraloc(izo, ime) + 1.6769))
                        NH3Gnst(nv) = max(min(1d1, NH3Gnst(nv)), 5d-1)
                    endif
                endif
            enddo ! loop on nv=1,nvegtype

 !----------------------------------------------------------------
! Final loop over landuse, deposited species and vegetation types           
              do nl = 1, nlduse
                do nd = 1, ndepo

                    Rb = factRb(nd) / ustar

                    do nv = 1, nvegtype

                        !  Surface conductances and total resistance

                        Gsto = O3Gsto(nv) / (factD(nd) + Rm(nd) * O3Gsto(nv))
                        Gnst = 1d-5 * dHx(nd) * SO2Gnst(nv) + df0(nd) * O3Gnst(nv)
                        if(inNH3.ne.0.and.indepo(inNH3).eq.nd) Gnst = NH3Gnst(nv)
                        
                        Rc = dun / (Gsto + Gnst)

                        if(inHNO3.ne.0.and.indepo(inHNO3).eq.nd) Rc = max(1d-2, Rlow)
                        

                        !  Calculation of velocity deposition

                        depoloc(nd, izo, ime) = depoloc(nd, izo, ime) &
                                + fveg(izo, ime, nv, nl) * dland(izo, ime, nl) &
                                / (Ra + Rb + Rc) / hlayloc(izo, ime, 1)
                    end do
                end do ! loop over ndepo

            end do ! loop over nlduse
        end do ! loop over nzonal
    end do ! loop over nmerid


end subroutine depvel
