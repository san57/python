subroutine photorates

    use chimere_consts
    use worker_common

    implicit none

    !*****************************************************************************
    integer :: nt, np, k
    integer :: iz
    integer :: k1, k2
    integer :: izo, ime, ivert
    real(kind = 8) :: z1, z2
    real(kind = 8) :: cz, at, ph1, ph2
    real(kind = 8) :: wz1, wz2

    real(kind = 8), dimension(:), allocatable :: wvav

    !*****************************************************************************
    allocate(wvav(nlevphot))

    !  Photolysis rates interpolation in vertical and zenith angle

    do ime = 1, nmerid
        do izo = 1, nzonal
            !  Interpolation of zenith angle
            cz = zeniloc(izo, ime)
            if(cz.le.0) then
                iz = 1
                wz2 = 0.
                wz1 = 1.
            else
                do nt = 1, ntabuzen
                    if(cz.gt.zetaref(nt).and.cz.le.zetaref(nt + 1)) then
                        iz = nt
                        wz2 = (cz - zetaref(iz)) / (zetaref(iz + 1) - zetaref(iz))
                        wz1 = dun - wz2
                    end if
                end do
            end if
            do ivert = 1, nverti
                !  Weights for vertical averages
                z2 = hlayloc(izo, ime, ivert)
                if(ivert.eq.1) then
                    z1 = 0.
                else
                    z1 = hlayloc(izo, ime, ivert - 1)
                end if
                call vertav(nlevphot, altiphot, z1, z2, k1, k2, wvav)
                do np = 1, nphot
                    at = atteloc(izo, ime)
                    ph1 = 0.
                    ph2 = 0.
                    do k = k1, k2
                        ph1 = ph1 + wvav(k) * photoj(iz, k, np)
                        ph2 = ph2 + wvav(k) * photoj(iz + 1, k, np)
                    end do
                    phrate(np, izo, ime, ivert) = at * (wz1 * ph1 + wz2 * ph2)
                enddo
            enddo
        enddo
    enddo
    deallocate(wvav)

end subroutine photorates
