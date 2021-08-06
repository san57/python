subroutine locvalues

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

    ! DEEP CONVECTION IMPLEMENTATION BEGIN[1]
    real(kind = 8), parameter :: epsdeepconv = 1.d-4
    integer :: nvcloudtop(nzonal, nmerid) ! Top of the cloud
    integer :: nvdettop(nzonal, nmerid)   ! Start of the downdraft
    ! DEEP CONVECTION IMPLEMENTATION END[1]
    !**********************************************************************************
    ! airmloc has been evaluated in master and sent (hopefully) to worker
    ! before call to this routine.

    !print*,'  Interpolation weights'
    w1 = thour
    w2 = dun - w1
    ih = ihour(ihourrun)
    idt = idtyp(ihourrun)

    !print*,'  Interpolation of physics'
    !print*,'  1: 2D variables'

    hghtloc(:, :) = w2 * hght(:, :, 1) + w1 * hght(:, :, 2)
    atteloc(:, :) = w2 * atte(:, :, 1) + w1 * atte(:, :, 2)
    tem2loc(:, :) = w2 * tem2(:, :, 1) + w1 * tem2(:, :, 2)
    ustaloc(:, :) = w2 * usta(:, :, 1) + w1 * usta(:, :, 2)
    aerrloc(:, :) = w2 * aerr(:, :, 1) + w1 * aerr(:, :, 2)
    obukloc(:, :) = w2 * obuk(:, :, 1) + w1 * obuk(:, :, 2)
    wstaloc(:, :) = w2 * wsta(:, :, 1) + w1 * wsta(:, :, 2)
    topcloc(:, :) = w2 * topc(:, :, 1) + w1 * topc(:, :, 2)
    srehloc(:, :) = w2 * sreh(:, :, 1) + w1 * sreh(:, :, 2)

    !print*,'  2: 3D variables'

    ! Caution : some meteo variables are defined on a domain wider than nmerid*nzonal
    ! lmbb: use LIMITED array assignment

    airmloc(1:nzonal, 1:nmerid, 1:nverti) = w2 * airm(1:nzonal, 1:nmerid, 1:nverti, 1) + w1 * airm(1:nzonal, 1:nmerid, 1:nverti, 2)
    winzloc(1:nzonal, 1:nmerid, 1:nverti) = w2 * winz(1:nzonal, 1:nmerid, 1:nverti, 1) + w1 * winz(1:nzonal, 1:nmerid, 1:nverti, 2)
    winmloc(1:nzonal, 1:nmerid, 1:nverti) = w2 * winm(1:nzonal, 1:nmerid, 1:nverti, 1) + w1 * winm(1:nzonal, 1:nmerid, 1:nverti, 2)
    temploc(:, :, :) = w2 * temp(:, :, :, 1) + w1 * temp(:, :, :, 2)
    sphuloc(:, :, :) = w2 * sphu(:, :, :, 1) + w1 * sphu(:, :, :, 2)
    clwcloc(:, :, :) = w2 * clwc(:, :, :, 1) + w1 * clwc(:, :, :, 2)
    hlayloc(:, :, :) = w2 * hlay(:, :, :, 1) + w1 * hlay(:, :, :, 2)
    kzzzloc(:, :, :) = w2 * kzzz(:, :, :, 1) + w1 * kzzz(:, :, :, 2)

    !print*,' DEEP CONVECTION IMPLEMENTATION BEGIN[2]'
    if(ideepconv.ne.0)then

        ! Updrought and downdraught mass fluxes, entrainment and detrainment fluxes
        ! ideep is to know if we are in a convective column or not
        do ime = 1, nmerid
            do izo = 1, nzonal
                ideep(izo, ime) = 0
                do ivert = 1, nverti
                    if(dpeu(izo, ime, ivert, 1).gt.epsdeepconv.and.dpeu(izo, ime, ivert, 2).gt.epsdeepconv)then
                        ideep(izo, ime) = 1
                        dpeuloc(izo, ime, ivert) = w2 * dpeu(izo, ime, ivert, 1) + w1 * dpeu(izo, ime, ivert, 2)
                        dpedloc(izo, ime, ivert) = w2 * dped(izo, ime, ivert, 1) + w1 * dped(izo, ime, ivert, 2)
                        dpduloc(izo, ime, ivert) = w2 * dpdu(izo, ime, ivert, 1) + w1 * dpdu(izo, ime, ivert, 2)
                        dpddloc(izo, ime, ivert) = w2 * dpdd(izo, ime, ivert, 1) + w1 * dpdd(izo, ime, ivert, 2)
                    else
                        dpeuloc(izo, ime, ivert) = dzero
                        dpedloc(izo, ime, ivert) = dzero
                        dpduloc(izo, ime, ivert) = dzero
                        dpddloc(izo, ime, ivert) = dzero
                    endif
                enddo
            enddo
        enddo

        !  Calculation of the updraught and downdraught mass fluxes
        do ivert = 1, nverti
            do ime = 1, nmerid
                do izo = 1, nzonal
                    if(ideep(izo, ime).eq.1) then
                        if(ivert.eq.1)then
                            flxuloc(izo, ime, ivert) = max(dzero, dpeuloc(izo, ime, ivert) - dpduloc(izo, ime, ivert))
                            flxdloc(izo, ime, ivert) = min(dzero, dpedloc(izo, ime, ivert) - dpddloc(izo, ime, ivert))
                        else
                            flxuloc(izo, ime, ivert) = max(dzero, flxuloc(izo, ime, ivert - 1) + dpeuloc(izo, ime, ivert) - dpduloc(izo, ime, ivert))
                            flxdloc(izo, ime, ivert) = min(dzero, flxdloc(izo, ime, ivert - 1) + dpedloc(izo, ime, ivert) - dpddloc(izo, ime, ivert))
                        endif
                    endif
                enddo
            enddo
        enddo
        !  Calculation of cloud top and detrainement starting level
        do ime = 1, nmerid
            do izo = 1, nzonal
                if(ideep(izo, ime).eq.1) then
                    nvcloudtop(izo, ime) = 1
                    do ivert = nverti, 1, -1
                        if(dpduloc(izo, ime, ivert).gt.dzero) then
                            nvcloudtop(izo, ime) = ivert
                            goto 1005
                        endif
                    enddo
                    1005      continue
                    dpduloc(izo, ime, nvcloudtop(izo, ime)) = dpduloc(izo, ime, nvcloudtop(izo, ime)) &
                            + flxuloc(izo, ime, nvcloudtop(izo, ime))
                    flxuloc(izo, ime, nvcloudtop(izo, ime)) = dzero
                    nvdettop(izo, ime) = 1
                    do ivert = nverti, 1, -1
                        if(dpedloc(izo, ime, ivert).gt.dzero) then
                            nvdettop(izo, ime) = ivert
                            goto 1006
                        endif
                    enddo
                    1006      continue
                    dpedloc(izo, ime, nvdettop(izo, ime)) = dpedloc(izo, ime, nvdettop(izo, ime)) &
                            - flxdloc(izo, ime, nvdettop(izo, ime))
                    flxdloc(izo, ime, nvdettop(izo, ime)) = dzero
                    do ivert = nvcloudtop(izo, ime) + 1, nverti
                        flxuloc(izo, ime, ivert) = dzero
                        dpeuloc(izo, ime, ivert) = dzero
                        dpduloc(izo, ime, ivert) = dzero
                    enddo
                    do ivert = nvdettop(izo, ime) + 1, nverti
                        flxdloc(izo, ime, ivert) = dzero
                        dpedloc(izo, ime, ivert) = dzero
                        dpddloc(izo, ime, ivert) = dzero
                    enddo
                endif ! of ideep=1
            enddo
        enddo
    endif ! of ideepconv=1
    !print*,' DEEP CONVECTION IMPLEMENTATION END[2]'

    !print*,'  Layer thicknesses definitions'
    do ivert = 1, nverti
        do ime = 1, nmerid
            do izo = 1, nzonal
                if (ivert.eq.1) then
                    thlayloc(izo, ime, ivert) = hlayloc(izo, ime, ivert)
                else
                    thlayloc(izo, ime, ivert) = hlayloc(izo, ime, ivert) - hlayloc(izo, ime, ivert - 1)
                end if
            end do
        end do
    end do
    call mpi_barrier(wrk_comm, ierr)

    if(dom_i==1) then
        airmloc(0, :, :) = airmloc(1, :, :)
        airmloc(-1, :, :) = airmloc(1, :, :)
        airmloc(-2, :, :) = airmloc(1, :, :)
        winzloc(0, :, :) = winzloc(1, :, :)
        winmloc(0, :, :) = winmloc(1, :, :)
        thlayloc(0, :, :) = thlayloc(1, :, :)
        if(nzdoms>1) then
            call exchange_halo_right(airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
            call exchange_halo_right(winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_right(winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_right(thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        end if
    end if

    if(dom_i==nzdoms) then
        airmloc(nzonal + 1, :, :) = airmloc (nzonal, :, :)
        airmloc(nzonal + 2, :, :) = airmloc (nzonal, :, :)
        airmloc(nzonal + 3, :, :) = airmloc (nzonal, :, :)
        winzloc(nzonal + 1, :, :) = winzloc (nzonal, :, :)
        winmloc(nzonal + 1, :, :) = winmloc (nzonal, :, :)
        thlayloc(nzonal + 1, :, :) = thlayloc(nzonal, :, :)
        if(nzdoms>1) then
            call exchange_halo_left (airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
            call exchange_halo_left (winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_left (winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_left (thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        end if
    end if

    if((dom_i>1).and.(dom_i<nzdoms)) then
        call exchange_halo_right(airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
        call exchange_halo_right(winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_right(winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_right(thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_left (airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
        call exchange_halo_left (winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_left (winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_left (thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
    end if

    if(dom_j==1) then
        airmloc(:, 0, :) = airmloc (:, 1, :)
        airmloc(:, -1, :) = airmloc (:, 1, :)
        airmloc(:, -2, :) = airmloc (:, 1, :)
        winzloc(:, 0, :) = winzloc (:, 1, :)
        winmloc(:, 0, :) = winmloc (:, 1, :)
        thlayloc(:, 0, :) = thlayloc(:, 1, :)
        if(nmdoms>1) then
            call exchange_halo_upper(airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
            call exchange_halo_upper(winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_upper(winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_upper(thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        end if
    end if

    if(dom_j==nmdoms) then
        airmloc(:, nmerid + 1, :) = airmloc (:, nmerid, :)
        airmloc(:, nmerid + 2, :) = airmloc (:, nmerid, :)
        airmloc(:, nmerid + 3, :) = airmloc (:, nmerid, :)
        winzloc(:, nmerid + 1, :) = winzloc (:, nmerid, :)
        winmloc(:, nmerid + 1, :) = winmloc (:, nmerid, :)
        thlayloc(:, nmerid + 1, :) = thlayloc(:, nmerid, :)
        if(nmdoms>1) then
            call exchange_halo_lower(airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
            call exchange_halo_lower(winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_lower(winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
            call exchange_halo_lower(thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        end if
    end if

    if((dom_j>1).and.(dom_j<nmdoms)) then
        call exchange_halo_upper(airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
        call exchange_halo_upper(winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_upper(winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_upper(thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_lower(airmloc, 3, -2, nzonalmax + 3, -2, nmeridmax + 3, 1, nverti + 1)
        call exchange_halo_lower(winzloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_lower(winmloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
        call exchange_halo_lower(thlayloc, 1, 0, nzonalmax + 1, 0, nmeridmax + 1, 1, nverti)
    end if

    call mpi_barrier(wrk_comm, ierr)
    airmloc(:, :, nverti + 1) = airmloc(:, :, nverti)
    presloc(1:nzonal, 1:nmerid, 1:nverti) = airmloc(1:nzonal, 1:nmerid, 1:nverti) * temploc(1:nzonal, 1:nmerid, 1:nverti) / 7.2868d16
    dtenloc(1:nzonal, 1:nmerid, 1:nverti) = (airm(1:nzonal, 1:nmerid, 1:nverti, 2) - airm(1:nzonal, 1:nmerid, 1:nverti, 1)) / 36d2

    do ivert = 1, nverti
        do ime = 1, nmerid
            do izo = 1, nzonal
                !  Calculate relative humidity : Teten's formula
                clwckgkg = an * clwcloc(izo, ime, ivert) / (airmloc(izo, ime, ivert) * 29d0)
                if(clwckgkg.le.1d-6) then ! 1d-6 in unit of clwckgkg [kg/kg]
                    incloud(izo, ime, ivert) = 0
                else
                    incloud(izo, ime, ivert) = 1
                endif

            end do
        end do
    end do

    !  Interpolation of biogenic emissions
    !  Anthropic emissions are NOT interpolated
    if(tl.eq.1) emisbloc_tl(1:nemisb, 1:nzonal, 1:nmerid) = w2 * emisb_tl(1:nemisb, 1:nzonal, 1:nmerid, 1) + w1 * emisb_tl(1:nemisb, 1:nzonal, 1:nmerid, 2)
    emisbloc(1:nemisb, 1:nzonal, 1:nmerid) = w2 * emisb(1:nemisb, 1:nzonal, 1:nmerid, 1) + w1 * emisb(1:nemisb, 1:nzonal, 1:nmerid, 2)

contains

    !***********************************************************************
    subroutine exchange_halo_right(var, nz, i1, i2, j1, j2, k1, k2)
        integer :: nz, i1, i2, j1, j2, k1, k2
        real(kind = 8), dimension(i1:i2, j1:j2, k1:k2) :: var
        real(kind = 8), allocatable, dimension(:, :, :) :: send_buf, recv_buf
        integer :: send_req, recv_req, ierr

        allocate(send_buf(nz, nmerid, nverti))
        allocate(recv_buf(nz, nmerid, nverti))
        call mpi_irecv(&
                recv_buf, &
                nz * nmerid * nverti, &
                mpi_double_precision, &
                rank + 1, &
                11000 + rank + 1, &
                wrk_comm, &
                recv_req, &
                ierr                                 &
                )
        send_buf = var(nzonal - nz + 1:nzonal, 1:nmerid, 1:nverti)
        call mpi_issend(&
                send_buf, &
                nz * nmerid * nverti, &
                mpi_double_precision, &
                rank + 1, &
                10000 + rank, &
                wrk_comm, &
                send_req, &
                ierr                                 &
                )

        call mpi_wait(recv_req, mpi_status_ignore, ierr)
        var(nzonal + 1:nzonal + nz, 1:nmerid, 1:nverti) = recv_buf
        call mpi_wait(send_req, mpi_status_ignore, ierr)

        deallocate(send_buf)
        deallocate(recv_buf)

    end subroutine exchange_halo_right

    !***********************************************************************
    subroutine exchange_halo_left(var, nz, i1, i2, j1, j2, k1, k2)
        integer :: nz, i1, i2, j1, j2, k1, k2
        real(kind = 8), dimension(i1:i2, j1:j2, k1:k2) :: var
        real(kind = 8), allocatable, dimension(:, :, :) :: send_buf, recv_buf
        integer :: send_req, recv_req, ierr

        allocate(send_buf(nz, nmerid, nverti))
        allocate(recv_buf(nz, nmerid, nverti))
        call mpi_irecv(&
                recv_buf, &
                nz * nmerid * nverti, &
                mpi_double_precision, &
                rank - 1, &
                10000 + rank - 1, &
                wrk_comm, &
                recv_req, &
                ierr                                 &
                )

        send_buf = var(1:nz, 1:nmerid, 1:nverti)
        call mpi_issend(&
                send_buf, &
                nz * nmerid * nverti, &
                mpi_double_precision, &
                rank - 1, &
                11000 + rank, &
                wrk_comm, &
                send_req, &
                ierr                                 &
                )
        call mpi_wait(recv_req, mpi_status_ignore, ierr)
        var(1 - nz:0, 1:nmerid, 1:nverti) = recv_buf
        call mpi_wait(send_req, mpi_status_ignore, ierr)

        deallocate(send_buf)
        deallocate(recv_buf)

    end subroutine exchange_halo_left


    !***********************************************************************
    subroutine exchange_halo_upper(var, nm, i1, i2, j1, j2, k1, k2)
        integer :: nm, i1, i2, j1, j2, k1, k2
        real(kind = 8), dimension(i1:i2, j1:j2, k1:k2) :: var
        real(kind = 8), allocatable, dimension(:, :, :) :: send_buf, recv_buf
        integer :: send_req, recv_req, ierr

        allocate(send_buf(nzonal, nm, nverti))
        allocate(recv_buf(nzonal, nm, nverti))
        call mpi_irecv(&
                recv_buf, &
                nzonal * nm * nverti, &
                mpi_double_precision, &
                rank + nzdoms, &
                11000 + rank + nzdoms, &
                wrk_comm, &
                recv_req, &
                ierr                                 &
                )

        send_buf = var(1:nzonal, nmerid - nm + 1:nmerid, 1:nverti)
        call mpi_issend(&
                send_buf, &
                nzonal * nm * nverti, &
                mpi_double_precision, &
                rank + nzdoms, &
                10000 + rank, &
                wrk_comm, &
                send_req, &
                ierr                                 &
                )
        call mpi_wait(recv_req, mpi_status_ignore, ierr)
        var(1:nzonal, nmerid + 1:nmerid + nm, 1:nverti) = recv_buf
        call mpi_wait(send_req, mpi_status_ignore, ierr)

        deallocate(send_buf)
        deallocate(recv_buf)

    end subroutine exchange_halo_upper


    !***********************************************************************
    subroutine exchange_halo_lower(var, nm, i1, i2, j1, j2, k1, k2)
        integer :: nm, i1, i2, j1, j2, k1, k2
        real(kind = 8), dimension(i1:i2, j1:j2, k1:k2) :: var
        real(kind = 8), allocatable, dimension(:, :, :) :: send_buf, recv_buf
        integer :: send_req, recv_req, ierr

        allocate(send_buf(nzonal, nm, nverti))
        allocate(recv_buf(nzonal, nm, nverti))
        call mpi_irecv(&
                recv_buf, &
                nzonal * nm * nverti, &
                mpi_double_precision, &
                rank - nzdoms, &
                10000 + rank - nzdoms, &
                wrk_comm, &
                recv_req, &
                ierr                                 &
                )

        send_buf = var(1:nzonal, 1:nm, 1:nverti)
        call mpi_issend(&
                send_buf, &
                nzonal * nm * nverti, &
                mpi_double_precision, &
                rank - nzdoms, &
                11000 + rank, &
                wrk_comm, &
                send_req, &
                ierr                                 &
                )
        call mpi_wait(recv_req, mpi_status_ignore, ierr)
        var(1:nzonal, 1 - nm:0, 1:nverti) = recv_buf
        call mpi_wait(send_req, mpi_status_ignore, ierr)

        deallocate(send_buf)
        deallocate(recv_buf)

    end subroutine exchange_halo_lower

end subroutine locvalues
