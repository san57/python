module twostep_mod

    !  One time-step of the TWOSTEP algorithm [Version with fixed time step]
    !  INPUT : CONC     Array of concentrations
    !          CONCO    Previous array
    !          CLIPCONC Clipping value to preserve positiveness
    !  OUTPUT: CONC     Array of concentartions
    !          CONCO    Copy of CONC at beginning

    use chimere_consts
    use worker_common
    use worker_message_subs

    implicit none


    !*****************************************************************************************
    private
    integer :: ni, ns, nitgsef
    integer :: ispec, izo, ime, ivert, ivex, iex, ild
    real(kind = 8), dimension(:, :, :, :), allocatable :: cconc
    real(kind = 8), dimension(:, :, :, :), allocatable :: ccons
    real(kind = 8) :: ratloss, co, source, sink

    integer :: ierr
    integer :: send_right_req, send_left_req, recv_right_req, recv_left_req

    real(kind = 8), allocatable, dimension(:, :) :: send_line_buf
    real(kind = 8), allocatable, dimension(:, :) :: recv_line_buf
    real(kind = 8), allocatable, dimension(:, :, :, :) :: send_right_buf
    real(kind = 8), allocatable, dimension(:, :, :, :) :: send_left_buf
    real(kind = 8), allocatable, dimension(:, :, :, :) :: recv_right_buf
    real(kind = 8), allocatable, dimension(:, :, :, :) :: recv_left_buf

    public :: twostep

    !*****************************************************************************************
contains

    subroutine twostep
        allocate(send_line_buf(nspectot, 3))
        allocate(recv_line_buf(nspectot, 3))
        allocate(send_right_buf(nspectot, 3, nmerid, nverti))
        allocate(send_left_buf(nspectot, 3, nmerid, nverti))
        allocate(recv_right_buf(nspectot, 3, nmerid, nverti))
        allocate(recv_left_buf(nspectot, 3, nmerid, nverti))
        allocate(cconc(nspec, -2:nzonal + 3, -2:nmerid + 3, 1:nverti + 1))
        allocate(ccons(nspec, -2:nzonal + 3, -2:nmerid + 3, 1:nverti + 1))

        do ivert = 1, nverti
            do ime = 1, nmerid
                do izo = 1, nzonal
                    do ispec = 1, nspec
                        cconc(ispec, izo, ime, ivert) = &
                                (4d0 * conc(ispec, izo, ime, ivert) - &
                                        conco(ispec, izo, ime, ivert)) / 3d0
                        if(useabsclipconc.eq.1) then
                            if (abs(cconc(ispec, izo, ime, ivert)) < clipconc) &
                                    cconc(ispec, izo, ime, ivert) = clipconc
                        else
                            if (cconc(ispec, izo, ime, ivert) < clipconc) &
                                    cconc(ispec, izo, ime, ivert) = clipconc
                        endif
                        ccons(ispec, izo, ime, ivert) = conc(ispec, izo, ime, ivert)
                        conc(ispec, izo, ime, ivert) = &
                                & 2d0 * conc(ispec, izo, ime, ivert) - &
                                & conco(ispec, izo, ime, ivert)
                        if(useabsclipconc.eq.1) then
                            if (abs(conc(ispec, izo, ime, ivert)) < clipconc) &
                                    conc(ispec, izo, ime, ivert) = clipconc
                        else
                            if (conc(ispec, izo, ime, ivert) < clipconc) &
                                    conc(ispec, izo, ime, ivert) = clipconc
                        endif
                        conco(ispec, izo, ime, ivert) = ccons(ispec, izo, ime, ivert)
                    end do
                end do
            end do
        end do
        
        ! conc has been modified. But not yet halos ...
        call update_left_and_right_halos

        call mpi_barrier(wrk_comm, ierr)

        !  Gauss-Seidel iterations
        if(ihourrun.lt.ihoursu) then
            nitgsef = nitgssu
        else
            nitgsef = nitgs
        endif

        do ni = 1, nitgsef          
            do ivert = 1, nverti
                do ime = 1, nmerid 
                    if(dom_i==1) then
                        ! This is the first process of the line. Next process on the right needs its results
                        ! new conc
                        do izo = 1, nzonal
                            do ns = 1, nspec
                                call prodloss(ns, izo, ime, ivert, source, sink)
                                ratloss = sink / conc(ns, izo, ime, ivert)
                                conc(ns, izo, ime, ivert) = (cconc(ns, izo, ime, ivert) + dtr2 * source) / &
                                        (dun + dtr2 * ratloss) 
                            end do
                        end do
                        call send_first_end_of_line
                    else
                        ! This is not the first process of the line. We have a left halo to update
                        call update_line_left_halo
                        ! new conc
                        do izo = 1, nzonal
                            do ns = 1, nspec
                                call prodloss(ns, izo, ime, ivert, source, sink)
                                ratloss = sink / conc(ns, izo, ime, ivert)
                                conc(ns, izo, ime, ivert) = &
                                        (cconc(ns, izo, ime, ivert) + dtr2 * source) / (dun + dtr2 * ratloss)
                            end do
                        end do
                        call send_current_end_of_line
                    end if
                end do ! ime=1,nmerid
            end do ! ivert=1,nverti

            call mpi_barrier(wrk_comm, ierr)
            call update_right_halo

        end do ! ni=1,nitgsef

        call mpi_barrier(wrk_comm, ierr)

        !  Accumulated variables
        do ivert = 1, nverti
            drydep(:, :, :) = drydep(:, :, :) + drydepi(:, :, :, ivert) * dtr
            wetdep(:, :, :) = wetdep(:, :, :) + wetdepi(:, :, :, ivert) * dtr
        end do

        deallocate(send_line_buf)
        deallocate(recv_line_buf)
        deallocate(send_right_buf)
        deallocate(send_left_buf)
        deallocate(recv_right_buf)
        deallocate(recv_left_buf)
        deallocate(cconc)
        deallocate(ccons)

    END subroutine twostep

    !******************************************************************
    subroutine update_left_and_right_halos

        if((dom_i==1).and.(nzdoms>1)) then
            !print *,rank,'receives new right halo from',rank+1
            call mpi_irecv(&
                    recv_right_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank + 1, &
                    800 + rank + 1, &
                    wrk_comm, &
                    recv_right_req, &
                    ierr                            &
                    )
            !print *,rank,'sends my right data to',rank+1
            send_right_buf = conc(:, nzonal - 2:nzonal, 1:nmerid, 1:nverti)
            call mpi_issend(&
                    send_right_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank + 1, &
                    700 + rank, &
                    wrk_comm, &
                    send_right_req, &
                    ierr                            &
                    )
            call mpi_wait(recv_right_req, mpi_status_ignore, ierr)
            conc(:, nzonal + 1:nzonal + 3, 1:nmerid, 1:nverti) = recv_right_buf
            call mpi_wait(send_right_req, mpi_status_ignore, ierr)
            !print *,'done 1',rank
        end if

        if((dom_i>1).and.(dom_i<nzdoms)) then
            !print *,rank,'receives new right halo from',rank+1
            call mpi_irecv(&
                    recv_right_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank + 1, &
                    800 + rank + 1, &
                    wrk_comm, &
                    recv_right_req, &
                    ierr                            &
                    )
            !print *,rank,'sends my right data to',rank+1
            send_right_buf = conc(:, nzonal - 2:nzonal, 1:nmerid, 1:nverti)
            call mpi_issend(&
                    send_right_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank + 1, &
                    700 + rank, &
                    wrk_comm, &
                    send_right_req, &
                    ierr                            &
                    )
            !print *,rank,'receives new left halo from',rank-1
            call mpi_irecv(&
                    recv_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank - 1, &
                    700 + rank - 1, &
                    wrk_comm, &
                    recv_left_req, &
                    ierr                            &
                    )
            !print *,rank,'sends my left data to',rank-1
            send_left_buf = conc(:, 1:3, 1:nmerid, 1:nverti)
            call mpi_issend(&
                    send_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank - 1, &
                    800 + rank, &
                    wrk_comm, &
                    send_left_req, &
                    ierr                            &
                    )
            !print *,'waiting 2',rank
            call mpi_wait(recv_right_req, mpi_status_ignore, ierr)
            conc(:, nzonal + 1:nzonal + 3, 1:nmerid, 1:nverti) = recv_right_buf
            call mpi_wait(recv_left_req, mpi_status_ignore, ierr)
            conc(:, -2:0, 1:nmerid, 1:nverti) = recv_left_buf
            call mpi_wait(send_right_req, mpi_status_ignore, ierr)
            call mpi_wait(send_left_req, mpi_status_ignore, ierr)
            !print *,'done 2'
        end if

        if((dom_i==nzdoms).and.(nzdoms>1)) then
            !print *,rank,'sends my left data to',rank-1
            send_left_buf = conc(:, 1:3, 1:nmerid, 1:nverti)
            !print *,rank,'receives new left halo from',rank-1
            call mpi_irecv(&
                    recv_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank - 1, &
                    700 + rank - 1, &
                    wrk_comm, &
                    recv_left_req, &
                    ierr                            &
                    )
            call mpi_issend(&
                    send_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank - 1, &
                    800 + rank, &
                    wrk_comm, &
                    send_left_req, &
                    ierr                            &
                    )
            call mpi_wait(recv_left_req, mpi_status_ignore, ierr)
            conc(:, -2:0, 1:nmerid, 1:nverti) = recv_left_buf
            call mpi_wait(send_left_req, mpi_status_ignore, ierr)
            !print *,'done 3',rank
        end if

    end subroutine update_left_and_right_halos

    !******************************************************************
    subroutine send_first_end_of_line

        if(nzdoms>1) then
            send_line_buf = conc(:, nzonal - 2:nzonal, ime, ivert)
            call mpi_send(send_line_buf, nspectot * 3, mpi_double_precision, rank + 1, &
                    1000 + 200 * (rank + 1) + ime, wrk_comm, ierr)
        end if

    end subroutine send_first_end_of_line

    !******************************************************************
    subroutine update_line_left_halo

        ! update left halo for each scan line, once the left line is done
        call mpi_recv(recv_line_buf, nspectot * 3, mpi_double_precision, rank - 1, &
                1000 + 200 * rank + ime, wrk_comm, mpi_status_ignore, ierr)
        conc(:, -2:0, ime, ivert) = recv_line_buf

    end subroutine update_line_left_halo

    !******************************************************************
    subroutine send_current_end_of_line

        if(dom_i<nzdoms) then
            send_line_buf = conc(:, nzonal - 2:nzonal, ime, ivert)
            call mpi_send(send_line_buf, nspectot * 3, mpi_double_precision, rank + 1, &
                    1000 + 200 * (rank + 1) + ime, wrk_comm, ierr)
        end if

    end subroutine send_current_end_of_line

    !******************************************************************
    subroutine update_right_halo

        if(dom_i<nzdoms) then
            !print *,rank,'receives new right halo from',rank+1
            call mpi_irecv(&
                    recv_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank + 1, &
                    600 + rank + 1, &
                    wrk_comm, &
                    recv_left_req, &
                    ierr                            &
                    )
            call mpi_wait(recv_left_req, mpi_status_ignore, ierr)
            conc(:, nzonal + 1:nzonal + 3, 1:nmerid, 1:nverti) = recv_left_buf
        end if

        if(dom_i>1) then
            !print *,rank,'sends my left data to',rank-1
            send_left_buf = conc(:, 1:3, 1:nmerid, 1:nverti)
            call mpi_issend(&
                    send_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank - 1, &
                    600 + rank, &
                    wrk_comm, &
                    send_left_req, &
                    ierr                            &
                    )
            !print *,rank,'send_left_req',send_left_req
            call mpi_wait(send_left_req, mpi_status_ignore, ierr)
        end if

    end subroutine update_right_halo

end module twostep_mod
