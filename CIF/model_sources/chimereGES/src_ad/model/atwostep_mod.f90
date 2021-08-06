module atwostep_mod

    !  Adjoint of the TWOSTEP algorithm
    !  One time-step of the TWOSTEP algorithm [Version with fixed time step]
    !  INPUT : CONC
    !          CONCO
    !  OUTPUT: CONC     Array of concentartions
    !          CONCO    Copy of CONC at beginning

    use chimere_consts        !GD 20-04-10 to be check
    use worker_common         !GD 20-04-10 to be check

    implicit none
    include 'mpif.h'

    !*****************************************************************************************
    private
    integer :: ni, ns, nitgsef
    integer :: ispec, izo, ime, ivert, ivex, iex, ild
    real(kind = 8), dimension(:, :, :, :), allocatable :: cconc
    real(kind = 8), dimension(:, :, :, :), allocatable :: conc0
    real(kind = 8), dimension(:, :, :, :), allocatable :: conc1
    real(kind = 8), dimension(:, :, :, :), allocatable :: sink0
    real(kind = 8), dimension(:, :, :, :), allocatable :: acconc
    real(kind = 8) :: ratloss, co, source, sink, useful
    real(kind = 8) :: asource, asink, aco

    integer, dimension(mpi_status_size) :: status
    integer :: ierr
    integer :: send_right_req, send_left_req, recv_right_req, recv_left_req
    integer :: asend_right_req, asend_left_req, arecv_right_req, arecv_left_req

    real(kind = 8), allocatable, dimension(:, :) :: send_line_buf
    real(kind = 8), allocatable, dimension(:, :) :: asend_line_buf
    real(kind = 8), allocatable, dimension(:, :) :: asend_line_buf_quat
    real(kind = 8), allocatable, dimension(:, :) :: asend_line_buf_cinq
    real(kind = 8), allocatable, dimension(:, :) :: recv_line_buf
    real(kind = 8), allocatable, dimension(:, :) :: arecv_line_buf
    real(kind = 8), allocatable, dimension(:, :) :: arecv_line_buf_ter, asend_line_buf_ter
    real(kind = 8), allocatable, dimension(:, :) :: arecv_line_buf_quat
    real(kind = 8), allocatable, dimension(:, :) :: arecv_line_buf_cinq
    real(kind = 8), allocatable, dimension(:, :, :, :) :: send_right_buf, asend_right_buf
    real(kind = 8), allocatable, dimension(:, :, :, :) :: send_left_buf
    real(kind = 8), allocatable, dimension(:, :, :, :) :: recv_right_buf
    real(kind = 8), allocatable, dimension(:, :, :, :) :: recv_left_buf, arecv_left_buf

    public :: atwostep   !Catherine

    !*****************************************************************************************
contains

    subroutine atwostep

        allocate(send_line_buf(nspectot, 3))
        allocate(asend_line_buf(nspectot, 3))
        allocate(asend_line_buf_ter(nspectot, 3))
        allocate(asend_line_buf_quat(nspectot, 3))
        allocate(asend_line_buf_cinq(nspectot, 3))
        allocate(recv_line_buf(nspectot, 3))
        allocate(arecv_line_buf(nspectot, 3))
        allocate(arecv_line_buf_ter(nspectot, 3))
        allocate(arecv_line_buf_quat(nspectot, 3))
        allocate(arecv_line_buf_cinq(nspectot, 3))
        allocate(send_right_buf(nspectot, 3, nmerid, nverti))
        allocate(asend_right_buf(nspectot, 3, 1, 1))
        allocate(send_left_buf(nspectot, 3, nmerid, nverti))
        allocate(recv_right_buf(nspectot, 3, nmerid, nverti))
        allocate(recv_left_buf(nspectot, 3, nmerid, nverti))
        allocate(arecv_left_buf(nspectot, 3, 1, 1))
        allocate(cconc(nspec, -2:nzonal + 3, -2:nmerid + 3, 1:nverti + 1))
        allocate(acconc(nspec, -2:nzonal + 3, -2:nmerid + 3, 1:nverti + 1))    !Catherine
        allocate(conc0(nspec, nzonal + 3, nmerid, nverti))    !Catherine+IP plus large
        allocate(conc1(nspec, nzonal, nmerid, nverti))    !Catherine
        allocate(sink0(nspec, nzonal, nmerid, nverti))    !Catherine


        !print*,'  First stage : direct calculations',rank
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
                        conc0(ispec, izo, ime, ivert) = 2d0 * conc(ispec, izo, ime, ivert) - &
                                conco(ispec, izo, ime, ivert)
                        if(useabsclipconc.eq.1) then
                            if (abs(conc0(ispec, izo, ime, ivert)) < clipconc) &
                                    conc0(ispec, izo, ime, ivert) = clipconc
                        else
                            if (conc0(ispec, izo, ime, ivert) < clipconc) &               ! Catherine
                                    conc0(ispec, izo, ime, ivert) = clipconc                 !GD 20-04-10 change conc to conc0
                        endif
                        conco(ispec, izo, ime, ivert) = conc(ispec, izo, ime, ivert)
                        conc(ispec, izo, ime, ivert) = conc0(ispec, izo, ime, ivert)       !GD 20-04-10 added compared !to seq
                    end do
                end do
            end do
        end do

        !print*,' conc has been modified. But not yet halos ...',rank
        call update_left_and_right_halos

        ! IP ici, conc=conc0 dans les domaines + les halos
        ! IP mais on n'a pas les conc0 des halos en fait!
        ! IP on stocke les conc0 necessaires en plus
        do ivert = 1, nverti
            do ime = 1, nmerid
                do izo = nzonal, nzonal + 3
                    do ispec = 1, nspec
                        conc0(ispec, izo, ime, ivert) = conc(ispec, izo, ime, ivert)
                    enddo
                enddo
            enddo
        enddo

        call mpi_barrier(wrk_comm, ierr)

        !print*,'  Gauss-Seidel iterations',rank
        if(ihourrun.lt.ihoursu) then
            nitgsef = nitgssu
        else
            nitgsef = nitgs
        endif

        !GD 20-04-10 : to be checked, this version of the adjoint works only with nitgsef==1
        if(nitgsef.ne.1) print*, '**** NITGSEF MUST BE =1 IN THIS VERSION ****'
        !   do ni=1,nitgsef

        do ivert = 1, nverti
            do ime = 1, nmerid
                if(dom_i==1) then! This is the first process of the line. Next process on the right needs its results
                    ! new conc
                    do izo = 1, nzonal
                        do ns = 1, nspec
                            call prodloss(ns, izo, ime, ivert, source, sink)
                            sink0(ns, izo, ime, ivert) = sink         !GD 20-04-10
                            useful = dun / (conc0(ns, izo, ime, ivert) + dtr2 * sink)
                            conc(ns, izo, ime, ivert) = conc0(ns, izo, ime, ivert) * useful * (cconc(ns, izo, ime, ivert) + dtr2 * source)
                            conc1(ns, izo, ime, ivert) = conc(ns, izo, ime, ivert) !GD 20-04-10 save conc in conc1
                        end do
                    end do
                    call send_first_end_of_line
                else ! This is not the first process of the line. We have a left halo to update
                    call update_line_left_halo
                    ! new conc
                    do izo = 1, nzonal
                        do ns = 1, nspec
                            call prodloss(ns, izo, ime, ivert, source, sink)
                            sink0(ns, izo, ime, ivert) = sink                !GD 20-04-10
                            useful = dun / (conc0(ns, izo, ime, ivert) + dtr2 * sink)
                            conc(ns, izo, ime, ivert) = conc0(ns, izo, ime, ivert) * useful * (cconc(ns, izo, ime, ivert) + dtr2 * source)
                            conc1(ns, izo, ime, ivert) = conc(ns, izo, ime, ivert)   !GD 20-04-10 save conc in conc1
                        end do
                    end do
                    call send_current_end_of_line
                end if
            end do ! ime=1,nmerid
        end do ! ivert=1,nverti

        call mpi_barrier(wrk_comm, ierr)
        call update_right_halo

        !    end do ! ni=1,nitgsef

        call mpi_barrier(wrk_comm, ierr)

        !***********************************************************************************************!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   ADJOINT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !***********************************************************************************************!
       !GD 20-04-10 valid only for nitgsef=1?
        !print*,'ADJOINT',rank,nspec
        call mpi_barrier(wrk_comm, ierr)

        do ivert = nverti, 1, -1
            do ime = nmerid, 1, -1

                if(dom_i==1) then
                     
                    call asend_first_end_of_line
                    if(nzdoms>1)then
                        ! on met conc a la bonne valeur aussi dans le halo
                        conc(1:nspec, nzonal + 1:nzonal + 3, ime, ivert) = conc0(:, nzonal + 1:nzonal + 3, ime, ivert)
                    endif

                    do izo = nzonal, 1, -1
                        do ns = nspec, 1, -1
                            conc(ns, izo, ime, ivert) = conc0(ns, izo, ime, ivert)
                            useful = dun / (conc0(ns, izo, ime, ivert) + dtr2 * sink0(ns, izo, ime, ivert))
                            asource = useful * conc0(ns, izo, ime, ivert) * dtr2 * aconc(ns, izo, ime, ivert)
                            asink = -useful * dtr2 * aconc(ns, izo, ime, ivert) * conc1(ns, izo, ime, ivert)
                            acconc(ns, izo, ime, ivert) = useful * conc0(ns, izo, ime, ivert) * aconc(ns, izo, ime, ivert)
                            aconc(ns, izo, ime, ivert) = (-useful + dun / conc0(ns, izo, ime, ivert))&
                                    * conc1(ns, izo, ime, ivert) * aconc(ns, izo, ime, ivert)
                            call aprodloss(ns, izo, ime, ivert, asource, asink)

                        enddo
                    enddo

                else ! dom_i>1

                    call asend_current_end_of_line
                    if(nzdoms>1.and.dom_i<nzdoms)then
                        ! on met conc a la bonne valeur aussi dans le halo
                        conc(1:nspec, nzonal + 1:nzonal + 3, ime, ivert) = conc0(:, nzonal + 1:nzonal + 3, ime, ivert)
                    endif

                    do izo = nzonal, 1, -1
                        do ns = nspec, 1, -1
                            conc(ns, izo, ime, ivert) = conc0(ns, izo, ime, ivert)
                            useful = dun / (conc0(ns, izo, ime, ivert) + dtr2 * sink0(ns, izo, ime, ivert))
                            asource = useful * conc0(ns, izo, ime, ivert) * dtr2 * aconc(ns, izo, ime, ivert)
                            asink = -useful * dtr2 * aconc(ns, izo, ime, ivert) * conc1(ns, izo, ime, ivert)
                            acconc(ns, izo, ime, ivert) = useful * conc0(ns, izo, ime, ivert) * aconc(ns, izo, ime, ivert)
                            aconc(ns, izo, ime, ivert) = (-useful + dun / conc0(ns, izo, ime, ivert))&
                                    * conc1(ns, izo, ime, ivert) * aconc(ns, izo, ime, ivert)
                            call aprodloss(ns, izo, ime, ivert, asource, asink)
                        enddo
                    enddo

                    call aupdate_line_left_halo

                endif ! dom_i

            enddo ! ime
        enddo ! ivert

        call mpi_barrier(wrk_comm, ierr)

        call aupdate_right_halo

        do ivert = nverti, 1, -1
            do ime = nmerid, 1, -1
                do izo = nzonal, 1, -1
                    do ispec = nspec, 1, -1
                        if(conc0(ispec, izo, ime, ivert)==clipconc) &
                                aconc(ispec, izo, ime, ivert) = dzero
                        if(cconc(ispec, izo, ime, ivert)==clipconc)&
                                acconc(ispec, izo, ime, ivert) = dzero

                        aco = aconco(ispec, izo, ime, ivert)
                        aconco(ispec, izo, ime, ivert) = -acconc(ispec, izo, ime, ivert) / 3d0 &
                                - aconc(ispec, izo, ime, ivert)
                        aconc(ispec, izo, ime, ivert) = 4d0 * acconc(ispec, izo, ime, ivert) / 3d0 &
                                + 2d0 * aconc(ispec, izo, ime, ivert) + aco

                    enddo
                enddo
            enddo
        enddo
 
        deallocate(send_line_buf)
        deallocate(asend_line_buf)
        deallocate(asend_line_buf_ter)
        deallocate(asend_line_buf_quat)
        deallocate(asend_line_buf_cinq)
        deallocate(recv_line_buf)
        deallocate(arecv_line_buf)
        deallocate(arecv_line_buf_ter)
        deallocate(arecv_line_buf_quat)
        deallocate(arecv_line_buf_cinq)
        deallocate(send_right_buf)
        deallocate(asend_right_buf)
        deallocate(send_left_buf)
        deallocate(recv_right_buf)
        deallocate(recv_left_buf)
        deallocate(arecv_left_buf)
        deallocate(cconc)
        deallocate(acconc)
        deallocate(sink0)   !catherine
        deallocate(conc0)   !catherine
        deallocate(conc1)   !catherine

    END subroutine atwostep


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
            ! ON NE PASSE PAS ICI SI SEULEMENT DEUX TUILES!!!!
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
        ! ON NE PASSE PAS ICI SI SEULEMENT DEUX TUILES
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ADJOINT
    !******************************************************************
    ! routine pour adjoint transmission
    !******************************************************************
    subroutine asend_first_end_of_line

        if(nzdoms>1) then
            call mpi_recv(arecv_line_buf, nspectot * 3, mpi_double_precision, rank + 1, &
                    1000 + 200 * (rank + 1) + ime, wrk_comm, status, ierr)
            aconc(:, nzonal - 2:nzonal, ime, ivert) = aconc(:, nzonal - 2:nzonal, ime, ivert) + arecv_line_buf
        end if

    end subroutine asend_first_end_of_line

    !******************************************************************
    subroutine aupdate_line_left_halo

        asend_line_buf = aconc(:, -2:0, ime, ivert)
        call mpi_send(asend_line_buf, nspectot * 3, mpi_double_precision, rank - 1, &
                1000 + 200 * rank + ime, wrk_comm, ierr)
        aconc(:, -2:0, ime, ivert) = 0.d0

    end subroutine aupdate_line_left_halo

    !******************************************************************
    subroutine aupdate_line_right_halo

        if(dom_i<nzdoms) then

            asend_line_buf_ter = aconc(:, nzonal - 2:nzonal, ime - 1, ivert)
            call mpi_send(asend_line_buf_ter, nspectot * 3, mpi_double_precision, rank + 1, &
                    3000 + 200 * rank + ime, wrk_comm, ierr)
            aconc(:, nzonal - 2:nzonal, ime - 1, ivert) = 0.d0

            asend_line_buf_quat = aconc(:, nzonal + 1:nzonal + 3, ime, ivert)
            call mpi_send(asend_line_buf_quat, nspectot * 3, mpi_double_precision, rank + 1, &
                    4000 + 200 * rank + ime, wrk_comm, ierr)
            aconc(:, nzonal + 1:nzonal + 3, ime, ivert) = 0.d0

        end if

        if(dom_i>1) then

            call mpi_recv(arecv_line_buf_ter, nspectot * 3, mpi_double_precision, rank - 1, &
                    3000 + 200 * (rank - 1) + ime, wrk_comm, status, ierr)
            aconc(:, -2:0, ime - 1, ivert) = arecv_line_buf_ter

            call mpi_recv(arecv_line_buf_quat, nspectot * 3, mpi_double_precision, rank - 1, &
                    4000 + 200 * (rank - 1) + ime, wrk_comm, status, ierr)
            aconc(:, 1:3, ime, ivert) = arecv_line_buf_quat

        endif

    end subroutine aupdate_line_right_halo
    !*******************************************************************
    subroutine aupdate_verti

        if(dom_i<nzdoms) then
            do ns = 1, nspectot
                do izo = 1, 3
                    asend_line_buf_cinq(ns, izo) = aconc(ns, nzonal - 3 + izo, nmerid, ivert - 1)
                enddo
            enddo

            call mpi_send(asend_line_buf_cinq, nspectot * 3, mpi_double_precision, rank + 1, &
                    5000 + 200 * rank, wrk_comm, ierr)
            aconc(:, nzonal - 2:nzonal, nmerid, ivert - 1) = 0.d0

        end if

        if(dom_i>1) then

            call mpi_recv(arecv_line_buf_cinq, nspectot * 3, mpi_double_precision, rank - 1, &
                    5000 + 200 * (rank - 1), wrk_comm, status, ierr)
            do ns = 1, nspectot
                do izo = 1, 3
                    aconc(ns, -3 + izo, nmerid, ivert - 1) = arecv_line_buf_cinq(ns, izo)
                enddo
            enddo

        endif

    end subroutine  aupdate_verti

    !******************************************************************
    subroutine asend_current_end_of_line

        if(dom_i<nzdoms) then
            !   ON NE PASSE PAS ICI SI SEULEMENT DEUX TUILES
            call mpi_recv(arecv_line_buf, nspectot * 3, mpi_double_precision, rank + 1, &
                    1000 + 200 * (rank + 1) + ime, wrk_comm, status, ierr)
            aconc(:, nzonal - 2:nzonal, ime, ivert) = aconc(:, nzonal - 2:nzonal, ime, ivert) + arecv_line_buf
        end if

    end subroutine asend_current_end_of_line

    !******************************************************************
    subroutine aupdate_right_halo

        if(dom_i<nzdoms) then
            send_right_buf = aconc(:, nzonal + 1:nzonal + 3, 1:nmerid, 1:nverti)
            call mpi_issend(&
                    send_right_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank + 1, &
                    600 + rank + 1, &
                    wrk_comm, &
                    asend_right_req, &
                    ierr                            &
                    )
            call mpi_wait(asend_right_req, mpi_status_ignore, ierr)
            aconc(:, nzonal + 1:nzonal + 3, 1:nmerid, 1:nverti) = 0.d0
        end if

        if(dom_i>1) then
            call mpi_irecv(&
                    recv_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank - 1, &
                    600 + rank, &
                    wrk_comm, &
                    arecv_left_req, &
                    ierr                            &
                    )
            call mpi_wait(arecv_left_req, mpi_status_ignore, ierr)
            aconc(:, 1:3, 1:nmerid, 1:nverti) = aconc(:, 1:3, 1:nmerid, 1:nverti) + recv_left_buf
        end if

    end subroutine aupdate_right_halo
    !******************************************************************

    subroutine aupdate_left_halo

        if(dom_i<nzdoms) then
            send_right_buf = aconc(:, nzonal - 2:nzonal, 1:nmerid, 1:nverti)
            call mpi_issend(&
                    send_right_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank + 1, &
                    600 + rank + 1, &
                    wrk_comm, &
                    asend_right_req, &
                    ierr                            &
                    )
            call mpi_wait(asend_right_req, mpi_status_ignore, ierr)
        end if

        if(dom_i>1) then
            call mpi_irecv(&
                    recv_left_buf, &
                    nspectot * 3 * nmerid * nverti, &
                    mpi_double_precision, &
                    rank - 1, &
                    600 + rank, &
                    wrk_comm, &
                    arecv_left_req, &
                    ierr                            &
                    )
            call mpi_wait(arecv_left_req, mpi_status_ignore, ierr)
            aconc(:, -2:0, 1:nmerid, 1:nverti) = recv_left_buf
        end if

    end subroutine aupdate_left_halo
    !******************************************************************

end module atwostep_mod
