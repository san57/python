module worker_message_subs

  use worker_common
  use message_defs
  implicit none
  include 'mpif.h'


contains

  !*****************************************
  subroutine worker_recv_once
    implicit none
    integer :: ierr
    call mpi_recv(int_params, 1, int_params_typ, 0, ias_int_params, mpi_comm_world,mpi_status_ignore,ierr)
    call rcopy_int_params
    ! Now that we know dimensions, we can allocate to receive arrays
    call worker_allocall
    if (tl.eq.1) call worker_allocall_tl
    if (ad.eq.1) call worker_allocall_ad
    call mpi_recv(dbl_params, 1, dbl_params_typ, 0, ias_dbl_params, mpi_comm_world,mpi_status_ignore,ierr)
    call rcopy_dbl_params
    call recv_int_arrays
    call recv_char_arrays
    call recv_real_arrays
    call recv_dbl_arrays
    call recv_hourly_real_arrays
    if(tl.eq.1) call recv_hourly_real_arrays_tl
    call mpi_recv(nphour,1,mpi_integer,            0,ias_nphour,mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(dtr,   1,mpi_double_precision,   0,ias_dtr,   mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(dtr2,  1,mpi_double_precision,   0,ias_dtr2,  mpi_comm_world,mpi_status_ignore,ierr)
    call recv_output_species
  end subroutine worker_recv_once

  !*****************************************
  subroutine aworker_recv_once

    implicit none

    integer :: ierr
    real(kind=8),dimension(:,:,:,:,:),allocatable :: rbuf5
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf4

! equivalent de

    allocate(aemisb(0:nhourrun+1,nemisb,nzonal,nmerid))
    allocate(aemisaloc(0:nhourrun+1,nemisa,nzonal,nmerid,nlevemis))

    allocate(rbuf5(nhourrun+2,nemisa,nzonal,nmerid,nlevemis))
    call mpi_recv(rbuf5,(nhourrun+2)*nemisa*nzonal*nmerid*nlevemis&
    , mpi_double_precision,0, ias_aemisaloc, mpi_comm_world,mpi_status_ignore,ierr)
    aemisaloc(0:nhourrun+1,:,:,:,:) = rbuf5(:,:,:,:,:)
    deallocate(rbuf5)

    allocate(rbuf4(nhourrun+2,nemisb,nzonal,nmerid))
    call mpi_recv(rbuf4,(nhourrun+2)*nemisb*nzonal*nmerid,mpi_double_precision,0,&
            ias_aemisb,mpi_comm_world,mpi_status_ignore,ierr)
    aemisb(0:nhourrun+1,:,:,:) = rbuf4(:,:,:,:)
    deallocate(rbuf4)

  end subroutine aworker_recv_once

 !*****************************************
  subroutine aworker_send_once

    implicit none

    integer :: ierr
    real(kind=8),dimension(:,:,:,:,:),allocatable :: rbuf5
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf4

    allocate(rbuf5(nhourrun+2,nemisa,nzonal,nmerid,nlevemis))
    rbuf5(1:nhourrun+2,:,:,:,:)=aemisaloc(0:nhourrun+1,1:nemisa,1:nzonal,1:nmerid,1:nlevemis)
    call mpi_send(rbuf5,(nhourrun+2)*nemisa*nzonal*nmerid*nlevemis&
    , mpi_double_precision,0, iar_aemisaloc, mpi_comm_world,ierr)
    deallocate(rbuf5)

    allocate(rbuf4(nhourrun+2,nemisb,nzonal,nmerid))
    rbuf4(1:nhourrun+2,:,:,:)=aemisb(0:nhourrun+1,:,:,:)
    call mpi_send(rbuf4,(nhourrun+2)*nemisb*nzonal*nmerid,mpi_double_precision,0,&
            iar_aemisb,mpi_comm_world,ierr)
    deallocate(rbuf4)

  end subroutine aworker_send_once
!*****************************************

  subroutine worker_recv_hourly
    implicit none
    integer :: ierr

    call recv_hourly_real_arrays
    if (tl.eq.1) call recv_hourly_real_arrays_tl
    ! lm add reception of nphour,dtr,dtr2
    call mpi_recv(nphour,1,mpi_integer,            0,ias_nphour,mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(dtr,   1,mpi_double_precision,   0,ias_dtr,   mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(dtr2,  1,mpi_double_precision,   0,ias_dtr2,  mpi_comm_world,mpi_status_ignore,ierr)

  end subroutine worker_recv_hourly


  !*****************************************
  subroutine worker_recv_frac_hourly

    call recv_frac_hourly_int_arrays
    ! conc excluded

  end subroutine worker_recv_frac_hourly

  !*****************************************
  subroutine worker_send_conc(ns)

    integer :: ns

    integer :: ierr
    real(kind=8),allocatable,dimension(:,:,:) :: dbuf3

    call mpi_barrier(mpi_comm_world,ierr)

    allocate(dbuf3(nzonal,nmerid,nverti+1))
    dbuf3=conc(ns,1:nzonal,1:nmerid,:)
    call mpi_send(dbuf3, nzonal*nmerid*(nverti+1), mpi_double_precision, 0, iar_conc, mpi_comm_world,ierr)
    deallocate(dbuf3)
    allocate(dbuf3(nzonal,nmerid,nverti))
    dbuf3=conco(ns,1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3,nzonal*nmerid*nverti, mpi_double_precision, 0, iar_conco, mpi_comm_world,ierr)
    deallocate(dbuf3)
  end subroutine worker_send_conc

!*****************************************
  subroutine worker_send_conc_tl(ns)

    integer :: ns

    integer :: ierr
    real(kind=8),allocatable,dimension(:,:,:) :: dbuf6

    allocate(dbuf6(nzonal,nmerid,nverti+1))
    dbuf6=conc_tl(ns,1:nzonal,1:nmerid,:)
    call mpi_send(dbuf6, nzonal*nmerid*(nverti+1), mpi_double_precision, 0, iar_conc_tl, mpi_comm_world,ierr)
    deallocate(dbuf6)
    allocate(dbuf6(nzonal,nmerid,nverti))
    dbuf6=conco_tl(ns,1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf6, nzonal*nmerid*nverti, mpi_double_precision, 0, iar_conco_tl, mpi_comm_world,ierr)
    deallocate(dbuf6)

  end subroutine worker_send_conc_tl
 !*****************************************
  subroutine aworker_send_aconc

    integer :: ierr,ime,izo,ns,ivert
    real(kind=8),allocatable,dimension(:,:,:,:) :: dbuf4

    allocate(dbuf4(nspectot,nzonal+6,nmerid+6,nverti+1))
    do ns=1,nspectot
    do ivert=1,nverti+1
    do izo=1,nzonal+6
    do ime=1,nmerid+6
       dbuf4(ns,izo,ime,ivert)=aconc(ns,izo-3,ime-3,ivert)
    enddo
    enddo
    enddo
    enddo
!    dbuf4(:,1:nzonal+6,1:nmerid+6,:)=aconc(:,-2:nzonal+3,-2:nmerid+3,:)
    call mpi_send(dbuf4, nspectot*(nzonal+6)*(nmerid+6)*(nverti+1), mpi_double_precision, 0, iar_aconc, mpi_comm_world,ierr)
    deallocate(dbuf4)
    aconc(:,:,:,:)=0.d0

    allocate(dbuf4(nspectot,nzonal,nmerid,nverti))
    dbuf4=aconco(:,:,:,:)
    call mpi_send(dbuf4, nspectot*nzonal*nmerid*nverti, mpi_double_precision, 0, iar_aconco, mpi_comm_world,ierr)
    deallocate(dbuf4)
    aconco=0.d0

  end subroutine aworker_send_aconc

  !*****************************************
  subroutine worker_recv_ns(ns)
    implicit none
    integer :: ns
    integer :: ierr

    call mpi_recv(ns,1,mpi_integer,0,ias_ns,mpi_comm_world,mpi_status_ignore,ierr)
  end subroutine worker_recv_ns


  !*****************************************
  subroutine worker_recv_conc(ns)
    implicit none
    integer :: ns
    integer :: ierr
    real(kind=8),dimension(:,:,:),allocatable :: dbuf3
    allocate(dbuf3(-2:nzonal+3,-2:nmerid+3,nverti+1))
    call mpi_recv(dbuf3,                            &
         (nzonal+6)*(nmerid+6)*(nverti+1), &
         mpi_double_precision,                      &
         0,                                         &
         ias_conc,                                  &
         mpi_comm_world,mpi_status_ignore,ierr)
    conc(ns,-2:nzonal+3,-2:nmerid+3,:) = dbuf3
    deallocate(dbuf3)
    allocate(dbuf3(nzonal,nmerid,nverti))
    call mpi_recv(dbuf3,                            &
         nzonal*nmerid*nverti, &
         mpi_double_precision,                      &
         0,                                         &
         ias_conco,                                 &
         mpi_comm_world,mpi_status_ignore,ierr)
    conco(ns,1:nzonal,1:nmerid,1:nverti) = dbuf3
    deallocate(dbuf3)

  end subroutine worker_recv_conc
!*****************************************
  subroutine worker_recv_conc_tl(ns)
    implicit none
    integer :: ns
    integer :: ierr
    real(kind=8),dimension(:,:,:),allocatable :: dbuf6

    allocate(dbuf6(-2:nzonal+3,-2:nmerid+3,nverti+1 ))
    call mpi_recv(dbuf6,(nzonal+6)*(nmerid+6)*(nverti+1) &
                     ,mpi_double_precision,0,ias_conc_tl, &
                      mpi_comm_world,mpi_status_ignore,ierr)
    conc_tl(ns,-2:nzonal+3,-2:nmerid+3,:)=dbuf6
    deallocate(dbuf6)
    allocate(dbuf6(nzonal,nmerid,nverti))
    call mpi_recv(dbuf6,                            &
         nzonal*nmerid*nverti, &
         mpi_double_precision,                      &
         0,                                         &
         ias_conco_tl,                              &
         mpi_comm_world,mpi_status_ignore,ierr)
    conco_tl(ns,1:nzonal,1:nmerid,1:nverti) = dbuf6
    deallocate(dbuf6)

    !print*,'test send',rank,nmeridmax+3,nmerid,conc_tl(1,1,nmerid:nmeridmax+3)

  end subroutine worker_recv_conc_tl
   !*****************************************
  subroutine aworker_recv_aconc
    implicit none
    integer :: ierr
    real(kind=8),dimension(:,:,:,:),allocatable :: dbuf4

    allocate(dbuf4(nspectot,-2:nzonal+3,-2:nmerid+3,nverti+1))
    call mpi_recv(dbuf4,                            &
         nspectot*(nzonal+6)*(nmerid+6)*(nverti+1), &
         mpi_double_precision,                      &
         0,                                         &
         ias_aconc,                                  &
         mpi_comm_world,mpi_status_ignore,ierr)
    aconc(:,-2:nzonal+3,-2:nmerid+3,:) = &
       aconc(:,-2:nzonal+3,-2:nmerid+3,:) + dbuf4
    deallocate(dbuf4)

    allocate(dbuf4(nspectot,nzonal,nmerid,nverti))
    call mpi_recv(dbuf4,                            &
         nspectot*(nzonal)*(nmerid)*nverti, &
         mpi_double_precision,                      &
         0,                                         &
         ias_aconco,                                  &
         mpi_comm_world,mpi_status_ignore,ierr)
    aconco(:,1:nzonal,1:nmerid,:) = &
        aconco(:,1:nzonal,1:nmerid,:) + dbuf4
    deallocate(dbuf4)

  end subroutine aworker_recv_aconc
   !*****************************************
  subroutine worker_recv_conc_bounds(ns)
    implicit none
    integer :: ns

    integer :: ierr
    real(kind=8),allocatable,dimension(:,:,:) :: hbuf,vbuf
    real(kind=8),allocatable,dimension(:,:)   :: topbuf

    allocate(vbuf(3,         -2:nmerid+3,nverti+1))
    allocate(hbuf(-2:nzonal+3,3,         nverti+1))
    allocate(topbuf(nzonal,nmerid))

    !print*,rank,'wwwwwwwwwwwwwwww',dom_i,nzdoms
    if(dom_i==1) then
       !print*,rank,ns,' recv left boundary',nmerid+6,nverti+1
       call mpi_recv(vbuf,                             &
            3*(nmerid+6)*(nverti+1),          &
            mpi_double_precision,                      &
            0,                                         &
            ias_conc_bounds,                           &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc(ns,-2:0,-2:nmerid+3,:) = vbuf(:,:,:)
    end if

    if(dom_i==nzdoms) then
       !print*,rank,ns,' recv right boundary'
       call mpi_recv(vbuf,                             &
            3*(nmerid+6)*(nverti+1),          &
            mpi_double_precision,                      &
            0,                                         &
            ias_conc_bounds,                           &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc(ns,nzonal+1:nzonal+3,-2:nmerid+3,:) = vbuf(:,:,:)
    end if

    if(dom_j==1) then
       !print*,rank,ns,' recv lower boundary'
       call mpi_recv(hbuf,                          &
            (nzonal+6)*3*(nverti+1),       &
            mpi_double_precision,                   &
            0,                                      &
            ias_conc_bounds,                        &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc(ns,-2:nzonal+3,-2:0,:) = hbuf(:,1:3,:)
    end if

    if(dom_j==nmdoms) then
       !print*,rank,ns,' recv upper boundary'
       call mpi_recv(hbuf,                          &
            (nzonal+6)*3*(nverti+1),       &
            mpi_double_precision,                   &
            0,                                      &
            ias_conc_bounds,                        &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc(ns,-2:nzonal+3,nmerid+1:nmerid+3,:) = hbuf(:,:,:)
    end if

    !print*,rank,ns,' recv top boundary'
    call mpi_recv(topbuf,                        &
         nzonal*nmerid,                 &
         mpi_double_precision,                   &
         0,                                      &
         ias_conc_top,                           &
         mpi_comm_world,mpi_status_ignore,ierr)
    conc(ns,1:nzonal,1:nmerid,nverti+1) = topbuf(:,:)
    deallocate(hbuf)
    deallocate(vbuf)
    deallocate(topbuf)

  end subroutine worker_recv_conc_bounds

!*****************************************
  subroutine worker_recv_conc_bounds_tl(ns)
    implicit none

    integer :: ns

    integer :: ierr

    real(kind=8),allocatable,dimension(:,:,:) :: hbuf,vbuf
    real(kind=8),allocatable,dimension(:,:)   :: topbuf

    allocate(vbuf(3,         -2:nmerid+3,nverti+1))
    allocate(hbuf(-2:nzonal+3,3,         nverti+1))
    allocate(topbuf(nzonal,nmerid))

    if(dom_i==1) then
       ! recv left boundary
       call mpi_recv(vbuf,                             &
            3*(nmerid+6)*(nverti+1),          &
            mpi_double_precision,                      &
            0,                                         &
            ias_conc_bounds_tl,                           &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc_tl(ns,-2:0,-2:nmerid+3,:) = vbuf(:,:,:)
    end if

    if(dom_i==nzdoms) then
       ! recv right boundary
       call mpi_recv(vbuf,                             &
            3*(nmerid+6)*(nverti+1),          &
            mpi_double_precision,                      &
            0,                                         &
            ias_conc_bounds_tl,                           &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc_tl(ns,nzonal+1:nzonal+3,-2:nmerid+3,:) = vbuf(:,:,:)
    end if

    if(dom_j==1) then
       ! recv lower boundary
       call mpi_recv(hbuf,                          &
            (nzonal+6)*3*(nverti+1),       &
            mpi_double_precision,                   &
            0,                                      &
            ias_conc_bounds_tl,                        &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc_tl(ns,-2:nzonal+3,-2:0,:) = hbuf(:,1:3,:)
    end if

    if(dom_j==nmdoms) then
       ! recv upper boundary
       call mpi_recv(hbuf,                          &
            (nzonal+6)*3*(nverti+1),       &
            mpi_double_precision,                   &
            0,                                      &
            ias_conc_bounds_tl,                        &
            mpi_comm_world,mpi_status_ignore,ierr)
       conc_tl(ns,-2:nzonal+3,nmerid+1:nmerid+3,:) = hbuf(:,:,:)
    end if

        ! recv top boundary
    call mpi_recv(topbuf,                        &
         nzonal*nmerid,                 &
         mpi_double_precision,                   &
         0,                                      &
         ias_conc_top_tl,                           &
         mpi_comm_world,mpi_status_ignore,ierr)
    conc_tl(ns,1:nzonal,1:nmerid,nverti+1) = topbuf(:,:)

    deallocate(hbuf)
    deallocate(vbuf)
    deallocate(topbuf)

  end subroutine worker_recv_conc_bounds_tl
  !*****************************************
  subroutine aworker_send_aconc_bounds
    integer :: ierr
    real(kind=8),allocatable,dimension(:,:,:) :: topbuf
    real(kind=8),allocatable,dimension(:,:,:,:) :: hbuf,vbuf

    allocate(topbuf(nspectot,nzonal,nmerid))
    topbuf=aconc(:,1:nzonal,1:nmerid,nverti+1)
    call mpi_send(topbuf,nspectot*nzonal*nmerid, mpi_double_precision, 0, iar_aconc_top, mpi_comm_world,ierr)
    deallocate(topbuf)
! if(dom_i==1)write(*,1001)aconc(1,9,1,nverti+1)
!1001 format('WW ',e64.56)
    aconc(:,1:nzonal,1:nmerid,nverti+1)=0.d0

    if(dom_i==1) then ! left boundary
      allocate(vbuf(nspectot,3,nmerid+6,nverti+1))
      vbuf(:,:,:,:)=aconc(:,-2:0,-2:nmerid+3,:)
      call mpi_send(vbuf,                             &
            nspectot*3*(nmerid+6)*(nverti+1),          &
            mpi_double_precision,                      &
            0,                                         &
            iar_aconc_bounds,                           &
            mpi_comm_world,mpi_status_ignore,ierr)
      deallocate(vbuf)
      aconc(:,-2:0,-2:nmerid+3,:)=0.d0
    endif

    if(dom_i==nzdoms) then ! right boundary
       allocate(vbuf(nspectot,3,nmerid+6,nverti+1))
       vbuf(:,:,:,:)=aconc(:,nzonal+1:nzonal+3,-2:nmerid+3,:)
       call mpi_send(vbuf,                             &
            nspectot*3*(nmerid+6)*(nverti+1),          &
            mpi_double_precision,                      &
            0,                                         &
            iar_aconc_bounds,                           &
            mpi_comm_world,mpi_status_ignore,ierr)
       deallocate(vbuf)
       aconc(:,nzonal+1:nzonal+3,-2:nmerid+3,:)=0.d0
    end if

    if(dom_j==1) then ! lower boundary
       allocate(hbuf(nspectot,-2:nzonal+3,3,nverti+1))
       hbuf(:,:,:,:)=aconc(:, -2:nzonal+3, -2:0, :)
       call mpi_send(hbuf,                          &
            nspectot*(nzonal+6)*3*(nverti+1),       &
            mpi_double_precision,                   &
            0,                                      &
            iar_aconc_bounds,                        &
            mpi_comm_world,mpi_status_ignore,ierr)
       deallocate(hbuf)
       aconc(:, -2:nzonal+3, -2:0, :)=0.d0
    end if

    if(dom_j==nmdoms) then ! upper boundary
       allocate(hbuf(nspectot,-2:nzonal+3,3,nverti+1))
       hbuf(:,:,:,:)=aconc(:,-2:nzonal+3,nmerid+1:nmerid+3,:)
       call mpi_send(hbuf,                          &
            nspectot*(nzonal+6)*3*(nverti+1),       &
            mpi_double_precision,                   &
            0,                                      &
            iar_aconc_bounds,                        &
            mpi_comm_world,mpi_status_ignore,ierr)
       deallocate(hbuf)
       aconc(:,-2:nzonal+3,nmerid+1:nmerid+3,:)=0.d0
    end if

  end subroutine aworker_send_aconc_bounds
  !*****************************************
  subroutine worker_update_halo
    ! update upper and lower halos

    integer :: ierr
    integer :: send_upper_req,send_lower_req,recv_upper_req,recv_lower_req
    real(kind=8),allocatable,dimension(:,:,:,:) :: send_upper_buf
    real(kind=8),allocatable,dimension(:,:,:,:) :: send_lower_buf
    real(kind=8),allocatable,dimension(:,:,:,:) :: recv_upper_buf
    real(kind=8),allocatable,dimension(:,:,:,:) :: recv_lower_buf


    allocate(send_upper_buf(nspectot,nzonal,3,nverti))
    allocate(send_lower_buf(nspectot,nzonal,3,nverti))
    allocate(recv_upper_buf(nspectot,nzonal,3,nverti))
    allocate(recv_lower_buf(nspectot,nzonal,3,nverti))
    if((dom_j==1).and.(nmdoms>1)) then
       call mpi_irecv(                      &
            recv_upper_buf,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            500+rank+nzdoms,                &
            wrk_comm,                     &
            recv_upper_req,                 &
            ierr                            &
            )
       send_upper_buf=conc(:,1:nzonal,nmerid-2:nmerid,1:nverti)
       call mpi_issend(                     &
            send_upper_buf,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            400+rank,                       &
            wrk_comm,                     &
            send_upper_req,                 &
            ierr                            &
            )
       call mpi_wait(recv_upper_req,mpi_status_ignore,ierr)
       conc(:,1:nzonal,nmerid+1:nmerid+3,1:nverti) = recv_upper_buf
       call mpi_wait(send_upper_req,mpi_status_ignore,ierr)
    end if


    if((dom_j>1).and.(dom_j<nmdoms)) then

       call mpi_irecv(                      &
            recv_upper_buf,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            500+rank+nzdoms,                &
            wrk_comm,                     &
            recv_upper_req,                 &
            ierr                            &
            )
       send_upper_buf=conc(:,1:nzonal,nmerid-2:nmerid,1:nverti)
       call mpi_issend(                     &
            send_upper_buf,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            400+rank,                       &
            wrk_comm,                     &
            send_upper_req,                 &
            ierr                            &
            )

       call mpi_irecv(                                   &
            recv_lower_buf,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            400+rank-nzdoms,                             &
            wrk_comm,                                  &
            recv_lower_req,                              &
            ierr                                         &
            )
       send_lower_buf=conc(:,1:nzonal,1:3,1:nverti)
       call mpi_issend(                                  &
            send_lower_buf,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            500+rank,                                    &
            wrk_comm,                                  &
            send_lower_req,                              &
            ierr                                         &
            )
       call mpi_wait(recv_upper_req,mpi_status_ignore,ierr)
       conc(:,1:nzonal,nmerid+1:nmerid+3,1:nverti) = recv_upper_buf
       call mpi_wait(recv_lower_req,mpi_status_ignore,ierr)
       conc(:,1:nzonal,-2:0,1:nverti) = recv_lower_buf
       call mpi_wait(send_upper_req,mpi_status_ignore,ierr)
       call mpi_wait(send_lower_req,mpi_status_ignore,ierr)
    end if

    if((dom_j==nmdoms).and.(nmdoms>1)) then
       call mpi_irecv(                                   &
            recv_lower_buf,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            400+rank-nzdoms,                             &
            wrk_comm,                                  &
            recv_lower_req,                              &
            ierr                                         &
            )
       send_lower_buf=conc(:,1:nzonal,1:3,1:nverti)
       call mpi_issend(                                  &
            send_lower_buf,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            500+rank,                                    &
            wrk_comm,                                  &
            send_lower_req,                              &
            ierr                                         &
            )
       call mpi_wait(recv_lower_req,mpi_status_ignore,ierr)
       conc(:,1:nzonal,-2:0,1:nverti) = recv_lower_buf
       call mpi_wait(send_lower_req,mpi_status_ignore,ierr)
    end if
    deallocate(send_upper_buf)
    deallocate(send_lower_buf)
    deallocate(recv_upper_buf)
    deallocate(recv_lower_buf)

  end subroutine worker_update_halo

  !*****************************************
  subroutine worker_update_halo_tl
    ! update upper and lower halos

    integer :: ierr
    integer :: send_upper_req_tl,send_lower_req_tl,recv_upper_req_tl,recv_lower_req_tl
    real(kind=8),allocatable,dimension(:,:,:,:) :: send_upper_buf_tl
    real(kind=8),allocatable,dimension(:,:,:,:) :: send_lower_buf_tl
    real(kind=8),allocatable,dimension(:,:,:,:) :: recv_upper_buf_tl
    real(kind=8),allocatable,dimension(:,:,:,:) :: recv_lower_buf_tl


    allocate(send_upper_buf_tl(nspectot,nzonal,3,nverti))
    allocate(send_lower_buf_tl(nspectot,nzonal,3,nverti))
    allocate(recv_upper_buf_tl(nspectot,nzonal,3,nverti))
    allocate(recv_lower_buf_tl(nspectot,nzonal,3,nverti))
    if((dom_j==1).and.(nmdoms>1)) then
       call mpi_irecv(                      &
            recv_upper_buf_tl,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            550+rank+nzdoms,                &
            wrk_comm,                     &
            recv_upper_req_tl,                 &
            ierr                            &
            )
       send_upper_buf_tl=conc_tl(:,1:nzonal,nmerid-2:nmerid,1:nverti)
       call mpi_issend(                     &
            send_upper_buf_tl,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            450+rank,                       &
            wrk_comm,                     &
            send_upper_req_tl,                 &
            ierr                            &
            )
       call mpi_wait(recv_upper_req_tl,mpi_status_ignore,ierr)
       conc_tl(:,1:nzonal,nmerid+1:nmerid+3,1:nverti) = recv_upper_buf_tl
       call mpi_wait(send_upper_req_tl,mpi_status_ignore,ierr)
    end if


    if((dom_j>1).and.(dom_j<nmdoms)) then

       call mpi_irecv(                      &
            recv_upper_buf_tl,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            550+rank+nzdoms,                &
            wrk_comm,                     &
            recv_upper_req_tl,                 &
            ierr                            &
            )
       send_upper_buf_tl=conc_tl(:,1:nzonal,nmerid-2:nmerid,1:nverti)
       call mpi_issend(                     &
            send_upper_buf_tl,                 &
            nspectot*nzonal*3*nverti,       &
            mpi_double_precision,           &
            rank+nzdoms,                    &
            450+rank,                       &
            wrk_comm,                     &
            send_upper_req_tl,                 &
            ierr                            &
            )

       call mpi_irecv(                                   &
            recv_lower_buf_tl,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            450+rank-nzdoms,                             &
            wrk_comm,                                  &
            recv_lower_req_tl,                              &
            ierr                                         &
            )
       send_lower_buf_tl=conc_tl(:,1:nzonal,1:3,1:nverti)
       call mpi_issend(                                  &
            send_lower_buf_tl,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            550+rank,                                    &
            wrk_comm,                                  &
            send_lower_req_tl,                              &
            ierr                                         &
            )
       call mpi_wait(recv_upper_req_tl,mpi_status_ignore,ierr)
       conc_tl(:,1:nzonal,nmerid+1:nmerid+3,1:nverti) = recv_upper_buf_tl
       call mpi_wait(recv_lower_req_tl,mpi_status_ignore,ierr)
       conc_tl(:,1:nzonal,-2:0,1:nverti) = recv_lower_buf_tl
       call mpi_wait(send_upper_req_tl,mpi_status_ignore,ierr)
       call mpi_wait(send_lower_req_tl,mpi_status_ignore,ierr)
    end if

    if((dom_j==nmdoms).and.(nmdoms>1)) then
       call mpi_irecv(                                   &
            recv_lower_buf_tl,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            450+rank-nzdoms,                             &
            wrk_comm,                                  &
            recv_lower_req_tl,                              &
            ierr                                         &
            )
       send_lower_buf_tl=conc_tl(:,1:nzonal,1:3,1:nverti)
       call mpi_issend(                                  &
            send_lower_buf_tl,                              &
            nspectot*nzonal*3*nverti,                    &
            mpi_double_precision,                        &
            rank-nzdoms,                                 &
            550+rank,                                    &
            wrk_comm,                                  &
            send_lower_req_tl,                              &
            ierr                                         &
            )
       call mpi_wait(recv_lower_req_tl,mpi_status_ignore,ierr)
       conc_tl(:,1:nzonal,-2:0,1:nverti) = recv_lower_buf_tl
       call mpi_wait(send_lower_req_tl,mpi_status_ignore,ierr)
    end if
    deallocate(send_upper_buf_tl)
    deallocate(send_lower_buf_tl)
    deallocate(recv_upper_buf_tl)
    deallocate(recv_lower_buf_tl)

  end subroutine worker_update_halo_tl
  !*****************************************
  subroutine rcopy_int_params

    ! This routine explodes the array int_params into
    ! various well-known chimere scalars

    ! It use a macro to shorten code writing and reduce risk of typo errors
    ! The macro generates lines like :
    ! nverti = int_params(ip_nverti)
    !
    ! ip_nverti beeing an arbitrary but unique index defined in module message_defs
#define R_INT_PARAMS(par)             par = int_params(ip_##par)
    implicit none

    R_INT_PARAMS(nzonal)
    R_INT_PARAMS(nmerid)
    R_INT_PARAMS(dom_i)
    R_INT_PARAMS(dom_j)
    R_INT_PARAMS(wrk_comm)
    R_INT_PARAMS(imstart)
    R_INT_PARAMS(izstart)
    R_INT_PARAMS(imend)
    R_INT_PARAMS(izend)
    R_INT_PARAMS(nspec)
    R_INT_PARAMS(nspectot)
    R_INT_PARAMS(nreac)
    R_INT_PARAMS(nfam)
    R_INT_PARAMS(ndepo)
    R_INT_PARAMS(nreactamax)
    R_INT_PARAMS(ntabuzen)
    R_INT_PARAMS(nlevphot)
    R_INT_PARAMS(nphot)
    R_INT_PARAMS(nhourrun)
    R_INT_PARAMS(ntemps)
    R_INT_PARAMS(ntabmax)
    R_INT_PARAMS(nverti)
    R_INT_PARAMS(nphour_ref)
    R_INT_PARAMS(ichemstep)
    R_INT_PARAMS(usechemistry)
    R_INT_PARAMS(usedepos)
    R_INT_PARAMS(usewetdepos)
    R_INT_PARAMS(useemissions)
    R_INT_PARAMS(usetransmix)
    R_INT_PARAMS(useabsclipconc)
    R_INT_PARAMS(dryairout)
    R_INT_PARAMS(nvegtype)
    R_INT_PARAMS(nlduse)
    R_INT_PARAMS(nsaveconcs)
    R_INT_PARAMS(noutspec)
    R_INT_PARAMS(nzdoms)
    R_INT_PARAMS(nmdoms)
    R_INT_PARAMS(nzonalmax)
    R_INT_PARAMS(nmeridmax)
    R_INT_PARAMS(ideepconv)
    R_INT_PARAMS(nitgs)
    R_INT_PARAMS(nitgssu)
    R_INT_PARAMS(optemisb)
    R_INT_PARAMS(nemisa)
    R_INT_PARAMS(nemisb)
    R_INT_PARAMS(ndep)
    R_INT_PARAMS(nlevemis)
    R_INT_PARAMS(ihoursu)
    R_INT_PARAMS(fwd)
    R_INT_PARAMS(tl)
    R_INT_PARAMS(ad)
    R_INT_PARAMS(hpulse)
  end subroutine rcopy_int_params

  !*****************************************
  subroutine rcopy_dbl_params
    ! This routine explodes the array dbl_params into
    ! various well-known chimere scalars

    ! It use a macro to shorten code writing and reduce risk of typo errors
    ! The macro generates lines like :
    ! factorl = dbl_params(ip_factorl)
    !
    ! ip_factorl beeing an arbitrary but unique index defined in module message_defs

#define R_DBL_PARAMS(par)             par = dbl_params(ip_##par)
    implicit none

    R_DBL_PARAMS(soltim)
    R_DBL_PARAMS(djul)
    R_DBL_PARAMS(clipconc)
    R_DBL_PARAMS(psurf)
  end subroutine rcopy_dbl_params

  !*****************************************
  subroutine recv_int_arrays
    implicit none
    integer :: ierr
    integer,dimension(:),    allocatable :: species_varid,species_transp,species_transpv,species_bounddry
    integer,dimension(:,:,:),allocatable :: ibuf3

    call mpi_recv(inemisa,   nspec,           mpi_integer,0,ias_inemisa, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(inemisb,   nspec,           mpi_integer,0,ias_inemisb,   mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(kreacl,    nspectot,        mpi_integer,0,ias_kreacl,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(kreacp,    nspectot,        mpi_integer,0,ias_kreacp,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(ireacl,    nspectot*nreac,  mpi_integer,0,ias_ireacl,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(ireacp,    nspectot*nreac,  mpi_integer,0,ias_ireacp,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(nreactants,nreac,           mpi_integer,0,ias_nreactants,mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(irctt,     nreac*nreactamax,mpi_integer,0,ias_irctt,     mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(indepo,    nspec+1,         mpi_integer,0,ias_indepo,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(inwetd,    nspec*2,         mpi_integer,0,ias_inwetd,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(imonth,    nhourrun+1,      mpi_integer,0,ias_imonth,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(ityperate, nreac,           mpi_integer,0,ias_ityperate, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(nwetd,     2,               mpi_integer,0,ias_nwetd,     mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(iphoto,    nreac,           mpi_integer,0,ias_iphoto,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(ispecemip, 1000,            mpi_integer,0,ias_ispecemip, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(nelem      ,nfam,             mpi_integer,0,ias_nelem,         mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(ifam       ,nfam*nspectot*10, mpi_integer,0,ias_ifam,          mpi_comm_world,mpi_status_ignore,ierr)


    allocate(species_varid(nspectot))
    call mpi_recv(species_varid,nspectot,     mpi_integer,0,ias_species_varid, mpi_comm_world,mpi_status_ignore,ierr)
    species(:)%varid=species_varid(:)
    deallocate(species_varid)

    allocate(species_transp(nspectot))
    call mpi_recv(species_transp,nspectot,     mpi_integer,0,ias_species_transp, mpi_comm_world,mpi_status_ignore,ierr)
    species(:)%transp=species_transp(:)
    deallocate(species_transp)

    allocate(species_transpv(nspectot))
    call mpi_recv(species_transpv,nspectot,     mpi_integer,0,ias_species_transpv, mpi_comm_world,mpi_status_ignore,ierr)
    species(:)%transpv=species_transpv(:)
    deallocate(species_transpv)

    allocate(species_bounddry(nspectot))
    call mpi_recv(species_bounddry,nspectot,     mpi_integer,0,ias_species_bounddry, mpi_comm_world,mpi_status_ignore,ierr)
    species(:)%bounddry=species_bounddry(:)
    deallocate(species_bounddry)

    call mpi_recv(ihour,nhourrun+1,           mpi_integer,0,ias_ihour,   mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(idtyp,nhourrun+1,           mpi_integer,0,ias_idtyp,   mpi_comm_world,mpi_status_ignore,ierr)

  end subroutine recv_int_arrays


  !*****************************************
  subroutine recv_char_arrays
    implicit none
    integer :: ierr
    character(len=16),dimension(:),allocatable :: species_name

    allocate(species_name(nspectot))
    call mpi_recv(species_name, 16*nspectot,  mpi_character,0,ias_species_name, mpi_comm_world,mpi_status_ignore,ierr)
    species(:)%name=species_name(:)
    deallocate(species_name)

  end subroutine recv_char_arrays


  !*****************************************
  subroutine recv_real_arrays
    implicit none
    integer :: ierr

    call mpi_recv(altiphot,                &
         nlevphot,mpi_double_precision, &
         0,                                &
         ias_altiphot,                     &
         mpi_comm_world,                       &
         mpi_status_ignore,ierr)
    call mpi_recv(photoj,                  &
         ntabuzen*nlevphot*nphot, &
         mpi_double_precision,             &
         0,                                &
         ias_photoj,                       &
         mpi_comm_world,                       &
         mpi_status_ignore,ierr)
  end subroutine recv_real_arrays

  !*****************************************
   subroutine recv_output_species
    implicit none
    integer :: ierr
    real(kind=8),        dimension(nspectot) :: fspec
    integer,             dimension(nspectot) :: varid, iaddr
    character(len=splen),dimension(nspectot) :: name
    character(len=16),   dimension(nspectot) :: units

    call mpi_recv(fspec,nspectot,      mpi_double_precision, 0,ias_os_fscpe,mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(varid,nspectot,      mpi_integer,          0,ias_os_varid,mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(iaddr,nspectot,      mpi_integer,          0,ias_os_iaddr,mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(name, splen*nspectot,mpi_character,        0,ias_os_name, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(units,   16*nspectot,mpi_character,        0,ias_os_units,mpi_comm_world,mpi_status_ignore,ierr)

    species(:)%fspec = fspec(:)
    output_species(:)%varid = varid(:)
    output_species(:)%iaddr = iaddr(:)
    output_species(:)%name  = name(:)
    output_species(:)%units = units(:)

  end subroutine recv_output_species

    !*****************************************
subroutine recv_dbl_arrays
    implicit none
    integer :: ierr
    real(kind=8),dimension(:,:),    allocatable   :: dbuf2
    real(kind=8),dimension(:,:,:),  allocatable   :: dbuf3
    real(kind=8),dimension(:,:,:,:),allocatable   :: dbuf4
    real(kind=8),dimension(:,:,:,:,:),allocatable :: dbuf5

    call mpi_recv(stoi,    nspectot*nreac*ntemps,        mpi_double_precision,0,ias_stoi,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(deptmin, nvegtype,                     mpi_double_precision,0,ias_deptmin, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(deptopt, nvegtype,                     mpi_double_precision,0,ias_deptopt, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depvpd1, nvegtype,                     mpi_double_precision,0,ias_depvpd1, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depvpd2, nvegtype,                     mpi_double_precision,0,ias_depvpd2, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depegs,  nvegtype,                     mpi_double_precision,0,ias_depegs,  mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depsgs,  nvegtype,                     mpi_double_precision,0,ias_depsgs,  mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depegl,  nvegtype,                     mpi_double_precision,0,ias_depegl,  mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depsgl,  nvegtype,                     mpi_double_precision,0,ias_depsgl,  mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(deplai1, nvegtype,                     mpi_double_precision,0,ias_deplai1, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(deplai2, nvegtype,                     mpi_double_precision,0,ias_deplai2, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depphe0, nvegtype,                     mpi_double_precision,0,ias_depphe0, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depphe1, nvegtype,                     mpi_double_precision,0,ias_depphe1, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depphe2, nvegtype,                     mpi_double_precision,0,ias_depphe2, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(depalph, nvegtype,                     mpi_double_precision,0,ias_depalph, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(gmax,    nvegtype,                     mpi_double_precision,0,ias_gmax,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(fmin,    nvegtype,                     mpi_double_precision,0,ias_fmin,    mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(zcanopy, nvegtype,                     mpi_double_precision,0,ias_zcanopy, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(rgso3,   nvegtype,                     mpi_double_precision,0,ias_rgso3,   mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(so2rh,   nvegtype,                     mpi_double_precision,0,ias_so2rh,   mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(rgsso2,  nvegtype,                     mpi_double_precision,0,ias_rgsso2,  mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(factrb,  nspec,                        mpi_double_precision,0,ias_factrb,  mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(factd,   nspec,                        mpi_double_precision,0,ias_factd,   mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(rm,      nspec,                        mpi_double_precision,0,ias_rm,      mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(dhx,     nspec,                        mpi_double_precision,0,ias_dhx,     mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(df0,     nspec,                        mpi_double_precision,0,ias_df0,     mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(zetaref, ntabuzen,                  mpi_double_precision,0,ias_zetaref, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(tabrate, ntabmax*nreac,                mpi_double_precision,0,ias_tabrate, mpi_comm_world,mpi_status_ignore,ierr)
    call mpi_recv(tabtemp, ntemps,                       mpi_double_precision,0,ias_tabtemp, mpi_comm_world,mpi_status_ignore,ierr)
    allocate(dbuf2(nzonal+4,nmerid+4))
    call mpi_recv(dbuf2,   (nzonal+4)*(nmerid+4),        mpi_double_precision,0,ias_xsize,   mpi_comm_world,mpi_status_ignore,ierr)
    xsize(-1:nzonal+2,-1:nmerid+2) = dbuf2(:,:)
    call mpi_recv(dbuf2,   (nzonal+4)*(nmerid+4),        mpi_double_precision,0,ias_ysize,   mpi_comm_world,mpi_status_ignore,ierr)
    ysize(-1:nzonal+2,-1:nmerid+2) = dbuf2(:,:)
    deallocate(dbuf2)
    allocate(dbuf2(nzonal+2,nmerid+2))
    call mpi_recv(dbuf2,   (nzonal+2)*(nmerid+2),        mpi_double_precision,0,ias_xbasx,   mpi_comm_world,mpi_status_ignore,ierr)
    xbasx(0:nzonal+1,0:nmerid+1) = dbuf2(:,:)
    call mpi_recv(dbuf2,   (nzonal+2)*(nmerid+2),        mpi_double_precision,0,ias_xbasy,   mpi_comm_world,mpi_status_ignore,ierr)
    xbasy(0:nzonal+1,0:nmerid+1) = dbuf2(:,:)
    call mpi_recv(dbuf2,   (nzonal+2)*(nmerid+2),        mpi_double_precision,0,ias_ybasx,   mpi_comm_world,mpi_status_ignore,ierr)
    ybasx(0:nzonal+1,0:nmerid+1) = dbuf2(:,:)
    call mpi_recv(dbuf2,   (nzonal+2)*(nmerid+2),        mpi_double_precision,0,ias_ybasy,   mpi_comm_world,mpi_status_ignore,ierr)
    ybasy(0:nzonal+1,0:nmerid+1) = dbuf2(:,:)
    deallocate(dbuf2)

    allocate(dbuf2(nzonal,nmerid))
    call mpi_recv(dbuf2,   (nzonal)*(nmerid),            mpi_double_precision,0,ias_slati,   mpi_comm_world,mpi_status_ignore,ierr)
    slati(1:nzonal,1:nmerid) = dbuf2(:,:)
    call mpi_recv(dbuf2,   (nzonal)*(nmerid),            mpi_double_precision,0,ias_clati,   mpi_comm_world,mpi_status_ignore,ierr)
    clati(1:nzonal,1:nmerid) = dbuf2(:,:)
    call mpi_recv(dbuf2,   (nzonal)*(nmerid),            mpi_double_precision,0,ias_xlong,   mpi_comm_world,mpi_status_ignore,ierr)
    xlong(1:nzonal,1:nmerid) = dbuf2(:,:)
    call mpi_recv(dbuf2,   (nzonal)*(nmerid),            mpi_double_precision,0,ias_xlati,   mpi_comm_world,mpi_status_ignore,ierr)
    xlati(1:nzonal,1:nmerid) = dbuf2(:,:)
    deallocate(dbuf2)

    allocate(dbuf4(nzonal,nmerid,nvegtype,nlduse))
    call mpi_recv(dbuf4,    nzonal*nmerid*nvegtype*nlduse,mpi_double_precision,0,ias_fveg,    mpi_comm_world,mpi_status_ignore,ierr)
    fveg(1:nzonal,1:nmerid,:,:) = dbuf4(:,:,:,:)
    deallocate(dbuf4)

    allocate(dbuf3(nzonal,nmerid,nlduse))
    call mpi_recv(dbuf3,   nzonal*nmerid*nlduse,         mpi_double_precision,0,ias_dland,   mpi_comm_world,mpi_status_ignore,ierr)
    dland(1:nzonal,1:nmerid,1:nlduse) = dbuf3(:,:,:)
    deallocate(dbuf3)

  end subroutine recv_dbl_arrays


  !*****************************************
  subroutine recv_hourly_real_arrays
    implicit none
    integer :: ierr
    real(kind=8),dimension(:,:,:),  allocatable :: rbuf3
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf4

    allocate(rbuf4(nemisa,nzonal,nmerid,nlevemis))
    call mpi_recv(rbuf4,nemisa*nzonal*nmerid*nlevemis, mpi_double_precision,0, ias_emisaloc, mpi_comm_world,mpi_status_ignore,ierr)
    emisaloc(:,1:nzonal,1:nmerid,:) = rbuf4(:,:,:,:)
    deallocate(rbuf4)

    allocate(rbuf4(nemisb,nzonal,nmerid,2))
    call mpi_recv(rbuf4,nemisb*nzonal*nmerid*2,mpi_double_precision,0, ias_emisb,mpi_comm_world,mpi_status_ignore,ierr)
    emisb(:,1:nzonal,1:nmerid,:) = rbuf4(:,:,:,:)
    deallocate(rbuf4)

    allocate(rbuf4(nzonal,nmerid,nverti,2))
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2,mpi_double_precision,0, ias_temp, mpi_comm_world,mpi_status_ignore,ierr)
    temp(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision,0, ias_sphu, mpi_comm_world,mpi_status_ignore,ierr)
    sphu(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2,mpi_double_precision ,0, ias_airm, mpi_comm_world,mpi_status_ignore,ierr)
    airm(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision,0, ias_kzzz, mpi_comm_world,mpi_status_ignore,ierr)
    kzzz(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2,mpi_double_precision ,0, ias_clwc, mpi_comm_world,mpi_status_ignore,ierr)
    clwc(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
! lmbb deepconv
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision,0, ias_dpeu, mpi_comm_world,mpi_status_ignore,ierr)
    dpeu(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision,0, ias_dped, mpi_comm_world,mpi_status_ignore,ierr)
    dped(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision,0, ias_dpdu, mpi_comm_world,mpi_status_ignore,ierr)
    dpdu(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision,0, ias_dpdd, mpi_comm_world,mpi_status_ignore,ierr)
    dpdd(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)

    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2,  mpi_double_precision,0, ias_winz, mpi_comm_world,mpi_status_ignore,ierr)
    winz(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision ,0, ias_winm, mpi_comm_world,mpi_status_ignore,ierr)
    winm(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    call mpi_recv(rbuf4, nzonal*nmerid*nverti*2, mpi_double_precision ,0, ias_hlay, mpi_comm_world,mpi_status_ignore,ierr)
    hlay(1:nzonal,1:nmerid,:,:) = rbuf4(:,:,:,:)
    deallocate(rbuf4)

    allocate(rbuf3(nzonal,nmerid,2))
    call mpi_recv(rbuf3, nzonal*nmerid*2,   mpi_double_precision      ,0, ias_hght, mpi_comm_world,mpi_status_ignore,ierr)
    hght(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,     mpi_double_precision    ,0, ias_atte, mpi_comm_world,mpi_status_ignore,ierr)
    atte(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,   mpi_double_precision      ,0, ias_tem2, mpi_comm_world,mpi_status_ignore,ierr)
    tem2(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,     mpi_double_precision      ,0, ias_usta, mpi_comm_world,mpi_status_ignore,ierr)
    usta(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,   mpi_double_precision     ,0, ias_aerr, mpi_comm_world,mpi_status_ignore,ierr)
    aerr(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,    mpi_double_precision    ,0, ias_obuk, mpi_comm_world,mpi_status_ignore,ierr)
    obuk(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,    mpi_double_precision    ,0, ias_wsta, mpi_comm_world,mpi_status_ignore,ierr)
    wsta(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,    mpi_double_precision    ,0, ias_topc, mpi_comm_world,mpi_status_ignore,ierr)
    topc(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    call mpi_recv(rbuf3, nzonal*nmerid*2,     mpi_double_precision   ,0, ias_sreh, mpi_comm_world,mpi_status_ignore,ierr)
    sreh(1:nzonal,1:nmerid,:) = rbuf3(:,:,:)
    deallocate(rbuf3)

  end subroutine recv_hourly_real_arrays

! ******************************
  subroutine recv_hourly_real_arrays_tl
    implicit none
    integer :: ierr
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf5

    allocate(rbuf5(nemisa,nzonal,nmerid,nlevemis))
    call mpi_recv(rbuf5,nemisa*nzonal*nmerid*nlevemis, mpi_double_precision,0, ias_emisaloc_tl, mpi_comm_world,mpi_status_ignore,ierr)
    emisaloc_tl(:,1:nzonal,1:nmerid,:) = rbuf5(:,:,:,:)
    deallocate(rbuf5)

    allocate(rbuf5(nemisb,nzonal,nmerid,2))
    call mpi_recv(rbuf5,nemisb*nzonal*nmerid*2,mpi_double_precision,0, ias_emisb_tl,mpi_comm_world,mpi_status_ignore,ierr)
    emisb_tl(:,1:nzonal,1:nmerid,:) = rbuf5(:,:,:,:)
    deallocate(rbuf5)

  end subroutine recv_hourly_real_arrays_tl
  !*****************************************
  subroutine recv_frac_hourly_int_arrays
    implicit none
    integer :: ierr
    integer,dimension(:,:),allocatable :: ibuf2
    integer,dimension(:,:,:),allocatable :: ibuf3

    allocate(ibuf2(nzonal,nmerid))
    call mpi_recv(ibuf2, nzonal*nmerid, mpi_integer,0, ias_ideep, mpi_comm_world,mpi_status_ignore,ierr)
    ideep(1:nzonal,1:nmerid) = ibuf2(:,:)
    deallocate(ibuf2)

    allocate(ibuf3(nzonal,nmerid,nverti))
    call mpi_recv(ibuf3, nzonal*nmerid*nverti, mpi_integer,0, ias_incloud, mpi_comm_world,mpi_status_ignore,ierr)
    incloud(1:nzonal,1:nmerid,:) = ibuf3(:,:,:)
    deallocate(ibuf3)

  end subroutine recv_frac_hourly_int_arrays

 !*****************************************
  subroutine worker_send_toprint(toprint)
    implicit none
    real(kind=4),dimension(nzonal,nmerid,nverti) :: toprint

    integer :: ierr
    real(kind=4),allocatable,dimension(:,:,:) :: buf3

    allocate(buf3(nzonal,nmerid,nverti))
    buf3(:,:,:) = toprint(1:nzonal,1:nmerid,:)
    call mpi_send(buf3, nzonal*nmerid*nverti, mpi_real,0, iar_toprint, mpi_comm_world,ierr)
    deallocate(buf3)
    call mpi_barrier(mpi_comm_world,ierr)

  end subroutine worker_send_toprint

  !*****************************************
  subroutine worker_send_locvalues
    implicit none
    integer :: ierr
    real(kind=8),allocatable,dimension(:,:)   :: dbuf2
    real(kind=8),allocatable,dimension(:,:,:) :: dbuf3

    allocate(dbuf3(nzonal,nmerid,nverti))

    dbuf3(:,:,:) = winvloc(1:nzonal,1:nmerid,:)
    call mpi_send(dbuf3, nzonal*nmerid*nverti, mpi_double_precision,0, iar_winvloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = winxloc(1:nzonal,1:nmerid,:)
    call mpi_send(dbuf3, nzonal*nmerid*nverti, mpi_double_precision,0, iar_winxloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = sphuloc(1:nzonal,1:nmerid,:)
    call mpi_send(dbuf3, nzonal*nmerid*nverti, mpi_double_precision,0, iar_sphuloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = temploc(1:nzonal,1:nmerid,:)
    call mpi_send(dbuf3, nzonal*nmerid*nverti, mpi_double_precision,0, iar_temploc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = clwcloc(1:nzonal,1:nmerid,:)
    call mpi_send(dbuf3, nzonal*nmerid*nverti, mpi_double_precision,0, iar_clwcloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = airmloc(1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3,   nzonal*nmerid*nverti, mpi_double_precision,0, iar_airmloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = winzloc(1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3,   nzonal*nmerid*nverti, mpi_double_precision,0, iar_winzloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = winmloc(1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3,   nzonal*nmerid*nverti, mpi_double_precision,0, iar_winmloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = thlayloc(1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3,   nzonal*nmerid*nverti, mpi_double_precision,0, iar_thlayloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = hlayloc(1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3,   nzonal*nmerid*nverti, mpi_double_precision,0, iar_hlayloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = presloc(1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3, nzonal*nmerid*nverti, mpi_double_precision,0, iar_presloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = kzzzloc(1:nzonal,1:nmerid,1:nverti)
    call mpi_send(dbuf3, nzonal*nmerid*nverti, mpi_double_precision,0, iar_kzzzloc, mpi_comm_world,ierr)

    deallocate(dbuf3)


    allocate(dbuf3(nspec,nzonal,nmerid))

    dbuf3(:,:,:) = depoloc(:,1:nzonal,1:nmerid)
    call mpi_send(dbuf3, nspec*nzonal*nmerid, mpi_double_precision,0, iar_depoloc, mpi_comm_world,ierr)

    dbuf3(:,:,:) = drydep(:,1:nzonal,1:nmerid)
    call mpi_send(dbuf3,  nspec*nzonal*nmerid, mpi_double_precision,0, iar_drydep,  mpi_comm_world,ierr)

    dbuf3(:,:,:) = wetdep(:,1:nzonal,1:nmerid)
    call mpi_send(dbuf3,  nspec*nzonal*nmerid, mpi_double_precision,0, iar_wetdep,  mpi_comm_world,ierr)

    deallocate(dbuf3)


    allocate(dbuf2(nzonal,nmerid))

    dbuf2(:,:) = hghtloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_hghtloc, mpi_comm_world,ierr)

    dbuf2(:,:) = atteloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_atteloc, mpi_comm_world,ierr)

    dbuf2(:,:) = zeniloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_zeniloc, mpi_comm_world,ierr)

    dbuf2(:,:) = tem2loc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_tem2loc, mpi_comm_world,ierr)

    dbuf2(:,:) = ustaloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_ustaloc, mpi_comm_world,ierr)

    dbuf2(:,:) = aerrloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_aerrloc, mpi_comm_world,ierr)

    dbuf2(:,:) = obukloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_obukloc, mpi_comm_world,ierr)

    dbuf2(:,:) = wstaloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_wstaloc, mpi_comm_world,ierr)

    dbuf2(:,:) = topcloc(1:nzonal,1:nmerid)
    call mpi_send(dbuf2, nzonal*nmerid, mpi_double_precision,0, iar_topcloc, mpi_comm_world,ierr)

    deallocate(dbuf2)


  end subroutine worker_send_locvalues

  !******************************************
  subroutine worker_send_locvalues_tl
    implicit none
    integer :: ierr
    real(kind=8),allocatable,dimension(:,:,:) :: dbuf3

    allocate(dbuf3(nspec,nzonal,nmerid))
    dbuf3(:,:,:) = drydep_tl(:,1:nzonal,1:nmerid)
    call mpi_send(dbuf3,  nspec*nzonal*nmerid, mpi_double_precision,0, iar_drydep_tl,  mpi_comm_world,ierr)
    dbuf3(:,:,:) = wetdep_tl(:,1:nzonal,1:nmerid)
    call mpi_send(dbuf3,  nspec*nzonal*nmerid, mpi_double_precision,0, iar_wetdep_tl,  mpi_comm_world,ierr)
    deallocate(dbuf3)
  end subroutine worker_send_locvalues_tl

end module worker_message_subs
