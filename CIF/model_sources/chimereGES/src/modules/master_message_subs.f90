module master_message_subs
  
  use chimere_common
  use wholedomain_common
  use message_defs
  use chimere_consts
  
  implicit none
  include 'mpif.h'

contains
  
  !*****************************************
  subroutine init_mpi
    
    integer :: ierr,color,size,wrank
    
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    
    if(rank==0) then
      color=0
    else
      color=1
    end if
    call mpi_comm_split(mpi_comm_world,color,rank,wrk_comm,ierr)
    
    allocate(int_params(nb_int_params))
    allocate(dbl_params(nb_dbl_params))
    
    ! definition of derived types
    call mpi_type_contiguous(            &
            nb_int_params,                  &
            mpi_integer,                    &
            int_params_typ,                 &
            ierr                            &
            )
    call mpi_type_commit(int_params_typ,ierr)
    call mpi_type_contiguous(            &
            nb_dbl_params,                  &
            mpi_double_precision,           &
            dbl_params_typ,                 &
            ierr                            &
            )
    call mpi_type_commit(dbl_params_typ,ierr)
  
  
  end subroutine init_mpi
  
  
  !*****************************************
  subroutine master_init_mpi
    
    integer :: inter_comm
    integer ::ierr
    integer :: i,j,ip
    integer :: jj,kk
    
    allocate(dom(ndoms))
    ip=0
    do j=1,nmdoms
      do i=1,nzdoms
        ip=ip+1
        dom(ip)%i = i
        dom(ip)%j = j
        if(i==1) dom(ip)%izstart = 1
        if (i>1) dom(ip)%izstart = dom((j-1)*nzdoms+i-1)%izend + 1
        if (i<nzdoms) then
          dom(ip)%izend   = dom(ip)%izstart + nzonal_domain/nzdoms -1
        else
          dom(ip)%izend   = nzonal_domain
        end if
        dom(ip)%nzcount = dom(ip)%izend - dom(ip)%izstart + 1
        
        if(j==1) dom(ip)%imstart = 1
        if (j>1) dom(ip)%imstart = dom((j-2)*nzdoms+1)%imend + 1
        if (j<nmdoms) then
          dom(ip)%imend   = dom(ip)%imstart + nmerid_domain/nmdoms -1
        else
          dom(ip)%imend   = nmerid_domain
        end if
        dom(ip)%nmcount = dom(ip)%imend - dom(ip)%imstart + 1
      end do
    end do
    
    ! Optimisation along nzonal_domain
    jj = 0
    do while (dom(nzdoms)%nzcount > dom(1)%nzcount) ! size imbalance
      jj = jj + 1
      do kk = jj,nzdoms-1
        dom(kk)%izend = dom(kk)%izend + 1
        dom(kk)%nzcount = dom(kk)%izend - dom(kk)%izstart + 1
        dom(kk+1)%izstart   = dom(kk+1)%izstart   +1
      end do
      dom(nzdoms)%nzcount = dom(nzdoms)%izend - dom(nzdoms)%izstart + 1
    end do
    ! clone to other domains
    ip = 0
    do j=1,nmdoms-1
      do i=1,nzdoms
        ip=ip+1
        dom(ip+nzdoms)%izstart = dom(ip)%izstart
        dom(ip+nzdoms)%izend   = dom(ip)%izend
        dom(ip+nzdoms)%nzcount = dom(ip)%nzcount
      end do
    end do
    ! Optimisation along nzmerid
    jj = 0
    do while (dom(1+(nmdoms-1)*nzdoms)%nmcount > dom(1)%nmcount) ! size unbalance
      jj = jj + 1
      do kk = jj,nmdoms-1
        dom(1+(kk-1)*nzdoms)%imend   = dom(1+(kk-1)*nzdoms)%imend + 1
        dom(1+(kk-1)*nzdoms)%nmcount = dom(1+(kk-1)*nzdoms)%imend - dom(1+(kk-1)*nzdoms)%imstart + 1
        dom(1+kk*nzdoms)%imstart     = dom(1+kk*nzdoms)%imstart   +1
      end do
      dom(1+(nmdoms-1)*nzdoms)%nmcount = dom(1+(nmdoms-1)*nzdoms)%imend - dom(1+(nmdoms-1)*nzdoms)%imstart + 1
    end do
    ! clone to other domains
    ip = 0
    do j=1,nmdoms
      do i=2,nzdoms
        dom(i+(j-1)*nzdoms)%imstart = dom(1+(j-1)*nzdoms)%imstart
        dom(i+(j-1)*nzdoms)%imend   = dom(1+(j-1)*nzdoms)%imend
        dom(i+(j-1)*nzdoms)%nmcount = dom(1+(j-1)*nzdoms)%nmcount
      end do
    end do
    
    nzonalmax=0
    nmeridmax=0
    do i=1,ndoms
      if (dom(i)%nzcount .gt. nzonalmax) nzonalmax=    dom(i)%nzcount
      if (dom(i)%nmcount .gt. nmeridmax) nmeridmax=    dom(i)%nmcount
    end do
    
    print *
    print *,'     +++ CHIMERE RUNNING IN PARALLEL MODE +++'
    print *,'             MPI SUB-DOMAINS :'
    write(*,'("rank  izstart  izend  nzcount  imstart imend  nmcount     i       j")')
    write(*,'("--------------------------------------------------------------------")')
    ip=0
    do j=1,nmdoms
      do i=1,nzdoms
        ip=ip+1
        
        write (*,'(1x,i3,5x,i3,5x,i3,5x,i3,5x,i3,5x,i3,5x,i3,5x,i3,5x,i3)') ip,dom(ip)
      end do
    end do
    print *
  
  end subroutine master_init_mpi
  
  
  !*****************************************
  subroutine master_send_once
    integer :: ip
    integer ::  ierr
    
    call wcopy_dbl_params
    do ip=1,ndoms
      int_params(ip_nzonal)      = dom(ip)%nzcount
      int_params(ip_nmerid)      = dom(ip)%nmcount
      int_params(ip_dom_i)       = dom(ip)%i
      int_params(ip_dom_j)       = dom(ip)%j
      int_params(ip_wrk_comm)    = wrk_comm
      int_params(ip_imstart)	= dom(ip)%imstart
      int_params(ip_imend)	= dom(ip)%imend
      int_params(ip_izstart)	= dom(ip)%izstart
      int_params(ip_izend)	= dom(ip)%izend
      int_params(ip_nspec) = nspec
      int_params(ip_nspectot) = nspectot
      int_params(ip_nreac) = nreac
      int_params(ip_nfam) = nfam
      int_params(ip_ndepo) = ndepo
      int_params(ip_nreactamax) = nreactamax
      int_params(ip_ntabuzen)=ntabuzen
      int_params(ip_nlevphot)=nlevphot
      int_params(ip_nphot)= nphot
      int_params(ip_nhourrun) = nhourrun
      int_params(ip_ntemps) = ntemps
      int_params(ip_ntabmax) = ntabmax
      int_params(ip_nverti) = nverti
      int_params(ip_nphour_ref) = nphour_ref
      int_params(ip_ichemstep) = ichemstep
      int_params(ip_usechemistry) = usechemistry
      int_params(ip_usedepos) = usedepos
      int_params(ip_usewetdepos) = usewetdepos
      int_params(ip_useemissions) = useemissions
      int_params(ip_usetransmix) = usetransmix
      int_params(ip_useabsclipconc)=useabsclipconc
      int_params(ip_dryairout)=dryairout
      int_params(ip_nvegtype)=nvegtype
      int_params(ip_nlduse)=nlduse
      int_params(ip_nsaveconcs)=nsaveconcs
      int_params(ip_noutspec)=noutspec
      int_params(ip_nzdoms)=nzdoms
      int_params(ip_nmdoms)=nmdoms
      int_params(ip_nzonalmax)=nzonalmax
      int_params(ip_nmeridmax)=nmeridmax
      int_params(ip_ideepconv)=ideepconv
      int_params(ip_nitgs)=nitgs
      int_params(ip_nitgssu)=nitgssu
      int_params(ip_optemisb)=optemisb
      int_params(ip_nemisa)=nemisa
      int_params(ip_nemisb)=nemisb
      int_params(ip_ndep)=ndep
      int_params(ip_nlevemis)=nlevemis
      int_params(ip_ihoursu)=ihoursu
      int_params(ip_fwd)=fwd
      int_params(ip_tl)=tl
      int_params(ip_ad)=ad
      int_params(ip_hpulse)=hpulse
      call mpi_send(int_params,1,int_params_typ,ip,ias_int_params,mpi_comm_world,ierr)
      call mpi_send(dbl_params,1,dbl_params_typ,ip,ias_dbl_params,mpi_comm_world,ierr)
      call send_int_arrays(ip,dom(ip))
      call send_char_arrays(ip)
      call send_real_arrays(ip)
      call send_dbl_arrays(ip,dom(ip))
      call send_hourly_real_arrays(ip,dom(ip))
      if (tl.eq.1) call send_hourly_real_arrays_tl(ip,dom(ip))
      call mpi_send(nphour,    1,mpi_integer,         ip, ias_nphour, mpi_comm_world,ierr)
      call mpi_send(dtr,   1,mpi_double_precision,   ip, ias_dtr,    mpi_comm_world,ierr)
      call mpi_send(dtr2,  1,mpi_double_precision,   ip, ias_dtr2,   mpi_comm_world,ierr)
    end do
    call send_output_species
  
  end subroutine master_send_once
  
  
  !*****************************************
  subroutine master_recv_locvalues
    integer :: ip
    
    do ip=1,ndoms
      call recv_locvalues(ip,dom(ip))
      if (tl.eq.1) call recv_locvalues_tl(ip,dom(ip))
    end do
  
  end subroutine master_recv_locvalues
  
  
  !*****************************************
  subroutine master_send_hourly
    integer :: ip,ierr
    
    do ip=1,ndoms
      call send_hourly_real_arrays(ip,dom(ip))
      if (tl.eq.1) call send_hourly_real_arrays_tl(ip,dom(ip))
      call mpi_send(nphour,1,mpi_integer,            ip, ias_nphour, mpi_comm_world,ierr)
      call mpi_send(dtr,   1,mpi_double_precision,   ip, ias_dtr,    mpi_comm_world,ierr)
      call mpi_send(dtr2,  1,mpi_double_precision,   ip, ias_dtr2,   mpi_comm_world,ierr)
    enddo
  
  end subroutine master_send_hourly
  
  
  !*****************************************
  subroutine master_send_frac_hourly
    integer :: ip
    
    do ip=1,ndoms
      call send_frac_hourly_int_arrays(ip,dom(ip))
    end do
  
  end subroutine master_send_frac_hourly
  
  !*****************************************
  subroutine master_recv_toprint(toprint)
    implicit none
    real(kind=4),dimension(nzonal_domain,nmerid_domain,nverti) :: toprint
    integer :: ip,ierr
    
    do ip=1,ndoms
      call recv_toprint(ip,dom(ip),toprint)
    end do
    call mpi_barrier(mpi_comm_world,ierr)
  
  end subroutine master_recv_toprint
  
  
  !*****************************************
  subroutine master_recv_conc
    integer :: ip
    integer :: ierr
    integer,dimension(mpi_status_size) :: status
    real(kind=8),allocatable,dimension(:,:,:) :: dbuf3
    
    call mpi_barrier(mpi_comm_world,ierr)
    do ip=1,ndoms
      allocate(dbuf3(dom(ip)%nzcount,dom(ip)%nmcount,nverti+1))
      ! conc
      call mpi_recv(                                            &
              dbuf3,                                               &
              dom(ip)%nzcount*dom(ip)%nmcount*(nverti+1), &
              mpi_double_precision,ip,                             &
              iar_conc,                                            &
              mpi_comm_world,status,ierr                               &
              )
      conc(dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,:) = dbuf3
      deallocate(dbuf3)
      ! conco
      allocate(dbuf3(dom(ip)%nzcount,dom(ip)%nmcount,nverti))
      call mpi_recv(                                            &
              dbuf3,                                               &
              dom(ip)%nzcount*dom(ip)%nmcount*nverti, &
              mpi_double_precision,ip,                             &
              iar_conco,                                            &
              mpi_comm_world,status,ierr                               &
              )
      conco(dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,1:nverti) = dbuf3
      deallocate(dbuf3)
    end do
  
  end subroutine master_recv_conc
  !*****************************************
  subroutine master_recv_conc_tl
    integer :: ip
    integer :: ierr
    integer,dimension(mpi_status_size) :: status
    real(kind=8),allocatable,dimension(:,:,:) :: dbuf6
    
    do ip=1,ndoms
      allocate(dbuf6(dom(ip)%nzcount,dom(ip)%nmcount,nverti+1))
      ! conc
      call mpi_recv(                                            &
              dbuf6,                                               &
              dom(ip)%nzcount*dom(ip)%nmcount*(nverti+1), &
              mpi_double_precision,ip,                             &
              iar_conc_tl,                                            &
              mpi_comm_world,status,ierr                               &
              )
      conc_tl(dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,:) = dbuf6
      deallocate(dbuf6)
      allocate(dbuf6(dom(ip)%nzcount,dom(ip)%nmcount,nverti))
      ! conco
      call mpi_recv(                                            &
              dbuf6,                                               &
              dom(ip)%nzcount*dom(ip)%nmcount*nverti, &
              mpi_double_precision,ip,                             &
              iar_conco_tl,                                            &
              mpi_comm_world,status,ierr                               &
              )
      conco_tl(dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,:) = dbuf6
      deallocate(dbuf6)
    end do
  
  end subroutine master_recv_conc_tl
  !*****************************************
  subroutine master_send_ns(ns)
    implicit none
    integer :: ns
    integer :: ip, ierr
    
    do ip=1,ndoms
      call mpi_send(ns,1,mpi_integer,ip, ias_ns, mpi_comm_world,ierr)
    end do
  end subroutine master_send_ns
  
  
  
  
  !*****************************************
  subroutine master_send_conc
    integer :: ip,ierr
    real(kind=8),dimension(:,:,:),allocatable :: dbuf4
    
    do ip=1,ndoms
      allocate(dbuf4(dom(ip)%nzcount+6, dom(ip)%nmcount+6, nverti+1))
      
      ! conc
      dbuf4=conc(dom(ip)%izstart-3:dom(ip)%izend+3, dom(ip)%imstart-3:dom(ip)%imend+3, :)
      call mpi_send( &
              dbuf4,                                                       &
              (dom(ip)%nzcount+6)*(dom(ip)%nmcount+6)*(nverti+1), &
              mpi_double_precision,ip,                                     &
              ias_conc,                                                    &
              mpi_comm_world,ierr)
      deallocate(dbuf4)
      ! conco
      allocate(dbuf4(dom(ip)%nzcount, dom(ip)%nmcount, nverti))
      dbuf4=conco(dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,1:nverti)
      call mpi_send( &
              dbuf4,                                                       &
              dom(ip)%nzcount*dom(ip)%nmcount*nverti,                      &
              mpi_double_precision,ip,                                     &
              ias_conco,                                                   &
              mpi_comm_world,ierr)
      deallocate(dbuf4)
    
    end do
  
  end subroutine master_send_conc
  !*****************************************
  subroutine master_send_conc_tl
    integer :: ip,ierr
    real(kind=8),dimension(:,:,:),allocatable :: dbuf6
    
    do ip=1,ndoms
      allocate(dbuf6(dom(ip)%nzcount+6,dom(ip)%nmcount+6,nverti+1))
      ! conc
      dbuf6=conc_tl(dom(ip)%izstart-3:dom(ip)%izend+3,dom(ip)%imstart-3:dom(ip)%imend+3,1:nverti+1)
      call mpi_send( &
              dbuf6,  &
              (dom(ip)%nzcount+6)*(dom(ip)%nmcount+6)*(nverti+1), &
              mpi_double_precision,                               &
              ip,ias_conc_tl,                                     &
              mpi_comm_world,ierr)
      deallocate(dbuf6)
      ! conco
      allocate(dbuf6(dom(ip)%nzcount, dom(ip)%nmcount, nverti))
      dbuf6=conco_tl(dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,1:nverti)
      call mpi_send( &
              dbuf6,                                                       &
              dom(ip)%nzcount*dom(ip)%nmcount*nverti,                      &
              mpi_double_precision,ip,                                     &
              ias_conco_tl,                                                   &
              mpi_comm_world,ierr)
      deallocate(dbuf6)
    enddo
  
  end subroutine master_send_conc_tl
  !*****************************************
  subroutine master_mpi_finalize
    integer :: ierr
    
    call master_deallocall
    call mpi_barrier(mpi_comm_world,ierr)
    call mpi_finalize(ierr)
  
  end subroutine master_mpi_finalize
  
  
  !*****************************************
  subroutine master_send_conc_bounds
    
    integer :: ip,ierr
    integer :: izstart,izend,imstart,imend,nzcount,nmcount,idom,jdom
    integer :: vl_req,vr_req,hl_req,hu_req,top_req
    real(kind=8),allocatable,dimension(:,:,:) :: hlbuf,hubuf,vlbuf,vrbuf
    real(kind=8),allocatable,dimension(:,:)   :: topbuf
    
    do ip=1,ndoms
      
      izstart = dom(ip)%izstart
      izend   = dom(ip)%izend
      imstart = dom(ip)%imstart
      imend   = dom(ip)%imend
      nzcount = dom(ip)%nzcount
      nmcount = dom(ip)%nmcount
      idom    = dom(ip)%i
      jdom    = dom(ip)%j
      
      if(idom==1) then
        !print*,'to ',ip,' send left boundary',nmcount+6,nverti+1
        allocate(vlbuf(3, nmcount+6, nverti+1))
        vlbuf(:,:,:)=conc(izstart-3:izstart-1, imstart-3:imend+3, :)
        call mpi_issend( &
                vlbuf,                                           &
                3*(nmcount+6)*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds,                                 &
                mpi_comm_world,                                  &
                vl_req,                                          &
                ierr)
      end if
      
      if(idom==nzdoms) then
        !print*,'to ',ip,' send right boundary'
        allocate(vrbuf(3, nmcount+6, nverti+1))
        vrbuf(:,:,:)=conc(izend+1:izend+3, imstart-3:imend+3, :)
        call mpi_issend( &
                vrbuf,                                           &
                3*(nmcount+6)*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds,                                 &
                mpi_comm_world,                                  &
                vr_req,                                          &
                ierr)
      end if
      
      if(jdom==1) then
        !print*,'to ',ip,' send lower boundary'
        allocate(hlbuf(nzcount+6, 3, nverti+1))
        hlbuf(:,:,:)=conc(izstart-3:izend+3, imstart-3:imstart-1, :)
        call mpi_issend( &
                hlbuf,                                           &
                (nzcount+6)*3*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds,                                 &
                mpi_comm_world,                                  &
                hl_req,                                          &
                ierr)
      end if
      
      if(jdom==nmdoms) then
        !print*,'to ',ip,' send upper boundary'
        allocate(hubuf(nzcount+6, 3, nverti+1))
        hubuf(:,:,:)=conc(izstart-3:izend+3, imend+1:imend+3, :)
        call mpi_issend( &
                hubuf,                                           &
                (nzcount+6)*3*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds,                                 &
                mpi_comm_world,                                  &
                hu_req,                                          &
                ierr)
      end if
      
      ! print*,'to ',ip,' send top boundary'
      allocate(topbuf(dom(ip)%nzcount, dom(ip)%nmcount))
      topbuf(:,:)=conc(izstart:izend, imstart:imend, nverti+1)
      call mpi_issend( &
              topbuf,                                          &
              nzcount*nmcount,                        &
              mpi_double_precision,ip,                         &
              ias_conc_top,                                    &
              mpi_comm_world,                                  &
              top_req,                                         &
              ierr)
      if(idom==1) then
        call mpi_wait(vl_req,mpi_status_ignore,ierr)
        deallocate(vlbuf)
      end if
      
      if(idom==nzdoms) then
        call mpi_wait(vr_req,mpi_status_ignore,ierr)
        deallocate(vrbuf)
      end if
      
      if(jdom==1) then
        call mpi_wait(hl_req,mpi_status_ignore,ierr)
        deallocate(hlbuf)
      end if
      
      if(jdom==nmdoms) then
        call mpi_wait(hu_req,mpi_status_ignore,ierr)
        deallocate(hubuf)
      end if
      
      call mpi_wait(top_req,mpi_status_ignore,ierr)
      deallocate(topbuf)
    
    end do
  
  end subroutine master_send_conc_bounds
  !*****************************************
  subroutine amaster_send_aconc_bounds
    
    integer :: ip,ierr
    integer :: izstart,izend,imstart,imend,nzcount,nmcount,idom,jdom
    integer :: vl_req,vr_req,hl_req,hu_req,top_req
    real(kind=8),allocatable,dimension(:,:,:,:) :: hlbuf,hubuf,vlbuf,vrbuf
    real(kind=8),allocatable,dimension(:,:,:)   :: topbuf
    
    do ip=1,ndoms
      
      izstart = dom(ip)%izstart
      izend   = dom(ip)%izend
      imstart = dom(ip)%imstart
      imend   = dom(ip)%imend
      nzcount = dom(ip)%nzcount
      nmcount = dom(ip)%nmcount
      idom    = dom(ip)%i
      jdom    = dom(ip)%j
      
      
      if(idom==1) then
        ! send left boundary
        allocate(vlbuf(nspectot, 3, nmcount+6, nverti+1))
        vlbuf(:,:,:,:)=aconc(:, izstart-3:izstart-1, imstart-3:imend+3, :)
        call mpi_issend( &
                vlbuf,                                           &
                nspectot*3*(nmcount+6)*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_aconc_bounds,                                 &
                mpi_comm_world,                                  &
                vl_req,                                          &
                ierr)
        aconc(:, izstart-3:izstart-1, imstart-3:imend+3, :)=0.d0
      end if
      
      if(idom==nzdoms) then
        ! send right boundary
        allocate(vrbuf(nspectot, 3, nmcount+6, nverti+1))
        vrbuf(:,:,:,:)=aconc(:, izend+1:izend+3, imstart-3:imend+3, :)
        call mpi_issend( &
                vrbuf,                                           &
                nspectot*3*(nmcount+6)*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_aconc_bounds,                                 &
                mpi_comm_world,                                  &
                vr_req,                                          &
                ierr)
        aconc(:, izend+1:izend+3, imstart-3:imend+3, :)=0.d0
      end if
      
      if(jdom==1) then
        ! send lower boundary
        allocate(hlbuf(nspectot, nzcount+6, 3, nverti+1))
        hlbuf(:,:,:,:)=aconc(:, izstart-3:izend+3, imstart-3:imstart-1, :)
        call mpi_issend( &
                hlbuf,                                           &
                nspectot*(nzcount+6)*3*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_aconc_bounds,                                 &
                mpi_comm_world,                                  &
                hl_req,                                          &
                ierr)
        aconc(:, izstart-3:izend+3, imstart-3:imstart-1, :)=0.d0
      end if
      
      if(jdom==nmdoms) then
        ! send upper boundary
        allocate(hubuf(nspectot, nzcount+6, 3, nverti+1))
        hubuf(:,:,:,:)=aconc(:, izstart-3:izend+3, imend+1:imend+3, :)
        call mpi_issend( &
                hubuf,                                           &
                nspectot*(nzcount+6)*3*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_aconc_bounds,                                 &
                mpi_comm_world,                                  &
                hu_req,                                          &
                ierr)
        aconc(:, izstart-3:izend+3, imend+1:imend+3, :)=0.d0
      end if
      
      ! send top boundary
      allocate(topbuf(nspectot, dom(ip)%nzcount, dom(ip)%nmcount))
      topbuf(:,:,:)=aconc(:, izstart:izend, imstart:imend, nverti+1)
      call mpi_issend( &
              topbuf,                                          &
              nspectot*nzcount*nmcount,                        &
              mpi_double_precision,ip,                         &
              ias_aconc_top,                                    &
              mpi_comm_world,                                  &
              top_req,                                         &
              ierr)
      
      if(idom==1) then
        call mpi_wait(vl_req,mpi_status_ignore,ierr)
        deallocate(vlbuf)
      end if
      
      if(idom==nzdoms) then
        call mpi_wait(vr_req,mpi_status_ignore,ierr)
        deallocate(vrbuf)
      end if
      
      if(jdom==1) then
        call mpi_wait(hl_req,mpi_status_ignore,ierr)
        deallocate(hlbuf)
      end if
      
      if(jdom==nmdoms) then
        call mpi_wait(hu_req,mpi_status_ignore,ierr)
        deallocate(hubuf)
      end if
      
      call mpi_wait(top_req,mpi_status_ignore,ierr)
      deallocate(topbuf)
      aconc(:, izstart:izend, imstart:imend, nverti+1)=0.d0
    
    end do
  
  end subroutine amaster_send_aconc_bounds
  !*****************************************
  subroutine master_send_conc_bounds_tl
    
    integer :: ip,ierr
    integer :: izstart,izend,imstart,imend,nzcount,nmcount,idom,jdom
    integer :: vl_req_tl,vr_req_tl,hl_req_tl,hu_req_tl,top_req_tl
    real(kind=8),allocatable,dimension(:,:,:) :: hlbuf,hubuf,vlbuf,vrbuf
    real(kind=8),allocatable,dimension(:,:)   :: topbuf
    
    do ip=1,ndoms
      
      izstart = dom(ip)%izstart
      izend   = dom(ip)%izend
      imstart = dom(ip)%imstart
      imend   = dom(ip)%imend
      nzcount = dom(ip)%nzcount
      nmcount = dom(ip)%nmcount
      idom    = dom(ip)%i
      jdom    = dom(ip)%j
      
      
      if(idom==1) then
        ! send left boundary
        allocate(vlbuf(3, nmcount+6, nverti+1))
        vlbuf(:,:,:)=conc_tl(izstart-3:izstart-1, imstart-3:imend+3, :)
        call mpi_issend( &
                vlbuf,                                           &
                3*(nmcount+6)*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds_tl,                                 &
                mpi_comm_world,                                  &
                vl_req_tl,                                          &
                ierr)
      end if
      
      if(idom==nzdoms) then
        ! send right boundary
        allocate(vrbuf(3, nmcount+6, nverti+1))
        vrbuf(:,:,:)=conc_tl(izend+1:izend+3, imstart-3:imend+3, :)
        call mpi_issend( &
                vrbuf,                                           &
                3*(nmcount+6)*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds_tl,                                 &
                mpi_comm_world,                                  &
                vr_req_tl,                                          &
                ierr)
      end if
      
      if(jdom==1) then
        ! send lower boundary
        allocate(hlbuf(nzcount+6, 3, nverti+1))
        hlbuf(:,:,:)=conc_tl(izstart-3:izend+3, imstart-3:imstart-1, :)
        call mpi_issend( &
                hlbuf,                                           &
                (nzcount+6)*3*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds_tl,                                 &
                mpi_comm_world,                                  &
                hl_req_tl,                                          &
                ierr)
      end if
      
      if(jdom==nmdoms) then
        ! send upper boundary
        allocate(hubuf(nzcount+6, 3, nverti+1))
        hubuf(:,:,:)=conc_tl(izstart-3:izend+3, imend+1:imend+3, :)
        call mpi_issend( &
                hubuf,                                           &
                (nzcount+6)*3*(nverti+1),               &
                mpi_double_precision,ip,                         &
                ias_conc_bounds_tl,                                 &
                mpi_comm_world,                                  &
                hu_req_tl,                                          &
                ierr)
      end if
      
      ! send top boundary
      allocate(topbuf(dom(ip)%nzcount, dom(ip)%nmcount))
      topbuf(:,:)=conc_tl(izstart:izend, imstart:imend, nverti+1)
      call mpi_issend( &
              topbuf,                                          &
              nzcount*nmcount,                        &
              mpi_double_precision,ip,                         &
              ias_conc_top_tl,                                    &
              mpi_comm_world,                                  &
              top_req_tl,                                         &
              ierr)
      
      if(idom==1) then
        call mpi_wait(vl_req_tl,mpi_status_ignore,ierr)
        deallocate(vlbuf)
      end if
      
      if(idom==nzdoms) then
        call mpi_wait(vr_req_tl,mpi_status_ignore,ierr)
        deallocate(vrbuf)
      end if
      
      if(jdom==1) then
        call mpi_wait(hl_req_tl,mpi_status_ignore,ierr)
        deallocate(hlbuf)
      end if
      
      if(jdom==nmdoms) then
        call mpi_wait(hu_req_tl,mpi_status_ignore,ierr)
        deallocate(hubuf)
      end if
      
      call mpi_wait(top_req_tl,mpi_status_ignore,ierr)
      deallocate(topbuf)
    
    
    end do
  
  end subroutine master_send_conc_bounds_tl
  !*****************************************
  subroutine wcopy_dbl_params
    ! This routine gathers well-known chimere scalars
    ! into the array int_params 
    
    ! It use a macro to shorten code writing and reduce risk of typo errors
    ! The macro generates lines like :
    ! dbl_params(ip_factorl)=factorl
    !
    ! ip_factorl beeing an arbitrary but unique index defined in module message_defs
    
#define W_DBL_PARAMS(par)             dbl_params(ip_##par)= par
    implicit none
  
    W_DBL_PARAMS(soltim)
    W_DBL_PARAMS(djul)
    W_DBL_PARAMS(clipconc)
    W_DBL_PARAMS(psurf)
    
  end subroutine wcopy_dbl_params
  
  
  !*****************************************
  subroutine send_int_arrays(ip,dom)
    implicit none
    integer :: ip,izstart,izend,nzcount,imstart,imend,nmcount
    integer :: ierr
    type(dom_type) :: dom
    integer,dimension(:),allocatable :: species_varid,species_transp,species_transpv,species_bounddry
    integer,dimension(:,:,:),allocatable :: ibuf3
    
    call mpi_send(inemisa,    nspec,                mpi_integer, ip,ias_inemisa,   mpi_comm_world,ierr)
    call mpi_send(inemisb,    nspec,                mpi_integer, ip,ias_inemisb,   mpi_comm_world,ierr)
    call mpi_send(kreacl,     nspectot,             mpi_integer, ip,ias_kreacl,    mpi_comm_world,ierr)
    call mpi_send(kreacp,     nspectot,             mpi_integer, ip,ias_kreacp,    mpi_comm_world,ierr)
    call mpi_send(ireacl,     nspectot*nreac,       mpi_integer, ip,ias_ireacl,    mpi_comm_world,ierr)
    call mpi_send(ireacp,     nspectot*nreac,       mpi_integer, ip,ias_ireacp,    mpi_comm_world,ierr)
    call mpi_send(nreactants, nreac,                mpi_integer, ip,ias_nreactants,mpi_comm_world,ierr)
    call mpi_send(irctt,      nreac*nreactamax,     mpi_integer, ip,ias_irctt,     mpi_comm_world,ierr)
    call mpi_send(indepo,     nspec+1,              mpi_integer, ip,ias_indepo,    mpi_comm_world,ierr)
    call mpi_send(inwetd,     nspec*2,              mpi_integer, ip,ias_inwetd,    mpi_comm_world,ierr)
    call mpi_send(imonth,     nhourrun+1,           mpi_integer, ip,ias_imonth,    mpi_comm_world,ierr)
    call mpi_send(ityperate,  nreac,                mpi_integer, ip,ias_ityperate, mpi_comm_world,ierr)
    call mpi_send(nwetd,      2,                    mpi_integer, ip,ias_nwetd,     mpi_comm_world,ierr)
    call mpi_send(iphoto,     nreac,                mpi_integer, ip,ias_iphoto,    mpi_comm_world,ierr)
    call mpi_send(ispecemip , 1000,                 mpi_integer, ip,ias_ispecemip, mpi_comm_world,ierr)
    call mpi_send(nelem      ,nfam,                 mpi_integer, ip,ias_nelem,     mpi_comm_world,ierr)
    call mpi_send(ifam       ,nfam*nspectot*10,     mpi_integer, ip,ias_ifam,      mpi_comm_world,ierr)
    
    
    nzcount = dom%nzcount
    nmcount = dom%nmcount
    izstart = dom%izstart
    izend   = dom%izend
    imstart = dom%imstart
    imend   = dom%imend
    
    allocate(species_varid(nspectot))
    species_varid(:)=species(:)%varid
    call mpi_send(species_varid, nspectot, mpi_integer,ip,ias_species_varid,mpi_comm_world,ierr)
    deallocate(species_varid)
    
    allocate(species_transp(nspectot))
    species_transp(:)=species(:)%transp
    call mpi_send(species_transp, nspectot, mpi_integer,ip,ias_species_transp,mpi_comm_world,ierr)
    deallocate(species_transp)
    
    allocate(species_transpv(nspectot))
    species_transpv(:)=species(:)%transpv
    call mpi_send(species_transpv, nspectot, mpi_integer,ip,ias_species_transpv,mpi_comm_world,ierr)
    deallocate(species_transpv)
    
    allocate(species_bounddry(nspectot))
    species_bounddry(:)=species(:)%bounddry
    call mpi_send(species_bounddry, nspectot, mpi_integer,ip,ias_species_bounddry,mpi_comm_world,ierr)
    deallocate(species_bounddry)
    
    call mpi_send(ihour, nhourrun+1,  mpi_integer, ip,ias_ihour, mpi_comm_world,ierr)
    
    call mpi_send(idtyp ,nhourrun+1,  mpi_integer, ip,ias_idtyp, mpi_comm_world,ierr)
  
  end subroutine send_int_arrays
  
  
  !*****************************************
  subroutine send_char_arrays(ip)
    implicit none
    integer :: ip,ierr
    character(len=16),dimension(:),allocatable :: species_name
    
    allocate(species_name(nspectot))
    species_name(:)=species(:)%name
    call mpi_send(species_name, 16*nspectot, mpi_character,ip,ias_species_name,mpi_comm_world,ierr)
    deallocate(species_name)
  
  end subroutine send_char_arrays
  
  
  !*****************************************
  subroutine send_real_arrays(ip)
    implicit none
    integer :: ip,ierr
    
    call mpi_send(altiphot,nlevphot,mpi_double_precision, ip, ias_altiphot, mpi_comm_world,ierr)
    call mpi_send(photoj,ntabuzen*nlevphot*nphot, mpi_double_precision, ip, ias_photoj, mpi_comm_world, ierr)
  
  end subroutine send_real_arrays
  
  !*****************************************
  subroutine send_output_species
    implicit none
    integer :: ip
    integer :: ierr
    real(kind=8),        dimension(nspectot) :: fspec
    integer,             dimension(nspectot) :: varid, iaddr
    character(len=splen),dimension(nspectot) :: name
    character(len=16),   dimension(nspectot) :: units
    
    fspec(:) = species(:)%fspec
    varid(:) = output_species(:)%varid
    iaddr(:) = output_species(:)%iaddr
    name(:)  = output_species(:)%name
    units(:) = output_species(:)%units
    
    do ip=1,ndoms
      call mpi_send(fspec,nspectot,      mpi_double_precision,ip, ias_os_fscpe, mpi_comm_world,ierr)
      call mpi_send(varid,nspectot,      mpi_integer,         ip, ias_os_varid, mpi_comm_world,ierr)
      call mpi_send(iaddr,nspectot,      mpi_integer,         ip, ias_os_iaddr, mpi_comm_world,ierr)
      call mpi_send(name, splen*nspectot,mpi_character,       ip, ias_os_name,  mpi_comm_world,ierr)
      call mpi_send(units,16*nspectot,   mpi_character,       ip, ias_os_units, mpi_comm_world,ierr)
    
    end do
  
  end subroutine send_output_species
  
  
  !*****************************************
  subroutine send_dbl_arrays(ip,dom)
    type(dom_type) :: dom
    integer :: ip,izstart,izend,nzcount,imstart,imend,nmcount,ierr
    real(kind=8),dimension(:,:)    ,allocatable :: dbuf2
    real(kind=8),dimension(:,:,:)  ,allocatable :: dbuf3
    real(kind=8),dimension(:,:,:,:),allocatable :: dbuf4
    real(kind=8),dimension(:,:,:,:,:),allocatable :: dbuf5
    call mpi_send(stoi,    nspectot*nreac*ntemps,       mpi_double_precision,ip, ias_stoi,   mpi_comm_world,ierr)
    call mpi_send(deptmin, nvegtype,                    mpi_double_precision,ip,ias_deptmin, mpi_comm_world,ierr)
    call mpi_send(deptopt, nvegtype,                    mpi_double_precision,ip,ias_deptopt, mpi_comm_world,ierr)
    call mpi_send(depvpd1, nvegtype,                    mpi_double_precision,ip,ias_depvpd1, mpi_comm_world,ierr)
    call mpi_send(depvpd2, nvegtype,                    mpi_double_precision,ip,ias_depvpd2, mpi_comm_world,ierr)
    call mpi_send(depegs,  nvegtype,                    mpi_double_precision,ip,ias_depegs,  mpi_comm_world,ierr)
    call mpi_send(depsgs,  nvegtype,                    mpi_double_precision,ip,ias_depsgs,  mpi_comm_world,ierr)
    call mpi_send(depegl,  nvegtype,                    mpi_double_precision,ip,ias_depegl,  mpi_comm_world,ierr)
    call mpi_send(depsgl,  nvegtype,                    mpi_double_precision,ip,ias_depsgl,  mpi_comm_world,ierr)
    call mpi_send(deplai1, nvegtype,                    mpi_double_precision,ip,ias_deplai1, mpi_comm_world,ierr)
    call mpi_send(deplai2, nvegtype,                    mpi_double_precision,ip,ias_deplai2, mpi_comm_world,ierr)
    call mpi_send(depphe0, nvegtype,                    mpi_double_precision,ip,ias_depphe0, mpi_comm_world,ierr)
    call mpi_send(depphe1, nvegtype,                    mpi_double_precision,ip,ias_depphe1, mpi_comm_world,ierr)
    call mpi_send(depphe2, nvegtype,                    mpi_double_precision,ip,ias_depphe2, mpi_comm_world,ierr)
    call mpi_send(depalph, nvegtype,                    mpi_double_precision,ip,ias_depalph, mpi_comm_world,ierr)
    call mpi_send(gmax,    nvegtype,                    mpi_double_precision,ip,ias_gmax,    mpi_comm_world,ierr)
    call mpi_send(fmin,    nvegtype,                    mpi_double_precision,ip,ias_fmin,    mpi_comm_world,ierr)
    call mpi_send(zcanopy, nvegtype,                    mpi_double_precision,ip,ias_zcanopy, mpi_comm_world,ierr)
    call mpi_send(rgso3,   nvegtype,                    mpi_double_precision,ip,ias_rgso3,   mpi_comm_world,ierr)
    call mpi_send(so2rh,   nvegtype,                    mpi_double_precision,ip,ias_so2rh,   mpi_comm_world,ierr)
    call mpi_send(rgsso2,  nvegtype,                    mpi_double_precision,ip,ias_rgsso2,  mpi_comm_world,ierr)
    call mpi_send(factrb,  nspec,                       mpi_double_precision,ip,ias_factrb,  mpi_comm_world,ierr)
    call mpi_send(factd,   nspec,                       mpi_double_precision,ip,ias_factd,   mpi_comm_world,ierr)
    call mpi_send(rm,      nspec,                       mpi_double_precision,ip,ias_rm,      mpi_comm_world,ierr)
    call mpi_send(dhx,     nspec,                       mpi_double_precision,ip,ias_dhx,     mpi_comm_world,ierr)
    call mpi_send(df0,     nspec,                       mpi_double_precision,ip,ias_df0,     mpi_comm_world,ierr)
    call mpi_send(zetaref, ntabuzen,                 mpi_double_precision,ip,ias_zetaref, mpi_comm_world,ierr)
    call mpi_send(tabrate, ntabmax*nreac,               mpi_double_precision,ip,ias_tabrate, mpi_comm_world,ierr)
    call mpi_send(tabtemp, ntemps,                      mpi_double_precision,ip,ias_tabtemp, mpi_comm_world,ierr)
    nzcount=dom%nzcount
    nmcount=dom%nmcount
    izstart=dom%izstart
    izend  =dom%izend
    imstart=dom%imstart
    imend  =dom%imend
    
    allocate(dbuf2(nzcount+4,nmcount+4))
    dbuf2=xsize(izstart-2:izend+2,imstart-2:imend+2)
    call mpi_send(dbuf2,                                    &
            (nzcount+4)*(nmcount+4),                           &
            mpi_double_precision,ip,                           &
            ias_xsize,                                         &
            mpi_comm_world,ierr)
    dbuf2=ysize(izstart-2:izend+2,imstart-2:imend+2)
    call mpi_send(dbuf2,                                    &
            (nzcount+4)*(nmcount+4),                           &
            mpi_double_precision,ip,                           &
            ias_ysize,                                         &
            mpi_comm_world,ierr)
    deallocate(dbuf2)
    allocate(dbuf2(nzcount+2,nmcount+2))
    dbuf2=xbasx(izstart-1:izend+1,imstart-1:imend+1)
    call mpi_send(dbuf2,                                    &
            (nzcount+2)*(nmcount+2),                           &
            mpi_double_precision,ip,                           &
            ias_xbasx,                                         &
            mpi_comm_world,ierr)
    dbuf2=xbasy(izstart-1:izend+1,imstart-1:imend+1)
    call mpi_send(dbuf2,                                    &
            (nzcount+2)*(nmcount+2),                           &
            mpi_double_precision,ip,                           &
            ias_xbasy,                                         &
            mpi_comm_world,ierr)
    dbuf2=ybasx(izstart-1:izend+1,imstart-1:imend+1)
    call mpi_send(dbuf2,                                    &
            (nzcount+2)*(nmcount+2),                           &
            mpi_double_precision,ip,                           &
            ias_ybasx,                                         &
            mpi_comm_world,ierr)
    dbuf2=ybasy(izstart-1:izend+1,imstart-1:imend+1)
    call mpi_send(dbuf2,                                    &
            (nzcount+2)*(nmcount+2),                           &
            mpi_double_precision,ip,                           &
            ias_ybasy,                                         &
            mpi_comm_world,ierr)
    deallocate(dbuf2)
    
    allocate(dbuf2(nzcount,nmcount))
    dbuf2=slati(izstart:izend,imstart:imend)
    call mpi_send(dbuf2,                                    &
            nzcount*nmcount,                                   &
            mpi_double_precision,ip,                           &
            ias_slati,                                         &
            mpi_comm_world,ierr)
    dbuf2=clati(izstart:izend,imstart:imend)
    call mpi_send(dbuf2,                                    &
            nzcount*nmcount,                                   &
            mpi_double_precision,ip,                           &
            ias_clati,                                         &
            mpi_comm_world,ierr)
    dbuf2=xlong(izstart:izend,imstart:imend)
    call mpi_send(dbuf2,                                    &
            nzcount*nmcount,                                   &
            mpi_double_precision,ip,                           &
            ias_xlong,                                         &
            mpi_comm_world,ierr)
    dbuf2=xlati(izstart:izend,imstart:imend)
    call mpi_send(dbuf2,                                    &
            nzcount*nmcount,                                   &
            mpi_double_precision,ip,                           &
            ias_xlati,                                         &
            mpi_comm_world,ierr)
    deallocate(dbuf2)
    
    allocate(dbuf4(nzcount,nmcount,nvegtype,nlduse))
    dbuf4=fveg(izstart:izend,imstart:imend,:,:)
    call mpi_send(dbuf4,                                    &
            nzcount*nmcount*nvegtype*nlduse,                   &
            mpi_double_precision,ip,                           &
            ias_fveg,                                          &
            mpi_comm_world,ierr)
    deallocate(dbuf4)
    
    allocate(dbuf3(nzcount,nmcount,nlduse))
    !    print*,' dland'
    dbuf3=dland(izstart:izend,imstart:imend,:)
    call mpi_send(dbuf3,                                    &
            nzcount*nmcount*nlduse,                            &
            mpi_double_precision,ip,                           &
            ias_dland,                                         &
            mpi_comm_world,ierr)
    deallocate(dbuf3)
  
  
  end subroutine send_dbl_arrays
  
  
  
  !*****************************************
  subroutine send_hourly_real_arrays(ip,dom)
    type(dom_type) :: dom
    integer :: ip,izstart,izend,nzcount,imstart,imend,nmcount,ierr
    real(kind=8),dimension(:,:,:),  allocatable :: rbuf3
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf4
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf5
    
    izstart = dom%izstart
    izend   = dom%izend
    imstart = dom%imstart
    imend   = dom%imend
    nzcount = dom%nzcount
    nmcount = dom%nmcount
    
    allocate(rbuf5(nemisa,nzcount,nmcount,nlevemis))
    ! emisaloc
    rbuf5=emisaloc(:, izstart:izend, imstart:imend, :)
    call mpi_send(rbuf5,                                    &
            nemisa*nzcount*nmcount*nlevemis,                   &
            mpi_double_precision,ip,                           &
            ias_emisaloc,                                      &
            mpi_comm_world,ierr)
    deallocate(rbuf5)
    
    allocate(rbuf5(nemisb,nzcount,nmcount,2))
    ! emisb
    rbuf5=emisb(:, izstart:izend, imstart:imend, :)
    call mpi_send(rbuf5,                                    &
            nemisb*nzcount*nmcount*2,                          &
            mpi_double_precision,ip,                           &
            ias_emisb,                                         &
            mpi_comm_world,ierr)
    deallocate(rbuf5)
    
    allocate(rbuf4(nzcount,nmcount,nverti,2))
    ! temp
    rbuf4=temp(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_temp, mpi_comm_world,ierr)
    ! sphu
    rbuf4=sphu(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_sphu, mpi_comm_world,ierr)
    ! airm
    rbuf4=airm(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_airm, mpi_comm_world,ierr)
    ! kzzz
    rbuf4=kzzz(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_kzzz, mpi_comm_world,ierr)
    ! clwc
    rbuf4=clwc(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_clwc, mpi_comm_world,ierr)
    ! lmbb deepconv
    rbuf4=dpeu(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_dpeu, mpi_comm_world,ierr)
    rbuf4=dped(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_dped, mpi_comm_world,ierr)
    rbuf4=dpdu(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_dpdu, mpi_comm_world,ierr)
    rbuf4=dpdd(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_dpdd, mpi_comm_world,ierr)
    
    
    ! winz
    rbuf4=winz(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_winz, mpi_comm_world,ierr)
    ! winm
    rbuf4=winm(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_winm, mpi_comm_world,ierr)
    ! hlay
    rbuf4=hlay(izstart:izend, imstart:imend, :, :)
    call mpi_send(rbuf4,nzcount*nmcount*nverti*2,mpi_double_precision,ip, ias_hlay, mpi_comm_world,ierr)
    deallocate(rbuf4)
    
    allocate(rbuf3(nzcount,nmcount,2))
    ! hght
    rbuf3=hght(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_hght, mpi_comm_world,ierr)
    ! atte
    rbuf3=atte(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_atte, mpi_comm_world,ierr)
    ! tem2
    rbuf3=tem2(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_tem2, mpi_comm_world,ierr)
    ! usta
    rbuf3=usta(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_usta, mpi_comm_world,ierr)
    ! aerr
    rbuf3=aerr(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_aerr, mpi_comm_world,ierr)
    ! obuk
    rbuf3=obuk(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_obuk, mpi_comm_world,ierr)
    ! wsta
    rbuf3=wsta(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_wsta, mpi_comm_world,ierr)
    ! topc
    rbuf3=topc(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_topc, mpi_comm_world,ierr)
    ! sreh
    rbuf3=sreh(izstart:izend, imstart:imend, :)
    call mpi_send(rbuf3,nzcount*nmcount*2,mpi_double_precision,ip, ias_sreh, mpi_comm_world,ierr)
    deallocate(rbuf3)
  
  end subroutine send_hourly_real_arrays
  
  !****************************************
  subroutine send_hourly_real_arrays_tl(ip,dom)
    type(dom_type) :: dom
    integer :: ip,izstart,izend,nzcount,imstart,imend,nmcount,ierr
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf5
    
    izstart = dom%izstart
    izend   = dom%izend
    imstart = dom%imstart
    imend   = dom%imend
    nzcount = dom%nzcount
    nmcount = dom%nmcount
    
    allocate(rbuf5(nemisa,nzcount,nmcount,nlevemis))
    rbuf5=emisaloc_tl(:, izstart:izend, imstart:imend, :)
    call mpi_send(rbuf5,                                    &
            nemisa*nzcount*nmcount*nlevemis,                   &
            mpi_double_precision,ip,                           &
            ias_emisaloc_tl,                                      &
            mpi_comm_world,ierr)
    deallocate(rbuf5)
    
    allocate(rbuf5(nemisb,nzcount,nmcount,2))
    rbuf5=emisb_tl(:, izstart:izend, imstart:imend, :)
    call mpi_send(rbuf5,                                    &
            nemisb*nzcount*nmcount*2,                          &
            mpi_double_precision,ip,                           &
            ias_emisb_tl,                                         &
            mpi_comm_world,ierr)
    deallocate(rbuf5)
  
  end subroutine send_hourly_real_arrays_tl
  
  !*****************************************
  subroutine send_frac_hourly_int_arrays(ip,dom)
    type(dom_type) :: dom
    integer :: ip,izstart,izend,nzcount,imstart,imend,nmcount,ierr
    integer,dimension(:,:),allocatable :: ibuf2
    integer,dimension(:,:,:),allocatable :: ibuf3
    
    izstart = dom%izstart
    izend   = dom%izend
    imstart = dom%imstart
    imend   = dom%imend
    nzcount = dom%nzcount
    nmcount = dom%nmcount
    
    allocate(ibuf2(nzcount,nmcount))
    ibuf2=ideep(izstart:izend,imstart:imend)
    call mpi_send(ibuf2,                          &
            nzcount*nmcount,                         &
            mpi_integer,ip,                          &
            ias_ideep,                               &
            mpi_comm_world,ierr)
    deallocate(ibuf2)
    
    
    allocate(ibuf3(nzcount,nmcount,nverti))
    ! incloud
    ibuf3=incloud(izstart:izend,imstart:imend,:)
    call mpi_send(ibuf3,                          &
            nzcount*nmcount*nverti,                  &
            mpi_integer,ip,                          &
            ias_incloud,                             &
            mpi_comm_world,ierr)
    deallocate(ibuf3)
  
  end subroutine send_frac_hourly_int_arrays
  
  !*****************************************
  subroutine recv_toprint(ip,dom,toprint)
    implicit none
    integer :: ip
    type(dom_type) :: dom
    real(kind=4),dimension(nzonal_domain,nmerid_domain,nverti) :: toprint
    
    real(kind=4),dimension(:,:,:),allocatable :: buf3
    integer :: izstart,izend,nzcount,imstart,imend,nmcount,ierr
    integer,dimension(mpi_status_size) :: status
    
    izstart = dom%izstart
    izend   = dom%izend
    imstart = dom%imstart
    imend   = dom%imend
    nzcount = dom%nzcount
    nmcount = dom%nmcount
    allocate(buf3(nzcount,nmcount,nverti))
    call mpi_recv(buf3,nzcount*nmcount*nverti,mpi_real,ip, iar_toprint, mpi_comm_world,status,ierr)
    toprint(izstart:izend, imstart:imend, :)=buf3
    deallocate(buf3)
  
  end subroutine recv_toprint
  
  
  !*****************************************
  subroutine recv_locvalues(ip,dom)
    implicit none
    type(dom_type) :: dom
    integer :: ip,izstart,izend,nzcount,imstart,imend,nmcount,ierr
    integer,dimension(mpi_status_size) :: status
    real(kind=8),dimension(:,:),allocatable :: dbuf2
    real(kind=8),dimension(:,:,:),allocatable :: dbuf3
    
    izstart = dom%izstart
    izend   = dom%izend
    imstart = dom%imstart
    imend   = dom%imend
    nzcount = dom%nzcount
    nmcount = dom%nmcount
    
    allocate(dbuf3(nzcount,nmcount,nverti))
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_winvloc, mpi_comm_world,status,ierr)
    winvloc(izstart:izend, imstart:imend, :)=dbuf3
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_winxloc, mpi_comm_world,status,ierr)
    winxloc(izstart:izend, imstart:imend, :)=dbuf3
    
    ! sphuloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_sphuloc, mpi_comm_world,status,ierr)
    sphuloc(izstart:izend, imstart:imend, :)=dbuf3
    ! temploc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_temploc, mpi_comm_world,status,ierr)
    temploc(izstart:izend, imstart:imend, :)=dbuf3
    ! clwcloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_clwcloc, mpi_comm_world,status,ierr)
    clwcloc(izstart:izend, imstart:imend, :)=dbuf3
    
    ! airmloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_airmloc, mpi_comm_world,status,ierr)
    airmloc(izstart:izend, imstart:imend, 1:nverti)=dbuf3
    
    ! winzloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_winzloc, mpi_comm_world,status,ierr)
    winzloc(izstart:izend, imstart:imend, :)=dbuf3
    
    ! winmloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_winmloc, mpi_comm_world,status,ierr)
    winmloc(izstart:izend, imstart:imend, :)=dbuf3
    
    ! thlayloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_thlayloc, mpi_comm_world,status,ierr)
    thlayloc(izstart:izend, imstart:imend, :)=dbuf3
    
    ! hlayloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_hlayloc, mpi_comm_world,status,ierr)
    hlayloc(izstart:izend, imstart:imend, :)=dbuf3
    
    ! presloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_presloc, mpi_comm_world,status,ierr)
    presloc(izstart:izend, imstart:imend, :)=dbuf3
    ! kzzzloc
    call mpi_recv(dbuf3,nzcount*nmcount*nverti,mpi_double_precision,ip, iar_kzzzloc, mpi_comm_world,status,ierr)
    kzzzloc(izstart:izend, imstart:imend,:)=dbuf3
    deallocate(dbuf3)
    
    allocate(dbuf3(nspec,nzcount,nmcount))
    ! depoloc
    call mpi_recv(dbuf3,nspec*nzcount*nmcount,mpi_double_precision,ip, iar_depoloc,   mpi_comm_world,status,ierr)
    depoloc(:,izstart:izend, imstart:imend)=dbuf3
    ! drydep
    call mpi_recv(dbuf3,nspec*nzcount*nmcount,mpi_double_precision,ip, iar_drydep,    mpi_comm_world,status,ierr)
    drydep(:,izstart:izend, imstart:imend)=dbuf3
    ! wetdep
    call mpi_recv(dbuf3,nspec*nzcount*nmcount,mpi_double_precision,ip, iar_wetdep,    mpi_comm_world,status,ierr)
    wetdep(:,izstart:izend, imstart:imend) = dbuf3
    deallocate(dbuf3)
    
    allocate(dbuf2(nzcount,nmcount))
    ! hghtloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_hghtloc, mpi_comm_world,status,ierr)
    hghtloc(izstart:izend, imstart:imend)=dbuf2
    ! atteloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_atteloc, mpi_comm_world,status,ierr)
    atteloc(izstart:izend, imstart:imend)=dbuf2
    ! zeniloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_zeniloc, mpi_comm_world,status,ierr)
    zeniloc(izstart:izend, imstart:imend)=dbuf2
    ! tem2loc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_tem2loc, mpi_comm_world,status,ierr)
    tem2loc(izstart:izend, imstart:imend)=dbuf2
    ! ustaloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_ustaloc, mpi_comm_world,status,ierr)
    ustaloc(izstart:izend, imstart:imend)=dbuf2
    ! aerrloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_aerrloc, mpi_comm_world,status,ierr)
    aerrloc(izstart:izend, imstart:imend)=dbuf2
    ! obukloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_obukloc, mpi_comm_world,status,ierr)
    obukloc(izstart:izend, imstart:imend)=dbuf2
    ! wstaloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_wstaloc, mpi_comm_world,status,ierr)
    wstaloc(izstart:izend, imstart:imend)=dbuf2
    ! topcloc
    call mpi_recv(dbuf2,nzcount*nmcount,mpi_double_precision,ip, iar_topcloc, mpi_comm_world,status,ierr)
    topcloc(izstart:izend, imstart:imend)=dbuf2
    deallocate(dbuf2)
  
  
  end subroutine recv_locvalues
  
  !**********************
  subroutine recv_locvalues_tl(ip,dom)
    implicit none
    type(dom_type) :: dom
    integer :: ip,izstart,izend,nzcount,imstart,imend,nmcount,ierr
    integer,dimension(mpi_status_size) :: status
    real(kind=8),dimension(:,:,:),allocatable :: dbuf3
    
    izstart = dom%izstart
    izend   = dom%izend
    imstart = dom%imstart
    imend   = dom%imend
    nzcount = dom%nzcount
    nmcount = dom%nmcount
    
    allocate(dbuf3(nspec,nzcount,nmcount))
    call mpi_recv(dbuf3,nspec*nzcount*nmcount,mpi_double_precision,ip, iar_drydep_tl,    mpi_comm_world,status,ierr)
    drydep_tl(:,izstart:izend, imstart:imend)=dbuf3
    ! wetdep
    call mpi_recv(dbuf3,nspec*nzcount*nmcount,mpi_double_precision,ip, iar_wetdep_tl,    mpi_comm_world,status,ierr)
    wetdep_tl(:,izstart:izend, imstart:imend) = dbuf3
    deallocate(dbuf3)
  end subroutine recv_locvalues_tl
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine amaster_send_aconc
    integer :: ip,ierr
    real(kind=8),dimension(:,:,:,:),allocatable :: dbuf4
    integer:: izstart,izend,imstart,imend,i,j
    
    do ip=1,ndoms
      allocate(dbuf4(nspectot, dom(ip)%nzcount+6, dom(ip)%nmcount+6, nverti+1))
      ! aconc
      i=dom(ip)%i
      j=dom(ip)%j
      izstart=dom(ip)%izstart
      izend=dom(ip)%izend
      imstart=dom(ip)%imstart
      imend=dom(ip)%imend
      dbuf4=aconcsave(:, izstart+(i-1)*6:izend+(i*6),imstart+(j-1)*6:imend+(j*6), :)
      call mpi_send( &
              dbuf4,                                                       &
              nspectot*(dom(ip)%nzcount+6)*(dom(ip)%nmcount+6)*(nverti+1), &
              mpi_double_precision,ip,                                     &
              ias_aconc,                                                    &
              mpi_comm_world,ierr)
      deallocate(dbuf4)
      
      
      allocate(dbuf4(nspectot, dom(ip)%nzcount, dom(ip)%nmcount, nverti))
      ! aconco
      dbuf4=aconco(:, dom(ip)%izstart:dom(ip)%izend, dom(ip)%imstart:dom(ip)%imend, :)
      call mpi_send( &
              dbuf4,                                                       &
              nspectot*(dom(ip)%nzcount)*(dom(ip)%nmcount)*nverti, &
              mpi_double_precision,ip,                                     &
              ias_aconco,                                                    &
              mpi_comm_world,ierr)
      deallocate(dbuf4)
    
    end do
    
    aconc(:,:,:,:)=0.d0
    aconco(:,:,:,:)=0.d0
  
  end subroutine amaster_send_aconc
  !******************************************
  subroutine amaster_send_once
    
    integer :: ip,nzcount,nmcount,izstart,imstart,izend,imend,ierr
    real(kind=8),dimension(:,:,:,:,:),allocatable :: rbuf5
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf4
    
    do ip=1,ndoms
      
      nzcount = dom(ip)%nzcount
      nmcount = dom(ip)%nmcount
      izstart = dom(ip)%izstart
      izend   = dom(ip)%izend
      imstart = dom(ip)%imstart
      imend   = dom(ip)%imend
      
      allocate(rbuf5(0:nhourrun+1,nemisa,nzcount,nmcount,nlevemis))
      rbuf5=aemisaloc(:,:, izstart:izend, imstart:imend, :)
      call mpi_send(rbuf5,                                    &
              (nhourrun+2)*nemisa*nzcount*nmcount*nlevemis,      &
              mpi_double_precision,ip,                           &
              ias_aemisaloc,                                      &
              mpi_comm_world,ierr)
      
      deallocate(rbuf5)
      
      allocate(rbuf4(0:nhourrun+1,nemisb,nzcount,nmcount))
      rbuf4=aemisb(:,:, izstart:izend, imstart:imend)
      call mpi_send(rbuf4,                                    &
              (nhourrun+2)*nemisb*nzcount*nmcount,             &
              mpi_double_precision,ip,                           &
              ias_aemisb,                                         &
              mpi_comm_world,ierr)
      deallocate(rbuf4)
    
    enddo
  
  end subroutine amaster_send_once
  !******************************************
  subroutine amaster_recv_once
    
    integer :: ip,nzcount,nmcount,izstart,imstart,izend,imend,ierr
    real(kind=8),dimension(:,:,:,:,:),allocatable :: rbuf5
    real(kind=8),dimension(:,:,:,:),allocatable :: rbuf4
    integer,dimension(mpi_status_size) :: status
    
    do ip=1,ndoms
      
      nzcount = dom(ip)%nzcount
      nmcount = dom(ip)%nmcount
      izstart = dom(ip)%izstart
      izend   = dom(ip)%izend
      imstart = dom(ip)%imstart
      imend   = dom(ip)%imend
      
      allocate(rbuf5(nhourrun+2,nemisa,nzcount,nmcount,nlevemis))
      call mpi_recv(rbuf5,                                    &
              (nhourrun+2)*nemisa*nzcount*nmcount*nlevemis,      &
              mpi_double_precision,ip,                           &
              iar_aemisaloc,                                      &
              mpi_comm_world,status,ierr)
      aemisaloc(0:nhourrun+1,:,izstart:izend,imstart:imend,:)  =&
              rbuf5(1:nhourrun+2,:,1:nzcount,1:nmcount,:)
      deallocate(rbuf5)
      
      allocate(rbuf4(nhourrun+2,nemisb,nzcount,nmcount))
      call mpi_recv(rbuf4,                                    &
              (nhourrun+2)*nemisb*nzcount*nmcount,             &
              mpi_double_precision,ip,                           &
              iar_aemisb,                                         &
              mpi_comm_world,status,ierr)
      aemisb(0:nhourrun+1,:, izstart:izend, imstart:imend) = &
              rbuf4(1:nhourrun+2,:,1:nzcount,1:nmcount)
      deallocate(rbuf4)
    
    enddo
  
  end subroutine amaster_recv_once
  
  !******************************
  subroutine amaster_recv_aconc
    integer :: ip,ime,izo,ns,ivert
    integer :: ierr
    integer,dimension(mpi_status_size) :: status
    real(kind=8),allocatable,dimension(:,:,:,:) :: dbuf4
    integer:: izstart,izend,imstart,imend,i,j
    
    do ip=1,ndoms
      allocate(dbuf4(nspectot,dom(ip)%nzcount+6,dom(ip)%nmcount+6,nverti+1))
      ! aconc
      i=dom(ip)%i
      j=dom(ip)%j
      izstart=dom(ip)%izstart
      izend=dom(ip)%izend
      imstart=dom(ip)%imstart
      imend=dom(ip)%imend
      call mpi_recv(                                            &
              dbuf4,                                               &
              nspectot*(dom(ip)%nzcount+6)*(dom(ip)%nmcount+6)*(nverti+1), &
              mpi_double_precision,ip,                             &
              iar_aconc,                                            &
              mpi_comm_world,status,ierr                               &
              )
      
      do ns=1,nspectot
        do ivert=1,nverti+1
          do izo=1,dom(ip)%nzcount+6
            do ime=1,dom(ip)%nmcount+6
              aconcsave(ns,izstart+(i-1)*6+izo-1,imstart+(j-1)*6+ime-1,ivert)=dbuf4(ns,izo,ime,ivert)
            enddo
          enddo
        enddo
      enddo
      !       aconcsave(:,izstart+(i-1)*6:izend+(i*6),imstart+(j-1)*6:imend+(j*6),:) =  dbuf4
      deallocate(dbuf4)
      allocate(dbuf4(nspectot,dom(ip)%nzcount,dom(ip)%nmcount,nverti+1))
      ! aconco
      call mpi_recv(                                            &
              dbuf4,                                               &
              nspectot*dom(ip)%nzcount*dom(ip)%nmcount*nverti, &
              mpi_double_precision,ip,                             &
              iar_aconco,                                            &
              mpi_comm_world,status,ierr                               &
              )
      aconco(:,dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,:) = &
              aconco(:,dom(ip)%izstart:dom(ip)%izend,dom(ip)%imstart:dom(ip)%imend,:) + dbuf4
      deallocate(dbuf4)
    end do
  
  end subroutine amaster_recv_aconc
  
  !*********************************
  subroutine amaster_recv_aconc_bounds
    integer :: ip,ierr,izstart,izend,imstart,imend,nzcount,nmcount,idom,jdom
    integer,dimension(mpi_status_size) :: status
    real(kind=8),allocatable,dimension(:,:,:) :: topbuf
    real(kind=8),allocatable,dimension(:,:,:,:) :: hlbuf,hubuf,vlbuf,vrbuf
    
    do ip=1,ndoms
      izstart = dom(ip)%izstart
      izend   = dom(ip)%izend
      imstart = dom(ip)%imstart
      imend   = dom(ip)%imend
      nzcount = dom(ip)%nzcount
      nmcount = dom(ip)%nmcount
      idom    = dom(ip)%i
      jdom    = dom(ip)%j
      
      allocate(topbuf(nspectot,nzcount,nmcount))
      call mpi_recv(                                            &
              topbuf,                                               &
              nspectot*nzcount*nmcount, &
              mpi_double_precision,ip,                             &
              iar_aconc_top,                                            &
              mpi_comm_world,status,ierr                               &
              )
      !   if(ip==1)write(*,1002)aconc(1,9,1,nverti+1),topbuf(1,9,1)
      aconc(:,izstart:izend,imstart:imend,nverti+1)= &
              aconc(:,izstart:izend,imstart:imend,nverti+1) + topbuf
      deallocate(topbuf)
      !  if(ip==1)write(*,1001)aconc(1,9,1,nverti+1)
      !1001 format('RECV ',e64.56)
      !1002 format('AV ',2(e64.56,' '))
      
      if(idom==1) then ! left boundary
        allocate(vlbuf(nspectot, 3, nmcount+6, nverti+1))
        call mpi_recv(                                            &
                vlbuf,                                               &
                nspectot*3*(nmcount+6)*(nverti+1), &
                mpi_double_precision,ip,                             &
                iar_aconc_bounds,                                            &
                mpi_comm_world,status,ierr                               &
                )
        aconc(:,-2:0,imstart-3:imend+3,:) = &
                aconc(:,-2:0,imstart-3:imend+3,:) + vlbuf(:,:,:,:)
        deallocate(vlbuf)
      endif
      
      if(idom==nzdoms) then ! right boundary
        allocate(vrbuf(nspectot, 3, nmcount+6, nverti+1))
        call mpi_recv(                                            &
                vrbuf,                                               &
                nspectot*3*(nmcount+6)*(nverti+1), &
                mpi_double_precision,ip,                             &
                iar_aconc_bounds,                                            &
                mpi_comm_world,status,ierr                               &
                )
        aconc(:,nzonal_domain+1:nzonal_domain+3,imstart-3:imend+3,:) = &
                aconc(:,nzonal_domain+1:nzonal_domain+3,imstart-3:imend+3,:) + vrbuf(:,:,:,:)
        deallocate(vrbuf)
      endif
      
      if(jdom==1) then ! lower boundary
        allocate(hlbuf(nspectot, nzcount+6, 3, nverti+1))
        call mpi_recv(                                            &
                hlbuf,                                               &
                nspectot*3*(nzcount+6)*(nverti+1), &
                mpi_double_precision,ip,                             &
                iar_aconc_bounds,                                            &
                mpi_comm_world,status,ierr                               &
                )
        aconc(:,izstart-3:izend+3,-2:0,:) = &
                aconc(:,izstart-3:izend+3,-2:0,:) + hlbuf(:,:,:,:)
        deallocate(hlbuf)
      endif
      
      if(jdom==nmdoms) then ! upper boundary
        allocate(hubuf(nspectot, nzcount+6, 3, nverti+1))
        call mpi_recv(                                            &
                hubuf,                                               &
                nspectot*3*(nzcount+6)*(nverti+1), &
                mpi_double_precision,ip,                             &
                iar_aconc_bounds,                                            &
                mpi_comm_world,status,ierr                               &
                )
        aconc(:,izstart-3:izend+3,nmerid_domain+1:nmerid_domain+3,:) = &
                aconc(:,izstart-3:izend+3,nmerid_domain+1:nmerid_domain+3,:) + hubuf(:,:,:,:)
        deallocate(hubuf)
      endif
    
    
    enddo
  
  end subroutine amaster_recv_aconc_bounds

end module master_message_subs
