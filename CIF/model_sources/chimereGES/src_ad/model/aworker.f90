subroutine aworker
  
  use chimere_consts
  use worker_common
  use message_defs
  use worker_message_subs
  use twostep_mod
  use atwostep_mod
  
  implicit none
  
  integer :: ir, np, nh, icpt
  integer :: ierr
  ! IP
  integer i, ns
  integer ime, izo, ivert
  
  CHARACTER*2 :: idstr
  character*4 :: idstr2
  character*10 :: str
  REAL(kind = 8) :: p_mid, ddp, psurf0, MMair
  real(kind = 8) :: fppb, fdry
  real(kind = 8), allocatable, dimension(:) :: stockthour ! for precision's sake
  ! fin IP
  
  call aworker_recv_once !initial aemis
  
  call mpi_comm_rank(wrk_comm, rank, ierr)
  
  ihourrun = ihourrun - 1 !pour bien commencer l'adjoint
  !  print*,'IHOURRUN ADJ worker',ihourrun
  
  call mpi_barrier(mpi_comm_world, ierr)
  
  ! recv aconc et aconco initiales
  call aworker_recv_aconc
  
  call mpi_barrier(mpi_comm_world, ierr)
  
  ! Loop on physical time steps
  do nh = nhourrun, 1, -1
    
    !    print*,'RANK=',rank,nh,nhourrun
    call mpi_barrier(mpi_comm_world, ierr)
    
    call worker_recv_hourly
    thour = 0.5d0 / dble(nphour)
    do np = 1, nphour
      djul = djul - dun / (nphour * nhourpday)
    enddo
    
    call mpi_barrier(mpi_comm_world, ierr)
    
    call awrite_concs_hour(3, ihourrun) ! c!atherine relecture conc IP: hour
    
    call awrite_concs_phys(1, rank) ! ouverture fichier ecriture conc IP phys car
    ! pas teste avec ichemstep != 1
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print*,' DIRECT INTEGRATION',rank,ihourrun
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    icpt = 0
    
    allocate(stockthour(nphour)) ! needed for precision
    
    do np = 1, nphour
      
      call mpi_barrier(mpi_comm_world, ierr)
      stockthour(np) = thour
      do ns = 1, nspec
        call worker_recv_conc_bounds(ns)
      end do
      
      call mpi_barrier(mpi_comm_world, ierr) !modif
      
      call worker_recv_frac_hourly
      call zenith
      call locvalues
      call physics2
      
      call mpi_barrier(mpi_comm_world, ierr)
      do ir = 1, ichemstep
        icpt = icpt + 1
        call awrite_concs_phys(2, icpt)
        call twostep
        call worker_update_halo ! between workers. Master is not involved.
      enddo
      
      ! IP
      ! RQ: pour le moment, on ne descend pas au niveau
      ! du pas de temps chimique pour la comparaison des obs
      ! RQ2: si on a des obs de familles, mettre ici a jour le calcul
      ! des familles (voir outprint)??
      
      if(any(tabobs(:, 1)==ihourrun)) then
        do i = 1, nobs
          if(tabobs(i, 1)==ihourrun .and.tabobs(i, 2)==np .and. tabobs(i, 6) >= 1) then
            ns = tabobs(i, 6)
            ivert = tabobs(i, 5)
            ime = tabobs(i, 3) - imstart + 1
            izo = tabobs(i, 4) - izstart + 1
            
            !tabobs(i,3)=nlat in the whole domain
            if (species(ns)%fspec > dzero) then
              ! conversion molec.cm-3 to the unit specified by the user
              tabobs(i, 7) = species(ns)%fspec * conc(ns, izo, ime, ivert)
            elseif(species(ns)%fspec == dzero) then
              fppb = 1d9 / airmloc(izo, ime, ivert)
              fdry = m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6) &
                      / (m_air * sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6 &
                              + m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6))
              tabobs(i, 7) = conc(ns, izo, ime, ivert) * fppb / fdry
            else
              ! conversion molec.cm-3 to ppb
              fppb = 1d9 / airmloc(izo, ime, ivert)
              tabobs(i, 7) = conc(ns, izo, ime, ivert) * fppb
            endif
            
            if (ivert.eq.1) then
              psurf0 = psurf
              if (psurf .le. presloc(izo, ime, ivert)) then
                MMair = 0.028976 !in kg.mol
                psurf0 = presloc(izo, ime, ivert) + g * airmloc(izo, ime, ivert) * hlayloc(izo, ime, ivert) / (avogadro * MMair)
              endif
              ddp = psurf0 - presloc(izo, ime, ivert)
              p_mid = psurf0 - ddp / 2.
            
            else
              ddp = presloc(izo, ime, ivert - 1) - presloc(izo, ime, ivert)
              p_mid = presloc(izo, ime, ivert - 1) - ddp / 2
            endif
            
            tabobs(i, 9) = p_mid                       !in Pa
            tabobs(i, 10) = ddp                        !in Pa
            tabobs(i, 11) = airmloc(izo, ime, ivert)   !in molec.cm-3
            tabobs(i, 12) = thlayloc(izo, ime, ivert)  !in cm
          endif! ihourrun,np,species
        enddo ! nobs
      endif! ihourrun
      ! end IP
      
      !  Timing update
      thour = thour + dun / dble(nphour)
      djul = djul + dun / (nphour * nhourpday)
    
    enddo ! np=1,nphour
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !print*,' Adjoint integration ',rank
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'nh=', nh !,' aconc=',aconc(1,:,nmerid,nverti)
    
    do np = nphour, 1, -1
      
      !  Timing update
      ! thour = thour - dun/dble(nphour)  ! catherine
      ! precision problem - > new formulation for thour
      !thour = 0.5d0/dble(nphour) + dble(np-1)/dble(nphour)
      thour = stockthour(np)
      djul = djul - dun / (nphour * nhourpday) !catherine IP a verifier
      call mpi_barrier(mpi_comm_world, ierr)
      do ns = 1, nspec
        call worker_recv_conc_bounds(ns)
      end do
      call mpi_barrier(mpi_comm_world, ierr)
      
      call worker_recv_frac_hourly
      call zenith
      call locvalues
      call physics2
      
      call mpi_barrier(mpi_comm_world, ierr)
      ! IP initalisation aconc
      ! RQ: pour le moment, on ne descend pas au niveau
      ! du pas de temps chimique pour la comparaison des obs  
      if(any(tabobs(:, 1)==ihourrun)) then
        do i = 1, nobs
          if(tabobs(i, 1)==ihourrun .and.tabobs(i, 2)==np .and. tabobs(i, 6) >= 1) then
            ns = tabobs(i, 6)
            ivert = tabobs(i, 5)
            ime = tabobs(i, 3) - imstart + 1
            izo = tabobs(i, 4) - izstart + 1

            if (species(ns)%fspec > dzero) then
              ! conversion unit specified by the user / molec.cm-3
              aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + &
                      species(ns)%fspec * tabobs(i, 13)
            elseif(species(ns)%fspec == dzero) then
              fppb = 1d9 / airmloc(izo, ime, ivert)
              fdry = m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6) &
                      / (m_air * sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6 &
                              + m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6))
              aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + fppb / fdry * tabobs(i, 13)
            else
              ! conversion ppb / molec.cm-3
              fppb = 1d9 / airmloc(izo, ime, ivert)
              aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) + fppb * tabobs(i, 13)
            endif
       
            if (ivert.eq.1) then
              psurf0 = psurf
              if (psurf .le. presloc(izo, ime, ivert)) then
                MMair = 0.028976 !in kg.mol
                psurf0 = presloc(izo, ime, ivert) + g * airmloc(izo, ime, ivert) * hlayloc(izo, ime, ivert) / (avogadro * MMair)
              endif
              ddp = psurf0 - presloc(izo, ime, ivert)
              p_mid = psurf0 - ddp / 2.
            
            else
              ddp = presloc(izo, ime, ivert - 1) - presloc(izo, ime, ivert)
              p_mid = presloc(izo, ime, ivert - 1) - ddp / 2
            endif
            
            tabobs(i, 9) = p_mid                       !in Pa
            tabobs(i, 10) = ddp                        !in Pa
            tabobs(i, 11) = airmloc(izo, ime, ivert)   !in molec.cm-3
            tabobs(i, 12) = thlayloc(izo, ime, ivert)  !in cm
          endif ! ihourrun,np
        enddo ! nobs
      endif
      !end IP
      aemisbloc = dzero ! iniadjloc
      ! ici dans adjinteg, aoutprint = aconc pour les familles + outprint
      ! mais necessite passages avec master donc pas fait pour le moment
      do ir = ichemstep, 1, -1
        call awrite_concs_phys(3, icpt)
        icpt = icpt - 1
        !! COMMENTE parce qu'on ne sait pas quoi en faire exactement
        !!!!!!!!!          call aworker_update_halo
        call atwostep
        !!!!!!!!!   	   call aworker_update_halo
      enddo
      
      call mpi_barrier(mpi_comm_world, ierr)
      call aphysics2 ! pour aconc inactive species
      call alocvalues ! valable pour aemisb
      call aworker_send_aconc_bounds
      call mpi_barrier(mpi_comm_world, ierr)
 
    end do ! np=nphour
    
    call awrite_concs_phys(4, 0)
    deallocate(stockthour)
    if(nh.gt.1)ihourrun = ihourrun - 1 ! sinon, on n'avance pas!
  
  end do ! nh=nhourrun,1,-1
  
  call aphysics2 ! pour aconc inactive species
  call mpi_barrier(mpi_comm_world, ierr)
  ! Save of all aconcentrations
  
  call aworker_send_aconc
  
  call mpi_barrier(mpi_comm_world, ierr)
  !IP envoi des tabobs au maitre
  call mpi_send(tabobs, nobs * 13, mpi_double_precision, 0, &
          iar_tabobs, mpi_comm_world, ierr)
  ! fin IP
  
  ! envoi des aemis au maitre aemisaloc(0:nhourrun+1,1:nemisa,1:nzonal,1:nmerid,1:nlevemis)
  call aworker_send_once

end subroutine aworker
