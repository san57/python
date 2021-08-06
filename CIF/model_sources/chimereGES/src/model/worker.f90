subroutine worker
  
  use chimere_consts
  use worker_common
  use message_defs
  use worker_message_subs
  use twostep_mod
  
  implicit none
  
  integer :: ir, np, nh
  integer :: ierr
  
  ! PYVAR
  integer i, ns
  integer ime, izo, ivert
  CHARACTER*2 :: idstr
  character*4 :: idstr2
  character*10 :: str
  REAL(kind = 8) :: p_mid, ddp, psurf0, MMair
  real(kind = 8) :: fppb, fdry
  ! end PYVAR
  
  !print*,rank,' initial time-step values',nphour_ref,ichemstep
  dtr = 3600d0 / (nphour_ref * ichemstep)
  dtr2 = 2d0 * dtr / 3d0
  nphour = nphour_ref
  
  ihourrun = 0
  thour = 0.5d0 / dble(nphour)
  
  call mpi_barrier(mpi_comm_world, ierr)
  
  call zenith
  call locvalues
  call physics2
  
  !print*,' send locvalues to master for printing of initial state'
  call worker_send_locvalues
  call mpi_barrier(mpi_comm_world, ierr)
  
  ! outprint
  call prep_outprint(1)
  call mpi_barrier(mpi_comm_world, ierr)
  if(ad.eq.1) call awrite_concs_hour(1, rank)
  
  !print*,'Loop on physical time steps'
  do nh = 1, nhourrun
    call mpi_barrier(mpi_comm_world, ierr)
    call worker_recv_hourly
    thour = 0.5d0 / dble(nphour)
    call mpi_barrier(mpi_comm_world, ierr)
    if(ad.eq.1) call awrite_concs_hour(2, ihourrun)
    ! pulse hour Elise Potier
    if (nh.ge.hpulse) then
      do np = 1, nphour
        call mpi_barrier(mpi_comm_world, ierr)
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
            endif ! ihourrun,np,species
          enddo  ! nobs
        
        endif ! ihourrun
        ! end IP
        
        !  Timing update
        thour = thour + dun / dble(nphour)
        djul = djul + dun / (nphour * nhourpday)
      
      enddo ! np=1,nphour
    endif ! hpulse
    ihourrun = ihourrun + 1
    ! Printout. Master needs locvalues and conc
    call worker_send_locvalues
    call prep_outprint(nh) ! also sends conc to master
  enddo ! nh=1,nhourrun
  
  !print*,'IP sending tabobs to the master'
  if (fwd.eq.1) call mpi_send(tabobs, nobs * 12, mpi_double_precision, 0, &
          iar_tabobs, mpi_comm_world, ierr)
  if(fwd.eq.0)call mpi_send(tabobs, nobs * 13, mpi_double_precision, 0, &
          iar_tabobs, mpi_comm_world, ierr)
  
  ! fin IP

end subroutine worker
