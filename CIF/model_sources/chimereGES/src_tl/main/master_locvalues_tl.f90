subroutine master_locvalues_tl
  
  use chimere_consts
  use chimere_common
  use master_message_subs
  use wholedomain_common
  implicit none
  
  ! lateral boundary conditions are the concentration values for
  ! the points surrrounding the nzonal_domain*nmerid_domain rectangle
  !*****************************************************************************************
  integer :: ns, ilat, ivert, izo, ime
  real(kind = 8) :: w1, w2, fdry
  !*****************************************************************************************
  w1 = thour
  w2 = dun - w1
  
  ! I need airmloc for conc calculation. I could get it from workers.
  ! But it is faster to recalculate it here, since I only need boundary values
  
  do ivert = 1, nverti
    do izo = 1, nzonal_domain, nzonal_domain - 1
      do ime = 1, nmerid_domain
        airmloc(izo, ime, ivert) = w2 * airm(izo, ime, ivert, 1) + w1 * airm(izo, ime, ivert, 2)
        sphuloc(izo, ime, ivert) = w2 * sphu(izo, ime, ivert, 1) + w1 * sphu(izo, ime, ivert, 2)
      end do
    end do
    do ime = 1, nmerid_domain, nmerid_domain - 1
      do izo = 1, nzonal_domain
        airmloc(izo, ime, ivert) = w2 * airm(izo, ime, ivert, 1) + w1 * airm(izo, ime, ivert, 2)
        sphuloc(izo, ime, ivert) = w2 * sphu(izo, ime, ivert, 1) + w1 * sphu(izo, ime, ivert, 2)
      end do
    end do
  end do
  
  do ime = 1, nmerid_domain
    do izo = 1, nzonal_domain
      airmloc(izo, ime, nverti + 1) = w2 * airm(izo, ime, nverti, 1) + w1 * airm(izo, ime, nverti, 2)
      sphuloc(izo, ime, nverti + 1) = w2 * sphu(izo, ime, nverti, 1) + w1 * sphu(izo, ime, nverti, 2)
    end do
  end do
  
  airmloc(0, :, :) = airmloc(1, :, :)
  airmloc(nzonal_domain + 1, :, :) = airmloc(nzonal_domain, :, :)
  airmloc(:, 0, :) = airmloc(:, 1, :)
  airmloc(:, nmerid_domain + 1, :) = airmloc(:, nmerid_domain, :)
  
  sphuloc(0, :, :) = sphuloc(1, :, :)
  sphuloc(nzonal_domain + 1, :, :) = sphuloc(nzonal_domain, :, :)
  sphuloc(:, 0, :) = sphuloc(:, 1, :)
  sphuloc(:, nmerid_domain + 1, :) = sphuloc(:, nmerid_domain, :)
  
  do ns = 1, nspec
    ! lat
    ilat = 0
    do ivert = 1, nverti
      do izo = 0, nzonal_domain + 1, nzonal_domain + 1
        do ime = 1, nmerid_domain
          ilat = ilat + 1
          conc_tl(izo, ime, ivert) = w2 * boundlat_tl(ilat, ns, 1) + w1 * boundlat_tl(ilat, ns, 2)
          conc(izo, ime, ivert) = w2 * boundlat(ilat, ns, 1) + w1 * boundlat(ilat, ns, 2)
          conc_tl(izo, ime, ivert) = conc_tl(izo, ime, ivert) * airmloc(izo, ime, ivert) * 1d-9
          conc(izo, ime, ivert) = conc(izo, ime, ivert) * airmloc(izo, ime, ivert) * 1d-9
          
          if(species(ns)%bounddry.eq.1) then
            fdry = m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6) &
                / (m_air * sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6 &
                    + m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6))
            conc(izo, ime, ivert) = conc(izo, ime, ivert) * fdry
            conc_tl(izo, ime, ivert) = conc_tl(izo, ime, ivert) * fdry
          end if
        end do
      end do
      do ime = 0, nmerid_domain + 1, nmerid_domain + 1
        do izo = 1, nzonal_domain
          ilat = ilat + 1
          conc_tl(izo, ime, ivert) = w2 * boundlat_tl(ilat, ns, 1) + w1 * boundlat_tl(ilat, ns, 2)
          conc(izo, ime, ivert) = w2 * boundlat(ilat, ns, 1) + w1 * boundlat(ilat, ns, 2)
          conc_tl(izo, ime, ivert) = conc_tl(izo, ime, ivert) * airmloc(izo, ime, ivert) * 1d-9
          conc(izo, ime, ivert) = conc(izo, ime, ivert) * airmloc(izo, ime, ivert) * 1d-9
          if(species(ns)%bounddry.eq.1) then
            fdry = m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6) &
                / (m_air * sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6 &
                    + m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6))
            conc(izo, ime, ivert) = conc(izo, ime, ivert) * fdry
            conc_tl(izo, ime, ivert) = conc_tl(izo, ime, ivert) * fdry
          end if
        end do
      end do
    end do
    
    ! top
    do ime = 1, nmerid_domain
      do izo = 1, nzonal_domain
        conc_tl(izo, ime, nverti + 1) = w2 * boundtop_tl(izo, ime, ns, 1) + w1 * boundtop_tl(izo, ime, ns, 2)
        conc(izo, ime, nverti + 1) = w2 * boundtop(izo, ime, ns, 1) + w1 * boundtop(izo, ime, ns, 2)
        conc_tl(izo, ime, nverti + 1) = conc_tl(izo, ime, nverti + 1) * airmloc(izo, ime, nverti + 1) * 1d-9
        conc(izo, ime, nverti + 1) = conc(izo, ime, nverti + 1) * airmloc(izo, ime, nverti + 1) * 1d-9
        if(species(ns)%bounddry.eq.1) then
          fdry = m_h2o * (1 - sphuloc(izo, ime, nverti + 1) / airmloc(izo, ime, nverti + 1) / 1.6) &
              / (m_air * sphuloc(izo, ime, nverti + 1) / airmloc(izo, ime, nverti + 1) / 1.6 &
                  + m_h2o * (1 - sphuloc(izo, ime, nverti + 1) / airmloc(izo, ime, nverti + 1) / 1.6))
          conc(izo, ime, nverti + 1) = conc(izo, ime, nverti + 1) * fdry
          conc_tl(izo, ime, nverti + 1) = conc_tl(izo, ime, nverti + 1) * fdry
        end if
      end do
    end do
    
    conc_tl(-1, :, :) = conc_tl(0, :, :)
    conc_tl(-2, :, :) = conc_tl(0, :, :)
    conc_tl(nzonal_domain + 2, :, :) = conc_tl(nzonal_domain + 1, :, :)
    conc_tl(nzonal_domain + 3, :, :) = conc_tl(nzonal_domain + 1, :, :)
    conc(-1, :, :) = conc(0, :, :)
    conc(-2, :, :) = conc(0, :, :)
    conc(nzonal_domain + 2, :, :) = conc(nzonal_domain + 1, :, :)
    conc(nzonal_domain + 3, :, :) = conc(nzonal_domain + 1, :, :)
    conc_tl(:, -1, :) = conc_tl(:, 0, :)
    conc_tl(:, -2, :) = conc_tl(:, 0, :)
    conc_tl(:, nmerid_domain + 2, :) = conc_tl(:, nmerid_domain + 1, :)
    conc_tl(:, nmerid_domain + 3, :) = conc_tl(:, nmerid_domain + 1, :)
    conc(:, -1, :) = conc(:, 0, :)
    conc(:, -2, :) = conc(:, 0, :)
    conc(:, nmerid_domain + 2, :) = conc(:, nmerid_domain + 1, :)
    conc(:, nmerid_domain + 3, :) = conc(:, nmerid_domain + 1, :)
    
    ! send boundaries for the current species
    call master_send_conc_bounds_tl
    call master_send_conc_bounds
  
  end do ! do ns=1,nspec

end subroutine master_locvalues_tl
