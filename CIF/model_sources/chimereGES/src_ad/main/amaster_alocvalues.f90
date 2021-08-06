subroutine amaster_alocvalues
  
  use chimere_consts
  use chimere_common
  use wholedomain_common
  
  implicit none
  
  ! lateral boundary conditions are the concentration values for
  ! the points surrrounding the nzonal_domain*nmerid_domain rectangle
  !*****************************************************************************************
  integer :: ns, ilat, ivert, izo, ime
  real(kind = 8) :: w1, w2, fdry
  !real(kind=8),dimension(0:nzonal_domain+1,0:nmerid_domain+1,1:nverti+1) :: airmtmp
  !*****************************************************************************************
  w1 = thour
  w2 = dun - w1
  
  ! airmloc already recomputed by master_locvalues
  aconc(:, :, nmerid_domain + 1, :) = aconc(:, :, nmerid_domain + 1, :) + aconc(:, :, nmerid_domain + 3, :)
  aconc(:, :, nmerid_domain + 3, :) = dzero
  aconc(:, :, nmerid_domain + 1, :) = aconc(:, :, nmerid_domain + 1, :) + aconc(:, :, nmerid_domain + 2, :)
  aconc(:, :, nmerid_domain + 2, :) = dzero
  aconc(:, :, 0, :) = aconc(:, :, 0, :) + aconc(:, :, -2, :)
  aconc(:, :, -2, :) = dzero
  aconc(:, :, 0, :) = aconc(:, :, 0, :) + aconc(:, :, -1, :)
  aconc(:, :, -1, :) = dzero
  
  aconc(:, nzonal_domain + 1, :, :) = aconc(:, nzonal_domain + 1, :, :) + aconc(:, nzonal_domain + 3, :, :)
  aconc(:, nzonal_domain + 3, :, :) = dzero
  aconc(:, nzonal_domain + 1, :, :) = aconc(:, nzonal_domain + 1, :, :) + aconc(:, nzonal_domain + 2, :, :)
  aconc(:, nzonal_domain + 2, :, :) = dzero
  aconc(:, 0, :, :) = aconc(:, 0, :, :) + aconc(:, -2, :, :)
  aconc(:, -2, :, :) = dzero
  aconc(:, 0, :, :) = aconc(:, 0, :, :) + aconc(:, -1, :, :)
  aconc(:, -1, :, :) = dzero
  
  
  ! top
  do ns = 1, nspec
    do ime = 1, nmerid_domain
      do izo = 1, nzonal_domain
        if(species(ns)%bounddry.eq.1) then
          fdry = m_h2o * (1 - sphuloc(izo, ime, nverti + 1) / airmloc(izo, ime, nverti + 1) / 1.6) &
              / (m_air * sphuloc(izo, ime, nverti + 1) / airmloc(izo, ime, nverti + 1) / 1.6 &
                  + m_h2o * (1 - sphuloc(izo, ime, nverti + 1) / airmloc(izo, ime, nverti + 1) / 1.6))
          aconc(ns, izo, ime, nverti + 1) = aconc(ns, izo, ime, nverti + 1) * fdry
        end if
        
        aconc(ns, izo, ime, nverti + 1) = aconc(ns, izo, ime, nverti + 1) &
            * airmloc(izo, ime, nverti + 1) * 1d-9
        aboundtop(izo, ime, ns, ihourrun) = aboundtop(izo, ime, ns, ihourrun) &
            + w2 * aconc(ns, izo, ime, nverti + 1)
        aboundtop(izo, ime, ns, ihourrun + 1) = aboundtop(izo, ime, ns, ihourrun + 1)&
            + w1 * aconc(ns, izo, ime, nverti + 1)
        aconc(ns, izo, ime, nverti + 1) = dzero
      end do
    end do
  end do
  1001 format('TOP ', 3(i3.3, ' '), 6(e64.56, ' '))
  1002 format('AV ', 3(i3.3, ' '), 2(e64.56, ' '))
  !
  ! lat
  do ns = 1, nspec
    ilat = nlatbound_domain + 1
    do ivert = nverti, 1, -1
      do ime = nmerid_domain + 1, 0, -nmerid_domain - 1
        do izo = nzonal_domain, 1, -1
          ilat = ilat - 1
          
          if(species(ns)%bounddry.eq.1) then
            fdry = m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6) &
                / (m_air * sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6 &
                    + m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6))
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) * fdry
          end if
          aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert)&
              * airmloc(izo, ime, ivert) * 1d-9
          aboundlat(ilat, ns, ihourrun) = aboundlat(ilat, ns, ihourrun) &
              + w2 * aconc(ns, izo, ime, ivert)
          aboundlat(ilat, ns, ihourrun + 1) = aboundlat(ilat, ns, ihourrun + 1) &
              + w1 * aconc(ns, izo, ime, ivert)
          aconc(ns, izo, ime, ivert) = dzero
          if (abs(aboundlat(ilat, ns, ihourrun + 1)).lt.1d-100) then
            aboundlat(ilat, ns, ihourrun + 1) = dzero
          end if
          if (abs(aboundlat(ilat, ns, ihourrun)).lt.1d-100) then
            aboundlat(ilat, ns, ihourrun) = dzero
          end if
        end do
      end do
      do izo = nzonal_domain + 1, 0, -nzonal_domain - 1
        do ime = nmerid_domain, 1, -1
          ilat = ilat - 1
          
          if(species(ns)%bounddry.eq.1) then
            fdry = m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6) &
                / (m_air * sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6 &
                    + m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6))
            aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert) * fdry
          end if
          
          aconc(ns, izo, ime, ivert) = aconc(ns, izo, ime, ivert)&
              * airmloc(izo, ime, ivert) * 1d-9
          aboundlat(ilat, ns, ihourrun) = aboundlat(ilat, ns, ihourrun) &
              + w2 * aconc(ns, izo, ime, ivert)
          aboundlat(ilat, ns, ihourrun + 1) = aboundlat(ilat, ns, ihourrun + 1) &
              + w1 * aconc(ns, izo, ime, ivert)
          aconc(ns, izo, ime, ivert) = dzero
          if (abs(aboundlat(ilat, ns, ihourrun + 1)).lt.1d-100) then
            aboundlat(ilat, ns, ihourrun + 1) = dzero
          end if
          if (abs(aboundlat(ilat, ns, ihourrun)).lt.1d-100) then
            aboundlat(ilat, ns, ihourrun) = dzero
          end if
        end do
      end do
    end do
  end do
  
  !1002 format('LAT ',4(i3.3,' '),2(e64.56,' '))
end subroutine amaster_alocvalues
