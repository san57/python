subroutine ainiconc
  
  !  adjoint / Initialization of species concentrations
  
  use chimere_consts
  use chimere_common
  use wholedomain_common
  
  implicit none
  integer :: ip, ideb, ifin, jdeb, jfin
  integer :: izstart, izend, imstart, imend, i, j
  integer :: izo, ime, ivert, ns, iztot, imtot
  integer :: ispec
  real(kind=8) :: fdry
  real(kind=8),allocatable,dimension(:,:,:) :: fdry_array
  !******************************************************************************
  
  do ip = 1, ndoms
    i = dom(ip)%i
    j = dom(ip)%j
    izstart = dom(ip)%izstart
    izend = dom(ip)%izend
    imstart = dom(ip)%imstart
    imend = dom(ip)%imend
    ! only what is inside the domain
    if(i==1)ideb = izstart + 3
    if(i>1)ideb = izstart + (i - 1) * 6
    if(i<nzdoms)ifin = izend + (i - 1) * 6 + 6
    if(i==nzdoms)ifin = izend + (i - 1) * 6 + 3
    if(j==1)jdeb = imstart + 3
    if(j>1)jdeb = imstart + (j - 1) * 6
    if(j<nmdoms)jfin = imend + (j - 1) * 6 + 6
    if(j==nmdoms)jfin = imend + (j - 1) * 6 + 3
    do ns = 1, nspectot
      do ivert = 1, nverti
        ! Fetching inside parts of aconcsave for each sub-domain
        do iztot = ideb, ifin
          do imtot = jdeb, jfin
            izo = iztot - (i - 1) * 6 - 3
            ime = imtot - (j - 1) * 6 - 3
            fdry = 1.
            if(species(ns)%bounddry.eq.1) then
              fdry = m_h2o * (1 - sphu(izo, ime, ivert, 2) / airm(izo, ime, ivert, 2) / 1.6) &
                  / (m_air * sphu(izo, ime, ivert, 2) / airm(izo, ime, ivert, 2) / 1.6 &
                      + m_h2o * (1 - sphu(izo, ime, ivert, 2) / airm(izo, ime, ivert, 2) / 1.6))
            end if
            
            aconcini(ns, izo, ime, ivert) = aconcini(ns, izo, ime, ivert)&
                + aconcsave(ns, iztot, imtot, ivert) &
                    * airm(izo, ime, ivert, 2) * 1d-9 * fdry
          
          enddo
        enddo
      enddo
    enddo
  enddo
  
  ! Adding contribution from aconco
  allocate(fdry_array(nzonal_domain,nmerid_domain,nverti))
  do ispec = 1, nspectot
    fdry_array(:,:,:) = 1.
    if(species(ispec)%bounddry.eq.1) then
      fdry_array(:,:,:) = m_h2o * (1 - sphu(:,:,:, 2) / airm(:,:,:, 2) / 1.6) &
          / (m_air * sphu(:,:,:, 2) / airm(:,:,:, 2) / 1.6 &
              + m_h2o * (1 - sphu(:,:,:, 2) / airm(:,:,:, 2) / 1.6))
    end if
    aconcini(ispec, :, :, :) = aconcini(ispec, :, :, :)&
        + aconco(ispec, :, :, :) * airm(:, :, :, 2) * 1d-9 * fdry_array(:, :, :)
  end do
  !  aconc(:,1:nzonal,1:nmerid,1:nverti)=0.d0

end subroutine ainiconc
