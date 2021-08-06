subroutine prep_outprint_tl(nh)
    
    use chimere_consts
    use worker_common
    use message_defs
    use worker_message_subs
    
    implicit none
    
    !****************************************************************************
    integer :: izo, ime, ivert
    integer :: nh, nf, ne, no, ns, nsp
    integer :: nout
    real(kind = 8) :: fppb, fdry
    real(kind = 4), allocatable, dimension(:, :, :, :) :: toprint, toprint_tl
    !*****************
    allocate(toprint(nspectot, nzonal, nmerid, nverti))
    allocate(toprint_tl(nspectot, nzonal, nmerid, nverti))
    !print*,rank,'  Calculation of families'
    do nf = 1, nfam
        nsp = nspectot - nfam + nf
        do ivert = 1, nverti
            do ime = 1, nmerid
                do izo = 1, nzonal
                    conc_tl(nsp, izo, ime, ivert) = dzero
                    conc(nsp, izo, ime, ivert) = dzero
                    do ne = 1, nelem(nf)
                        conc_tl (nsp, izo, ime, ivert) = &
                                conc_tl(nsp, izo, ime, ivert) + conc_tl(ifam(nf, ne), izo, ime, ivert)
                        conc    (nsp, izo, ime, ivert) = &
                                conc(nsp, izo, ime, ivert) + conc(ifam(nf, ne), izo, ime, ivert)
                    end do
                end do
            end do
        end do
    end do
    
    !print*,rank,'  Transformation of concentrations into ppb or microg/m3'
    do ivert = 1, nverti
        do ime = 1, nmerid
            do izo = 1, nzonal
                fppb = 1d9 / airmloc(izo, ime, ivert)
                fdry = m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6) &
                        / (m_air * sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6 &
                                + m_h2o * (1 - sphuloc(izo, ime, ivert) / airmloc(izo, ime, ivert) / 1.6))
                do nf = 1, nfam
                    nsp = nspectot - nfam + nf
                    nout = 0
                    loop : do no = 1, noutspec
                        ns = output_species(no)%iaddr
                        if(ns.eq.nsp) then
                            nout = no
                            exit loop
                        endif
                    enddo loop
                
                enddo
                do no = 1, noutspec
                    ns = output_species(no)%iaddr
                    if (species(ns)%fspec > dzero) then
                        if(ns.le.nspec) then
                            toprint_tl(no, izo, ime, ivert) = &
                                    species(ns)%fspec * conc_tl(ns, izo, ime, ivert)
                            toprint(no, izo, ime, ivert) = &
                                    species(ns)%fspec * conc(ns, izo, ime, ivert)
                        endif
                    elseif(species(ns)%fspec == dzero) then
                        toprint(no, izo, ime, ivert) = fppb / fdry &
                                * conc(ns, izo, ime, ivert)
                        toprint_tl(no, izo, ime, ivert) = fppb / fdry &
                                * conc_tl(ns, izo, ime, ivert)
                    else
                        toprint_tl(no, izo, ime, ivert) = fppb * conc_tl(ns, izo, ime, ivert)
                        toprint(no, izo, ime, ivert) = fppb * conc(ns, izo, ime, ivert)
                    endif
                enddo
            enddo
        enddo
    enddo
    
    do no = 1, noutspec
        call worker_send_toprint(toprint(no, :, :, :))
        call worker_send_toprint(toprint_tl(no, :, :, :))
    end do
    
    if(mod(nh, nsaveconcs).eq.0) then
        do no = 1, nspectot
            call worker_send_conc(no)
            call worker_send_conc_tl(no)
        end do
    end if
    
    deallocate(toprint_tl)
    deallocate(toprint)

end subroutine prep_outprint_tl
