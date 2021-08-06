subroutine awrite_concs_hour(iop, idout)

    use worker_common

    implicit none
    !-------------------------------------------------------------------
    ! Management of concentrations tables at the PHYSICAL TIME STEP
    !
    ! 'iop' is a flag for file management
    !      [1] = opening files
    !      [2] = direct mode  (write results)
    !      [3] = adjoint mode (read  results) and backspace
    !      [4] = last integration, closing files
    !-------------------------------------------------------------------
    integer iop, ide, ns
    integer ivert, izo, ime, ispec
    integer ityp, idout, k, idtmp, idadj1
    integer recordlen
    character*132 fnconc
    character*2 idoutstr
    !-------------------------------------------------------------------
    ! first time step: opening temporary concentrations save file
    if(iop.eq.1)then
        write(idoutstr, '(i2.2)')idout
        fnconc = 'TMPconchour' // idoutstr
        !call opfi(ifnconch,fnconc,'u','n')
        recordlen = 2 * nspectot * (nverti * nmerid * nzonal + (nverti + 1) * (nmerid + 6) * (nzonal + 6))
        !print*,'RECORDLEN',recordlen
        call opfi(ifnconch, fnconc, 'd', 'n', recordlen)
    endif

    ! direct print
    if(iop.eq.2)then
        !	  write(ifnconch)ihourrun,&
        write(ifnconch, rec = ihourrun + 1)&
                ((((conco(ispec, izo, ime, ivert), ivert = 1, nverti)&
                        , ime = 1, nmerid), izo = 1, nzonal), ispec = 1, nspectot), &
                ((((conc(ispec, izo, ime, ivert), ivert = 1, nverti + 1), &
                        ime = -2, nmerid + 3), izo = -2, nzonal + 3), ispec = 1, nspectot)
    endif

    ! adjoint read: loop to ensure to read the right table
    if(iop.eq.3)then
        !         rewind(ifnconch)
        !         do k=1,1000000000
        !            read(ifnconch)idtmp
        !	    !print*,idtmp
        !       	    if(idtmp.gt.idout)then
        !	       backspace(ifnconch)
        !	       backspace(ifnconch)
        !	    endif
        !	    if(idtmp.eq.idout)then
        !	       backspace(ifnconch)
        !               read(ifnconch)idadj1,&
        read(ifnconch, rec = idout + 1)&
                ((((conco(ispec, izo, ime, ivert), ivert = 1, nverti)&
                        , ime = 1, nmerid), izo = 1, nzonal), ispec = 1, nspectot), &
                ((((conc(ispec, izo, ime, ivert), ivert = 1, nverti + 1), &
                        ime = -2, nmerid + 3), izo = -2, nzonal + 3), ispec = 1, nspectot)
        !               goto 10
        !	    endif
        !	 enddo
        !10       continue
    endif


    ! Last adjoint time step: closing files
    if(iop.eq.4)close(ifnconch)

    !--------------------------------------------------
    return
end subroutine awrite_concs_hour
