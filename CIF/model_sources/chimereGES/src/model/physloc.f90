subroutine physloc

    !  Calculates various physical quantities

    use chimere_consts
    use worker_common

    implicit none

    ! local variables
    integer :: l, izo, ime, ivert, ns
    real(kind = 8) :: clwckgkg

    !  print*,rank,'ww Sets concentrations of inactive species'
    conc(nspec + 1, 1:nzonal, 1:nmerid, 1:nverti) = airmloc(1:nzonal, 1:nmerid, 1:nverti)
    conc(nspec + 2, 1:nzonal, 1:nmerid, 1:nverti) = 0.2095d0 * airmloc(1:nzonal, 1:nmerid, 1:nverti)
    conc(nspec + 3, 1:nzonal, 1:nmerid, 1:nverti) = 0.7905d0 * airmloc(1:nzonal, 1:nmerid, 1:nverti)
    conc(nspec + 4, 1:nzonal, 1:nmerid, 1:nverti) = sphuloc(1:nzonal, 1:nmerid, 1:nverti)

    !  print*,rank,'ww in or out of cloud'
    do ivert = 1, nverti
        do ime = 1, nmerid
            do izo = 1, nzonal
                clwckgkg = an * clwcloc(izo, ime, ivert) / (airmloc(izo, ime, ivert) * 29d0)
                if(clwckgkg.le.1d-6) then
                    incloud(izo, ime, ivert) = 0
                else
                    incloud(izo, ime, ivert) = 1
                endif
            enddo
        enddo
    enddo

    return
END subroutine physloc
