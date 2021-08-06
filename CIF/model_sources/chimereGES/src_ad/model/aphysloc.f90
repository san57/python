subroutine aphysloc

    !  Calculates various physical quantities

    use chimere_consts
    use worker_common

    implicit none


    !  Sets concentrations of inactive species
    aconc(nspec + 1, 1:nzonal, 1:nmerid, 1:nverti) = 0.d0
    aconc(nspec + 2, 1:nzonal, 1:nmerid, 1:nverti) = 0.d0
    aconc(nspec + 3, 1:nzonal, 1:nmerid, 1:nverti) = 0.d0
    aconc(nspec + 4, 1:nzonal, 1:nmerid, 1:nverti) = 0.d0

    return
END subroutine aphysloc
