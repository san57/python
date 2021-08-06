subroutine aphysics2

    !  Calculates second set of physical parameters (thereby not dependent of any species)
    !  Runs in parallelized segment

    use worker_common

    implicit none


    !*****************************************************************************
    !  Various physical calculations
    call aphysloc

end subroutine aphysics2
