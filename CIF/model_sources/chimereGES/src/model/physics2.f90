subroutine physics2

    !  Calculates second set of physical parameters (thereby not dependent of any species)
    !  Runs in parallelized segment

    use worker_common

    implicit none


    !*****************************************************************************
    !  Various physical calculations
    call physloc

    !  Mixing coefficients
    call mixing

    !  Horizontal transport
    call htransport

    !  Vertical transport
    call vtransport

    !  Deposition velocities
    if(ndep.ne.0)call depvel

    !  Photolysis rates
    if(ntabuzen.ne.0) call photorates

    !  Reaction rates
    call rates

end subroutine physics2
