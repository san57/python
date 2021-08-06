subroutine vertav(n, z, z1, z2, k1, k2, w)

    implicit none

    !*************************************************************************
    ! subroutine arguments
    integer :: n
    real(kind = 8), dimension(n) :: z
    real(kind = 8) :: z1
    real(kind = 8) :: z2
    integer :: k1
    integer :: k2
    real(kind = 8), dimension(n) :: w

    ! local variables
    integer :: k

    !**************************************************************************
    !  Finds lower and upper boundaries

    do k = 1, n - 1
        if(z(k).le.z1.and.z(k + 1).ge.z1) then
            k1 = k
            go to 1001
        endif
    enddo
    print *, '*** VERTAV.F: Z1 off the limits of Z'
    print *, '*** ', z1, z(1), z(n)
    stop
    1001 continue
    do k = 1, n - 1
        if(z(k).le.z2.and.z(k + 1).ge.z2) then
            k2 = k + 1
            go to 1002
        endif
    enddo
    print *, '*** VERTAV.F: Z2 off the limits of Z'
    print *, '*** ', z2, z(1), z(n)
    stop
    1002 continue

    !  Definition of weights

    do k = 1, n
        w(k) = 0.
    enddo

    if(k2 - k1.eq.1) then
        w(k1) = 0.5 * (2. * z(k2) - z1 - z2) * (z2 - z1) / (z(k2) - z(k1))
        w(k2) = 0.5 * (z1 + z2 - 2. * z(k1)) * (z2 - z1) / (z(k2) - z(k1))
    else
        w(k1) = 0.5 * (z(k1 + 1) - z1) * (z(k1 + 1) - z1) / (z(k1 + 1) - z(k1))
        w(k2) = 0.5 * (z2 - z(k2 - 1)) * (z2 - z(k2 - 1)) / (z(k2) - z(k2 - 1))
        w(k1 + 1) = w(k1 + 1) + 0.5 * (z(k1 + 1) - z1) * (z(k1 + 1) + z1 - 2 * z(k1)) / (z(k1 + 1) - z(k1))
        w(k2 - 1) = w(k2 - 1) + 0.5 * (z2 - z(k2 - 1)) * (2 * z(k2) - z(k2 - 1) - z2) / (z(k2) - z(k2 - 1))
        do k = k1 + 1, k2 - 2
            w(k) = w(k) + 0.5 * (z(k + 1) - z(k))
            w(k + 1) = w(k + 1) + 0.5 * (z(k + 1) - z(k))
        enddo
    endif

    do k = k1, k2
        w(k) = w(k) / (z2 - z1)
    enddo

end subroutine vertav
