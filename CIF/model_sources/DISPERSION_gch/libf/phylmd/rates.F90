subroutine comp_rates(temp, pmid, refjrates, rate)


    USE SPECIES_NAME
    USE CONSTANTS
    USE variables
    USE dimphy
    IMPLICIT NONE
    INCLUDE "dimensions.h"
    INCLUDE "paramet.h"

    !*****************************************************************************************
    REAL(kind=8),parameter :: Rg = 8.205d-2
    REAL(kind=8),parameter :: conv = 1d+3
    REAL(kind=8),parameter :: mm_HNO3 = 6.3d+1

    INTEGER :: i,ns,nr,nt, nrj
    INTEGER :: ity, ij
    
    REAL, DIMENSION(klon, klev)   :: mat_ones

    REAL, INTENT(out) :: rate(klon, klev, nreac)
    REAL, INTENT(in)  :: temp(klon, klev)
    REAL, INTENT(in)  :: pmid(klon, klev)
    REAL, INTENT(in)  :: refjrates(klon,klev,ijratesmax)

    !*****************************************************************************************
    rate = 0.
    mat_ones = 1.

    !  Computation of rates
    DO nr=1,nreac
        ity = idtyperate(nr)
        
        IF (ity == 0) THEN
            PRINT*, 'Chemical parsing problem'
        END IF
        
        !  Constant rates
        IF (ity == 1) THEN
            rate(:, :, nr) = tabrate(1,nr)

            !  Arrhenius simplified rates
        ELSE IF (ity == 2) THEN
            rate(:, :, nr) = tabrate(1,nr)*EXP(-tabrate(2,nr)/temp(:,:))

            !  Arrhenius complete rates
        ELSE IF (ity == 3) THEN
            rate(:, :, nr) = tabrate(1,nr)*EXP(-tabrate(2,nr)/temp(:,:))*(300d0/temp(:,:))**tabrate(3,nr)

            !  Pressure rates
        ELSE IF (ity == 4) THEN
            rate(:, :, nr) = tabrate(1,nr)*(tabrate(2,nr) * mat_ones(:, :) + tabrate(3,nr) * pmid(:,:) / p_ref)

            !  Photolysis rates
        ELSE IF (ity == 5) THEN
            DO ij = 1, ijratesmax
                IF (jrates(ij)%idreac == nr) THEN
                    nrj = ij
                END IF
            END DO
            rate(:, :, nr) = refjrates(:, :, nrj)

        END IF
    END DO

end subroutine comp_rates
