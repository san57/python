SUBROUTINE gr_ecrit_fi(nfield,nlon,iiim,jjmp1,ecrit,fi)

IMPLICIT none

! Transforms a variable from the writing grid to the physical grid

INTEGER nfield,nlon,iiim,jjmp1, jjjm
REAL fi(nlon,nfield), ecrit(iiim,jjmp1,nfield)
INTEGER i, j, n, ig

jjjm = jjmp1 - 1
DO n = 1, nfield
  fi(1,n) = ecrit(1,1,n)
  fi(nlon,n) = ecrit(1,jjmp1,n)
  DO j = 2, jjjm
    ig = 2+(j-2)*iiim
    DO i = 1, iiim
      fi(ig-1+i,n) = ecrit(i,j,n)
    ENDDO
  ENDDO
ENDDO
RETURN
END
