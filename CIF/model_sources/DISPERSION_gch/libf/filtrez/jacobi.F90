! ==========================================================================
! LMDZT - Atmospheric tracer transport, forward, tangent-linear and adjoint
!             for use within the PYVAR inversion system
!
! Copyright Laboratoire des Sciences du Climat et de l'Environnement (LSCE)
!           Unite mixte CEA-CNRS-UVSQ
! Initial version from the LMDZ.3.3 model, developed by IPSL
!
! Code manager:
! Frederic Chevallier, LSCE, frederic.chevallier@lsce.ipsl.fr
!
! This software is governed by the CeCILL  license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and,  more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
! ==================================================================

SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
  INTEGER,PARAMETER :: NMAX=400
  REAL :: A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
  IF (n>nmax) THEN
      PRINT*, 'n, nmax=', n, nmax
      PRINT*, 'Surdimensionnement insuffisant dans jacobi'
      CALL abort
  ENDIF
  DO IP=1,N
    DO IQ=1,N
      V(IP,IQ)=0.
    ENDDO
    V(IP,IP)=1.
  ENDDO
  DO IP=1,N
    B(IP)=A(IP,IP)
    D(IP)=B(IP)
    Z(IP)=0.
  ENDDO
  NROT=0
  DO I=1,50
    SM=0.
    DO IP=1,N-1
      DO IQ=IP+1,N
        SM=SM+ABS(A(IP,IQ))
      ENDDO
    ENDDO
    IF(SM==0.)RETURN
    IF(I<4)THEN
        TRESH=0.2*SM/N**2
    ELSE
        TRESH=0.
    ENDIF
    DO IP=1,N-1
      DO IQ=IP+1,N
        G=100.*ABS(A(IP,IQ))
        IF((I>4).AND.(ABS(D(IP))+G==ABS(D(IP))) &
           .AND.(ABS(D(IQ))+G==ABS(D(IQ))))THEN
            A(IP,IQ)=0.
        ELSE IF(ABS(A(IP,IQ))>TRESH)THEN
            H=D(IQ)-D(IP)
            IF(ABS(H)+G==ABS(H))THEN
                T=A(IP,IQ)/H
            ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA<0.)T=-T
            ENDIF
            C=1./SQRT(1+T**2)
            S=T*C
            TAU=S/(1.+C)
            H=T*A(IP,IQ)
            Z(IP)=Z(IP)-H
            Z(IQ)=Z(IQ)+H
            D(IP)=D(IP)-H
            D(IQ)=D(IQ)+H
            A(IP,IQ)=0.
            DO J=1,IP-1
              G=A(J,IP)
              H=A(J,IQ)
              A(J,IP)=G-S*(H+G*TAU)
              A(J,IQ)=H+S*(G-H*TAU)
            ENDDO
            DO J=IP+1,IQ-1
              G=A(IP,J)
              H=A(J,IQ)
              A(IP,J)=G-S*(H+G*TAU)
              A(J,IQ)=H+S*(G-H*TAU)
            ENDDO
            DO J=IQ+1,N
              G=A(IP,J)
              H=A(IQ,J)
              A(IP,J)=G-S*(H+G*TAU)
              A(IQ,J)=H+S*(G-H*TAU)
            ENDDO
            DO J=1,N
              G=V(J,IP)
              H=V(J,IQ)
              V(J,IP)=G-S*(H+G*TAU)
              V(J,IQ)=H+S*(G-H*TAU)
            ENDDO
            NROT=NROT+1
        ENDIF
      ENDDO
    ENDDO
    DO IP=1,N
      B(IP)=B(IP)+Z(IP)
      D(IP)=B(IP)
      Z(IP)=0.
    ENDDO
  ENDDO
  STOP '50 iterations should never happen'
  RETURN
END SUBROUTINE JACOBI
