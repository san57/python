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

SUBROUTINE inifilr
  !
  !    ... H. Upadhyaya, O.Sharma   ...
  !
  IMPLICIT NONE
  !
  !     version 3 .....
  
  !     Correction  le 28/10/97    P. Le Van .
  !  -------------------------------------------------------------------
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "parafilt.h"
  !  -------------------------------------------------------------------
  INCLUDE "comgeom.h"
  INCLUDE "coefils.h"
  INCLUDE "logic.h"
  INCLUDE "serre.h"
  
  REAL :: dlonu(iim),dlatu(jjm)
  REAL :: rlamda( iim ),  eignvl( iim )
  !
  
  REAL ::   lamdamax,pi,cof
  INTEGER :: i,j,modemax,imx,k,kf,ii
  REAL :: dymin,dxmin,colat0
  REAL :: eignft(iim,iim), coff
  REAL :: matriceun,matriceus,matricevn,matricevs,matrinvn,matrinvs
  COMMON/matrfil/matriceun(iim,iim,nfilun),matriceus(iim,iim,nfilus) &
     , matricevn(iim,iim,nfilvn),matricevs(iim,iim,nfilvs) &
     ,  matrinvn(iim,iim,nfilun),matrinvs (iim,iim,nfilus)
#ifdef CRAY
  INTEGER ::  ISMIN
  EXTERNAL :: ISMIN
  INTEGER :: iymin 
  INTEGER :: ixmineq
#endif
  EXTERNAL :: inifgn
  !
  ! ------------------------------------------------------------
  !   This routine computes the eigenfunctions of the laplacien
  !   on the stretched grid, and the filtering coefficients
  !      
  !  We designate:
  !   eignfn   eigenfunctions of the discrete laplacien
  !   eigenvl  eigenvalues
  !   jfiltn   indexof the last scalar line filtered in NH
  !   jfilts   index of the first line filtered in SH
  !   modfrst  index of the mode from WHERE modes are filtered
  !   modemax  maximum number of modes ( im )
  !   coefil   filtering coefficients ( lamda_max*COS(rlat)/lamda )
  !   sdd      SQRT( dx )
  !      
  !     the modes are filtered from modfrst to modemax
  !      
  !-----------------------------------------------------------
  !
  
  pi       = 2. * ASIN( 1. )
  
  DO i = 1,iim
    dlonu(i) = xprimu( i )
  ENDDO
  !
  CALL inifgn(eignvl)
  !
  PRINT *,' EIGNVL '
  PRINT 250,eignvl
250 FORMAT( 1x,5e13.6)
  !
  ! compute eigenvalues and eigenfunctions
  !
  !
  !.................................................................
  !
  !  compute the filtering coefficients for scalar lines and 
  !  meridional wind v-lines
  !
  !  we filter all those latitude lines WHERE coefil < 1
  !  NO FILTERING AT POLES
  !
  !  colat0 is to be used  when alpha (stretching coefficient)
  !  is set equal to zero for the regular grid CASE 
  !
  !    .......   Calcul  de  colat0   .........
  !     .....  colat0 = minimum de ( 0.5, min dy/ min dx )   ...
  !
  !
  DO j = 1,jjm
    dlatu( j ) = rlatu( j ) - rlatu( j+1 )
  ENDDO
  !
#ifdef CRAY
  iymin   = ISMIN( jjm, dlatu, 1 )
  ixmineq = ISMIN( iim, dlonu, 1 )
  dymin   = dlatu( iymin )
  dxmin   = dlonu( ixmineq )
#else
  dxmin   =  dlonu(1)
  DO  i  = 2, iim
    dxmin = MIN( dxmin,dlonu(i) )
  ENDDO
  dymin  = dlatu(1)
  DO j  = 2, jjm
    dymin = MIN( dymin,dlatu(j) )
  ENDDO
#endif
  !
  !
  colat0  =  MIN( 0.5, dymin/dxmin )
  !
  IF( .NOT.fxyhypb.AND.ysinus )  THEN
      colat0 = 0.6
      !         ...... a revoir  pour  ysinus !   .......
      alphax = 0.
  ENDIF
  !
  PRINT 50, colat0,alphax
50 FORMAT(/15x,' Inifilr colat0 alphax ',2e16.7)
  !
  IF(alphax==1. )  THEN
      PRINT *,' Inifilr  alphax doit etre  <  a 1.  Corriger '
      STOP
  ENDIF
  !
  lamdamax = iim / ( pi * colat0 * ( 1. - alphax ) )
  
  !c                        ... Correction  le 28/10/97  ( P.Le Van ) ..
  !
  DO i = 2,iim
    rlamda( i ) = lamdamax/ SQRT( ABS( eignvl(i) ) )
  ENDDO
  !
    
    DO j = 1,jjm
      DO i = 1,iim
        coefilu( i,j )  = 0.0
        coefilv( i,j )  = 0.0
        coefilu2( i,j ) = 0.0
        coefilv2( i,j ) = 0.0
      ENDDO
    ENDDO
        
    !
    !    ... Determination de jfiltnu,jfiltnv,jfiltsu,jfiltsv ....
    !    .........................................................
    !
    modemax = iim
    
    !ccc    imx = modemax - 4 * (modemax/iim)
    
    imx  = iim
    !
    PRINT *,' TRUNCATION AT ',imx
    !
    DO j = 2, jjm/2+1
      cof = COS( rlatu(j) )/ colat0
      IF ( cof < 1. ) THEN
          IF( rlamda(imx) * COS(rlatu(j) )<1. ) jfiltnu= j
      ENDIF
      
      cof = COS( rlatu(jjp1-j+1) )/ colat0
      IF ( cof < 1. ) THEN
          IF( rlamda(imx) * COS(rlatu(jjp1-j+1) )<1. ) &
             jfiltsu= jjp1-j+1
      ENDIF
    ENDDO
    !
    DO j = 1, jjm/2
      cof = COS( rlatv(j) )/ colat0
      IF ( cof < 1. ) THEN
          IF( rlamda(imx) * COS(rlatv(j) )<1. ) jfiltnv= j
      ENDIF
      
      cof = COS( rlatv(jjm-j+1) )/ colat0
      IF ( cof < 1. ) THEN
          IF( rlamda(imx) * COS(rlatv(jjm-j+1) )<1. ) &
             jfiltsv= jjm-j+1
      ENDIF
    ENDDO
    !                                 
    
    IF( jfiltnu<=0 .OR. jfiltnu> jjm/2 +1 )  THEN
        PRINT *,' jfiltnu en dehors des valeurs acceptables ' ,jfiltnu
        STOP
    ENDIF
    
    IF( jfiltsu<=0 .OR. jfiltsu>  jjm  +1 )  THEN
        PRINT *,' jfiltsu en dehors des valeurs acceptables ' ,jfiltsu
        STOP
    ENDIF
    
    IF( jfiltnv<=0 .OR. jfiltnv> jjm/2    )  THEN
        PRINT *,' jfiltnv en dehors des valeurs acceptables ' ,jfiltnv
        STOP
    ENDIF
    
    IF( jfiltsv<=0 .OR. jfiltsv>     jjm  )  THEN
        PRINT *,' jfiltsv en dehors des valeurs acceptables ' ,jfiltsv
        STOP
    ENDIF
    
    PRINT *,' jfiltnv jfiltsv jfiltnu jfiltsu ' , &
       jfiltnv,jfiltsv,jfiltnu,jfiltsu
    
    !                                 
    !   ... Determination de coefilu,coefilv,n=modfrstu,modfrstv ....
    !................................................................
    !
    !
    DO j = 1,jjm
      modfrstu( j ) = iim
      modfrstv( j ) = iim
    ENDDO
    !
    DO j = 2,jfiltnu
      DO k = 2,modemax
        cof = rlamda(k) * COS( rlatu(j) )
        IF ( cof < 1. ) EXIT
      ENDDO
      IF (k/=modemax+1) THEN
        modfrstu( j ) = k
        !
        kf = modfrstu( j )
        DO k = kf , modemax
          cof = rlamda(k) * COS( rlatu(j) )
          coefilu(k,j) = cof - 1.
          coefilu2(k,j) = cof*cof - 1.
        ENDDO
    ENDIF
  ENDDO
  !                                 
  !
  DO j = 1,jfiltnv
    !
    DO k = 2,modemax
      cof = rlamda(k) * COS( rlatv(j) )
      IF ( cof < 1. ) EXIT
    ENDDO
    IF (k/=modemax+1) THEN
        
        modfrstv( j ) = k
        !
        kf = modfrstv( j )
        DO k = kf , modemax
          cof = rlamda(k) * COS( rlatv(j) )
          coefilv(k,j) = cof - 1.
          coefilv2(k,j) = cof*cof - 1.
        ENDDO
    ENDIF
    !
  ENDDO
  !
  DO j = jfiltsu,jjm
    DO k = 2,modemax
      cof = rlamda(k) * COS( rlatu(j) )
      IF ( cof < 1. ) EXIT
    ENDDO
    IF (k/=modemax+1) THEN
        modfrstu( j ) = k
        !
        kf = modfrstu( j )
        DO k = kf , modemax
          cof = rlamda(k) * COS( rlatu(j) )
          coefilu(k,j) = cof - 1.
          coefilu2(k,j) = cof*cof - 1.
        ENDDO
    ENDIF
  ENDDO
  !                                 
  DO j = jfiltsv,jjm
    DO k = 2,modemax
      cof = rlamda(k) * COS( rlatv(j) )
      IF ( cof < 1. ) EXIT
    ENDDO
    IF (k/=modemax+1) THEN
        modfrstv( j ) = k
        !
        kf = modfrstv( j )
        DO k = kf , modemax
          cof = rlamda(k) * COS( rlatv(j) )
          coefilv(k,j) = cof - 1.
          coefilv2(k,j) = cof*cof - 1.
        ENDDO
    ENDIF
  ENDDO
  !
  
  IF(jfiltnv>=jjm/2 .OR. jfiltnu>=jjm/2)THEN
      
      IF(jfiltnv==jfiltsv)jfiltsv=1+jfiltnv
      IF(jfiltnu==jfiltsu)jfiltsu=1+jfiltnu
      
      PRINT *,'jfiltnv jfiltsv jfiltnu jfiltsu' , &
         jfiltnv,jfiltsv,jfiltnu,jfiltsu
  ENDIF
  
  PRINT *,'   Modes premiers  v  '
  PRINT 334,modfrstv
  PRINT *,'   Modes premiers  u  '
  PRINT 334,modfrstu
  
  
  IF( nfilun< jfiltnu )  THEN
      PRINT *,' le parametre nfilun utilise pour la matrice ', &
         ' matriceun  est trop petit ! ' 
      PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnu
      PRINT *,' Pour information, nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      STOP
  ENDIF
  IF( nfilun> jfiltnu+ 2 )  THEN
      PRINT *,' le parametre nfilun utilise pour la matrice ', &
         ' matriceun est trop grand ! Gachis de memoire ! ' 
      PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnu
      PRINT *,' Pour information, nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      !              STOP
  ENDIF
  IF( nfilus< jjm - jfiltsu +1 )  THEN
      PRINT *,' le parametre nfilus utilise pour la matrice ', &
         ' matriceus  est trop petit !  '
      PRINT *,' Le changer dans parafilt.h et le mettre a  ', &
         jjm - jfiltsu + 1
      PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      STOP
  ENDIF
  IF( nfilus> jjm - jfiltsu + 3 )  THEN
      PRINT *,' le parametre nfilus utilise pour la matrice ', &
         ' matriceus  est trop grand ! ' 
      PRINT *,' Le changer dans parafilt.h et le mettre a  ' , &
         jjm - jfiltsu + 1
      PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      !              STOP
  ENDIF
  IF( nfilvn< jfiltnv )  THEN
      PRINT *,' le parametre nfilvn utilise pour la matrice ', &
         ' matricevn  est trop petit ! '  
      PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnv
      PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      STOP
  ENDIF
  IF( nfilvn> jfiltnv+ 2 )  THEN
      PRINT *,' le parametre nfilvn utilise pour la matrice ', &
         ' matricevn est trop grand !  Gachis de memoire ! ' 
      PRINT *,'Le changer dans parafilt.h et le mettre a  ',jfiltnv
      PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      !              STOP
  ENDIF
  IF( nfilvs< jjm - jfiltsv +1 )  THEN
      PRINT *,' le parametre nfilvs utilise pour la matrice ', &
         ' matricevs  est trop petit !  Le changer dans parafilt.h '
      PRINT *,' Le changer dans parafilt.h et le mettre a  ' &
         , jjm - jfiltsv + 1
      PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      STOP
  ENDIF
  IF( nfilvs> jjm - jfiltsv + 3 )  THEN
      PRINT *,' le parametre nfilvs utilise pour la matrice ', &
         ' matricevs  est trop grand ! Gachis de memoire ! '
      PRINT *,' Le changer dans parafilt.h et le mettre a  ' &
         ,  jjm - jfiltsv + 1
      PRINT *,' Pour information , nfilun,nfilus,nfilvn,nfilvs ' &
         ,'doivent etre egaux successivement a  ',jfiltnu,jjm-jfiltsu+1 &
         ,jfiltnv,jjm-jfiltsv+1
      !              STOP
  ENDIF

  !  
  !   ...................................................................
  !
  !   ... Calcul de la matrice filtre 'matriceu'  pour les champs situes
  !                       sur la grille scalaire                 ........
  !   ...................................................................
  !
  DO j = 2, jfiltnu
    
    DO i=1,iim
      coff = coefilu(i,j)
      IF( i<modfrstu(j) ) coff = 0.
      DO k=1,iim
        eignft(i,k) = eignfnv(k,i) * coff
      ENDDO
    ENDDO
#ifdef CRAY
    CALL MXM( eignfnv,iim,eignft,iim,matriceun(1,1,j),iim )
#else
#ifdef BLAS
    CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
       eignfnv, iim, eignft, iim, 0.0, matriceun(1,1,j), iim)
#else
    DO k = 1, iim
      DO i = 1, iim
        matriceun(i,k,j) = 0.0
        DO ii = 1, iim
          matriceun(i,k,j) = matriceun(i,k,j) &
             + eignfnv(i,ii)*eignft(ii,k)
        ENDDO
      ENDDO
    ENDDO
#endif
#endif
    
  ENDDO
  
  DO j = jfiltsu, jjm
    
    DO i=1,iim
      coff = coefilu(i,j)
      IF( i<modfrstu(j) ) coff = 0.
      DO k=1,iim
        eignft(i,k) = eignfnv(k,i) * coff
      ENDDO
    ENDDO
#ifdef CRAY
    CALL MXM(eignfnv,iim,eignft,iim,matriceus(1,1,j-jfiltsu+1),iim)
#else
#ifdef BLAS 
    CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
       eignfnv, iim, eignft, iim, 0.0,  &
       matriceus(1,1,j-jfiltsu+1), iim)
#else
    DO k = 1, iim
      DO i = 1, iim
        matriceus(i,k,j-jfiltsu+1) = 0.0
        DO ii = 1, iim
          matriceus(i,k,j-jfiltsu+1) = matriceus(i,k,j-jfiltsu+1) &
             + eignfnv(i,ii)*eignft(ii,k)
        ENDDO
      ENDDO
    ENDDO
#endif
#endif
    
  ENDDO
  
  !   ...................................................................
  !
  !   ... Calcul de la matrice filtre 'matricev'  pour les champs situes
  !                       sur la grille   de V ou de Z           ........
  !   ...................................................................
  !
  DO j = 1, jfiltnv
    
    DO i = 1, iim
      coff = coefilv(i,j)
      IF( i<modfrstv(j) ) coff = 0.
      DO k = 1, iim
        eignft(i,k) = eignfnu(k,i) * coff
      ENDDO
    ENDDO
#ifdef CRAY
    CALL MXM( eignfnu,iim,eignft,iim,matricevn(1,1,j),iim )
#else
#ifdef BLAS
    CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
       eignfnu, iim, eignft, iim, 0.0, matricevn(1,1,j), iim)
#else
    DO k = 1, iim
      DO i = 1, iim
        matricevn(i,k,j) = 0.0
        DO ii = 1, iim
          matricevn(i,k,j) = matricevn(i,k,j) &
             + eignfnu(i,ii)*eignft(ii,k)
        ENDDO
      ENDDO
    ENDDO
#endif
#endif
    
  ENDDO
  
  DO j = jfiltsv, jjm
    
    DO i = 1, iim
      coff = coefilv(i,j)
      IF( i<modfrstv(j) ) coff = 0.
      DO k = 1, iim
        eignft(i,k) = eignfnu(k,i) * coff
      ENDDO
    ENDDO
#ifdef CRAY
    CALL MXM(eignfnu,iim,eignft,iim,matricevs(1,1,j-jfiltsv+1),iim)
#else
#ifdef BLAS
    CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
       eignfnu, iim, eignft, iim, 0.0,  &
       matricevs(1,1,j-jfiltsv+1), iim)
#else
    DO k = 1, iim
      DO i = 1, iim
        matricevs(i,k,j-jfiltsv+1) = 0.0
        DO ii = 1, iim
          matricevs(i,k,j-jfiltsv+1) = matricevs(i,k,j-jfiltsv+1) &
             + eignfnu(i,ii)*eignft(ii,k)
        ENDDO
      ENDDO
    ENDDO
#endif
#endif
    
  ENDDO
  
  !   ...................................................................
  !
  !   ... Calcul de la matrice filtre 'matrinv'  pour les champs situes
  !              sur la grille scalaire , pour le filtre inverse ........
  !   ...................................................................
  !
  DO j = 2, jfiltnu
    
    DO i = 1,iim
      coff = coefilu(i,j)/ ( 1. + coefilu(i,j) )
      IF( i<modfrstu(j) ) coff = 0.
      DO k=1,iim
        eignft(i,k) = eignfnv(k,i) * coff
      ENDDO
    ENDDO
#ifdef CRAY
    CALL MXM( eignfnv,iim,eignft,iim,matrinvn(1,1,j),iim )
#else
#ifdef BLAS
    CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
       eignfnv, iim, eignft, iim, 0.0, matrinvn(1,1,j), iim)
#else
    DO k = 1, iim
      DO i = 1, iim
        matrinvn(i,k,j) = 0.0
        DO ii = 1, iim
          matrinvn(i,k,j) = matrinvn(i,k,j) &
             + eignfnv(i,ii)*eignft(ii,k)
        ENDDO
      ENDDO
    ENDDO
#endif
#endif
    
  ENDDO
  
  DO j = jfiltsu, jjm
    
    DO i = 1,iim
      coff = coefilu(i,j) / ( 1. + coefilu(i,j) )
      IF( i<modfrstu(j) ) coff = 0.
      DO k=1,iim
        eignft(i,k) = eignfnv(k,i) * coff
      ENDDO
    ENDDO
#ifdef CRAY
    CALL MXM(eignfnv,iim,eignft,iim,matrinvs(1,1,j-jfiltsu+1),iim)
#else
#ifdef BLAS
    CALL SGEMM ('N', 'N', iim, iim, iim, 1.0, &
       eignfnv, iim, eignft, iim, 0.0, matrinvs(1,1,j-jfiltsu+1), iim)
#else
    DO k = 1, iim
      DO i = 1, iim
        matrinvs(i,k,j-jfiltsu+1) = 0.0
        DO ii = 1, iim
          matrinvs(i,k,j-jfiltsu+1) = matrinvs(i,k,j-jfiltsu+1) &
             + eignfnv(i,ii)*eignft(ii,k)
        ENDDO
      ENDDO
    ENDDO
#endif
#endif
    
  ENDDO
  
  !   ...................................................................
  
  !
334 FORMAT(1x,24i3)
755 FORMAT(1x,6f10.3,i3)
  
  RETURN
END SUBROUTINE inifilr
