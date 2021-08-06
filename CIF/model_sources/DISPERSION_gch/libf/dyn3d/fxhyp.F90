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

SUBROUTINE fxhyp ( xzoomdeg,grossism,dzoom,tau , &
   rlonm025,xprimm025,rlonv,xprimv,rlonu,xprimu,rlonp025,xprimp025, &
   champmin,champmax                                               )
  
  !      Auteur :  P. Le Van 
  
  IMPLICIT NONE
  
  !   Calcule les longitudes et derivees dans la grille du GCM pour une
  !    fonction f(x) a tangente  hyperbolique  .
  !
  !    grossism etant le grossissement ( = 2 si 2 fois, = 3 si 3 fois,etc.)
  !    dzoom  etant  la distance totale de la zone du zoom
  !    tau  la raideur de la transition de l'interieur a l'exterieur du zoom
  !
  !   On doit avoir grossism x dzoom <  pi ( radians )   , en longitude.
  !  ********************************************************************
  
  
  INTEGER,PARAMETER ::  nmax=30000, nmax2=2*nmax
  !
  LOGICAL,PARAMETER :: scal180=.TRUE.
  
  !      scal180 = .TRUE.  si on veut avoir le premier point scalaire pour   
  !      une grille reguliere ( grossism = 1.,tau=0.,clon=0. ) a -180. degres.
  !      sinon scal180 = .FALSE.
  
 include "dimensions.h"
 include "paramet.h"
  
  !     ......  arguments  d'entree   .......
  !
  REAL :: xzoomdeg,dzoom,tau,grossism
  
  !    ......   arguments  de  sortie  ......
  
  REAL :: rlonm025(iip1),xprimm025(iip1),rlonv(iip1),xprimv(iip1), &
     rlonu(iip1),xprimu(iip1),rlonp025(iip1),xprimp025(iip1)
  
  !     .... variables locales  ....
  !
  REAL*8 :: xlon(iip1),xprimm(iip1),xuv
  REAL*8 :: xtild(0:nmax2)
  REAL*8 :: fhyp(0:nmax2),ffdx,beta,Xprimt(0:nmax2)
  REAL*8 :: Xf(0:nmax2),xxpr(0:nmax2)
  REAL*8 :: xvrai(iip1),xxprim(iip1) 
  REAL*8 :: pi,depi,epsilon,xzoom,fa,fb
  REAL*8 :: Xf1, Xfi , a0,a1,a2,a3,xi2
  INTEGER :: i,it,ik,iter,ii,idif,ii1,ii2
  REAL*8 :: xi,xo1,xmoy,xlon2,fxm,Xprimin
  REAL*8 :: champmin,champmax,decalx
  INTEGER,save :: is2
  
  REAL*8 :: heavyside
  EXTERNAL :: coefpoly,heavyside
  
  pi       = 2. * ASIN(1.)
  depi     = 2. * pi
  epsilon  = 1.e-3
  xzoom    = xzoomdeg * pi/180. 
  !
  decalx   = .75
  IF( grossism==1..AND.scal180 )  THEN
      decalx   = 1.
  ENDIF
  
  WRITE(6,*) 'FXHYP scal180,decalx', scal180,decalx
  !
  IF( dzoom<1.)  THEN
      dzoom = dzoom * depi
  ELSEIF( dzoom< 25. ) THEN
      WRITE(6,*) ' Le param. dzoomy pour fxhyp est trop petit ! L aug & 
         &menter et relancer ! '
      STOP 1
  ELSE
      dzoom = dzoom * pi/180.
  ENDIF
  
  WRITE(6,*) ' xzoom( rad.),grossism,tau,dzoom (radians)'
  WRITE(6,24) xzoom,grossism,tau,dzoom
  
  DO i = 0, nmax2 
    xtild(i) = - pi + FLOAT(i) * depi /nmax2
  ENDDO
  
  DO i = nmax, nmax2
    
    fa  = tau*  ( dzoom/2.  - xtild(i) )
    fb  = xtild(i) *  ( pi - xtild(i) )
    
    IF( 200.* fb < - fa )   THEN
        fhyp ( i) = - 1.
    ELSEIF( 200. * fb < fa ) THEN
        fhyp ( i) =   1.
    ELSE
        IF( ABS(fa)<1.e-13.AND.ABS(fb)<1.e-13)  THEN
            IF(   200.*fb + fa<1.e-10 )  THEN
                fhyp ( i ) = - 1.
            ELSEIF( 200.*fb - fa<1.e-10 )  THEN
                fhyp ( i )  =   1.
            ENDIF
        ELSE
            fhyp ( i )  =  TANH ( fa/fb )
        ENDIF
    ENDIF
    
    IF ( xtild(i)== 0. )  fhyp(i) =  1.
    IF ( xtild(i)== pi )  fhyp(i) = -1.
    
  ENDDO
  
  !!  ....  Calcul  de  beta  ....
  !   ............................
  
  ffdx = 0.
  
  DO i = nmax +1,nmax2
    
    xmoy    = 0.5 * ( xtild(i-1) + xtild( i ) )
    fa  = tau*  ( dzoom/2.  - xmoy )
    fb  = xmoy *  ( pi - xmoy )
    
    IF( 200.* fb < - fa )   THEN
        fxm = - 1.
    ELSEIF( 200. * fb < fa ) THEN
        fxm =   1.
    ELSE
        IF( ABS(fa)<1.e-13.AND.ABS(fb)<1.e-13)  THEN
            IF(   200.*fb + fa<1.e-10 )  THEN
                fxm   = - 1.
            ELSEIF( 200.*fb - fa<1.e-10 )  THEN
                fxm   =   1.
            ENDIF
        ELSE
            fxm   =  TANH ( fa/fb )
        ENDIF
    ENDIF
    
    IF ( xmoy== 0. )  fxm  =  1.
    IF ( xmoy== pi )  fxm  = -1.
    
    ffdx = ffdx + fxm * ( xtild(i) - xtild(i-1) )
    
  ENDDO
  
  beta  = ( grossism * ffdx - pi ) / ( ffdx - pi )
  
  IF( 2.*beta - grossism< 0.)  THEN
      WRITE(6,*) ' **  Attention ! La valeur beta calculee dans la rou &
         &tine fxhyp est mauvaise ! '
      WRITE(6,*)'Modifier les valeurs de  grossismx ,tau ou dzoomx ', &
         ' et relancer ! ***  '
      CALL ABORT
  ENDIF
  !
  !   .....  calcul  de  Xprimt   .....
  !
  
  DO i = nmax, nmax2
    Xprimt(i) = beta  + ( grossism - beta ) * fhyp(i)
  ENDDO
  !   
  DO i =  nmax+1, nmax2
    Xprimt( nmax2 - i ) = Xprimt( i )
  ENDDO
  !
  
  !   .....  Calcul  de  Xf     ........
  
  Xf(0) = - pi
  
  DO i =  nmax +1, nmax2
    
    xmoy    = 0.5 * ( xtild(i-1) + xtild( i ) )
    fa  = tau*  ( dzoom/2.  - xmoy )
    fb  = xmoy *  ( pi - xmoy )
    
    IF( 200.* fb < - fa )   THEN
        fxm = - 1.
    ELSEIF( 200. * fb < fa ) THEN
        fxm =   1.
    ELSE
        fxm =  TANH ( fa/fb )
    ENDIF
    
    IF ( xmoy== 0. )  fxm =  1.
    IF ( xmoy== pi )  fxm = -1.
    xxpr(i)    = beta + ( grossism - beta ) * fxm
    
  ENDDO
  
  DO i = nmax+1, nmax2
    xxpr(nmax2-i+1) = xxpr(i)
  ENDDO
  
  DO i=1,nmax2
    Xf(i)   = Xf(i-1) + xxpr(i) * ( xtild(i) - xtild(i-1) )
  ENDDO
  
  
  !    *****************************************************************
  !
  
  !     .....  xuv = 0.   si  calcul  aux pts   scalaires   ........
  !     .....  xuv = 0.5  si  calcul  aux pts      U        ........
  !
  WRITE(6,18)
  !
  DO ik = 1, 4
    
    IF( ik==1 )        THEN
        xuv =  -0.25
    ELSE IF ( ik==2 )  THEN
        xuv =   0.
    ELSE IF ( ik==3 )  THEN
        xuv =   0.50
    ELSE IF ( ik==4 )  THEN
        xuv =   0.25
    ENDIF
    
    xo1   = 0.
    
    ii1=1
    ii2=iim
    IF(ik==1.AND.grossism==1.) THEN
        ii1 = 2 
        ii2 = iim+1
    ENDIF
    
    DO i = ii1, ii2
      
      xlon2 = - pi + (FLOAT(i) + xuv - decalx) * depi / FLOAT(iim) 
      
      Xfi    = xlon2
      !
      DO it =  nmax2,0,-1
        IF( Xfi>=Xf(it))  GO TO 350
      ENDDO
        
      it = 0
        
350   CONTINUE
        
      !    ......  Calcul de   Xf(xi)    ...... 
      !
      xi  = xtild(it)
      
      IF(it==nmax2)  THEN
          it       = nmax2 -1
          Xf(it+1) = pi
      ENDIF
      !  .....................................................................
      !
      !   Appel de la routine qui calcule les coefficients a0,a1,a2,a3 d'un
      !   polynome de degre 3  qui passe  par les points (Xf(it),xtild(it) )
      !          et (Xf(it+1),xtild(it+1) )
      
      CALL coefpoly ( Xf(it),Xf(it+1),Xprimt(it),Xprimt(it+1), &
         xtild(it),xtild(it+1),  a0, a1, a2, a3  )
      
      Xf1     = Xf(it)
      Xprimin = a1 + 2.* a2 * xi + 3.*a3 * xi *xi
      
      DO iter = 1,300
        xi = xi - ( Xf1 - Xfi )/ Xprimin
        
        IF( ABS(xi-xo1)<=epsilon)  GO TO 550
        xo1      = xi
        xi2      = xi * xi
        Xf1      = a0 +  a1 * xi +     a2 * xi2  +     a3 * xi2 * xi
        Xprimin  =       a1      + 2.* a2 *  xi  + 3.* a3 * xi2
      ENDDO
      WRITE(6,*) ' Pas de solution ***** ',i,xlon2,iter
      STOP 6
550   CONTINUE
      
      xxprim(i) = depi/ ( FLOAT(iim) * Xprimin )
      xvrai(i)  =  xi + xzoom
      
    ENDDO


    IF(ik==1.AND.grossism==1.)  THEN
        xvrai(1)    = xvrai(iip1)-depi
        xxprim(1)   = xxprim(iip1)
    ENDIF
    DO i = 1 , iim
      xlon(i)     = xvrai(i)
      xprimm(i)   = xxprim(i)
    ENDDO
    
    DO i = 1, iim -1
      IF( xvrai(i+1)< xvrai(i) )  THEN
          WRITE(6,*) ' PBS. avec rlonu(',i+1,') plus petit que rlonu(',i, &
             ')'
          STOP 7
      ENDIF
    ENDDO
    !
    !   ... Reorganisation  des  longitudes  pour les avoir  entre - pi et pi ..
    !   ........................................................................
    
    champmin =  1.e12
    champmax = -1.e12
    DO i = 1, iim
      champmin = MIN( champmin,xvrai(i) )
      champmax = MAX( champmax,xvrai(i) )
    ENDDO
    
    IF(champmin >=-pi-0.10.AND.champmax<=pi+0.10 )  THEN
        GO TO 1600
    ELSE
        WRITE(6,*) 'Reorganisation des longitudes pour avoir entre - pi', &
           ' et pi '
        !
        IF( xzoom<=0.)  THEN
            IF( ik== 1 )  THEN
                DO i = 1, iim
                  IF( xvrai(i)>= - pi )  GO TO 80
                ENDDO
                WRITE(6,*)  ' PBS. 1 !  Xvrai plus petit que  - pi ! '
                STOP 8
80              CONTINUE
                is2 = i
            ENDIF

            IF( is2/= 1 )  THEN
                DO ii = is2 , iim
                  xlon  (ii-is2+1) = xvrai(ii)
                  xprimm(ii-is2+1) = xxprim(ii)
                ENDDO
                DO ii = 1 , is2 -1
                  xlon  (ii+iim-is2+1) = xvrai(ii) + depi
                  xprimm(ii+iim-is2+1) = xxprim(ii) 
                ENDDO
            ENDIF
        ELSE 
            IF( ik==1 )  THEN
                DO i = iim,1,-1
                  IF( xvrai(i)<= pi ) GO TO 90
                ENDDO
                WRITE(6,*) ' PBS.  2 ! Xvrai plus grand  que   pi ! '
                STOP 9
90              CONTINUE
                is2 = i
            ENDIF
            idif = iim -is2
            DO ii = 1, is2
              xlon  (ii+idif) = xvrai(ii)
              xprimm(ii+idif) = xxprim(ii)
            ENDDO
            DO ii = 1, idif
              xlon (ii)  = xvrai (ii+is2) - depi
              xprimm(ii) = xxprim(ii+is2) 
            ENDDO
        ENDIF
    ENDIF
    !
    !     .........   Fin  de la reorganisation   ............................

1600 CONTINUE

    
    xlon  ( iip1)  = xlon(1) + depi
    xprimm( iip1 ) = xprimm (1 )
    
    DO i = 1, iim+1
      xvrai(i) = xlon(i)*180./pi
    ENDDO
    
    IF( ik==1 )  THEN
        DO i = 1,iim +1
          rlonm025(i) = xlon( i )
          xprimm025(i) = xprimm(i)
        ENDDO
        
    ELSE IF( ik==2 )  THEN
        DO i = 1,iim + 1
          rlonv(i) = xlon( i )
          xprimv(i) = xprimm(i)
        ENDDO
        
    ELSE IF( ik==3)   THEN
        DO i = 1,iim + 1
          rlonu(i) = xlon( i )
          xprimu(i) = xprimm(i)
        ENDDO
        
    ELSE IF( ik==4 )  THEN
        DO i = 1,iim + 1
          rlonp025(i) = xlon( i )
          xprimp025(i) = xprimm(i)
        ENDDO
        
    ENDIF
    
  ENDDO
  !
  WRITE(6,18)
  !
  !    ...........  fin  de la boucle  DO 5000      ............
  
  DO i = 1, iim
    xlon(i) = rlonv(i+1) - rlonv(i)
  ENDDO
  champmin =  1.e12
  champmax = -1.e12
  DO i = 1, iim
    champmin = MIN( champmin, xlon(i) )
    champmax = MAX( champmax, xlon(i) )
  ENDDO
  champmin = champmin * 180./pi
  champmax = champmax * 180./pi
  
18 FORMAT(/)
24 FORMAT(2x,'Parametres xzoom,gross,tau ,dzoom pour fxhyp ',4f8.3)
68 FORMAT(1x,7f9.2)
566 FORMAT(1x,7f9.4)
  
  RETURN
END SUBROUTINE fxhyp
