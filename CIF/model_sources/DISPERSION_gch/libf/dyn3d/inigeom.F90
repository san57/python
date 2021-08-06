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

SUBROUTINE inigeom
  !
  !     Auteur :  P. Le Van
  !
  !   ............      Version  du 01/04/2001     ........................
  !
  !  Calcul des elongations cuij1,.cuij4 , cvij1,..cvij4  aux memes en-
  !     endroits que les aires aireij1,..aireij4 .
  !
  !  Choix entre f(y) a derivee sinusoid. ou a derivee tangente hyperbol.
  !
  !
  IMPLICIT NONE
  !
 include "dimensions.h"
 include "paramet.h"
 include "comconst.h"
 include "comgeom2.h"
 include "serre.h"
 include "logic.h"
 include "comdissnew.h"

  !-----------------------------------------------------------------------
  !   ....  Variables  locales   ....
  !
  INTEGER :: i,j,itmax,itmay,iter
  REAL :: cvu(iip1,jjp1),cuv(iip1,jjm)
  REAL :: ai14,ai23,areaz,rlatp,rlatm,xprm,xprp,un4rad2,yprp,yprm
  REAL :: eps,x1,xo1,f,df,xdm,y1,yo1,ydm
  REAL :: coslatm,coslatp,radclatm,radclatp
  REAL :: cuij1(iip1,jjp1),cuij2(iip1,jjp1),cuij3(iip1,jjp1), &
     cuij4(iip1,jjp1)
  REAL :: cvij1(iip1,jjp1),cvij2(iip1,jjp1),cvij3(iip1,jjp1), &
     cvij4(iip1,jjp1)
  REAL :: rlonvv(iip1),rlatuu(jjp1)
  REAL,SAVE :: rlatu1(jjm),yprimu1(jjm),rlatu2(jjm),yprimu2(jjm) , &
     yprimv(jjm),yprimu(jjp1)
  REAL :: gamdi_gdiv, gamdi_grot, gamdi_h
  
  REAL,SAVE :: rlonm025(iip1),xprimm025(iip1), rlonp025(iip1), &
     xprimp025(iip1)

  REAL ::     SSUM
  EXTERNAL :: SSUM
  !
  !
  !   ------------------------------------------------------------------
  !   -                                                                -
  !   -    calcul des coeff. ( cu, cv , 1./cu**2,  1./cv**2  )         -
  !   -                                                                -
  !   ------------------------------------------------------------------
  !
  !      les coef. ( cu, cv ) permettent de passer des vitesses naturelles
  !      aux vitesses covariantes et contravariantes , ou vice-versa ...
  !
  !
  !     on a :  u (covariant) = cu * u (naturel)   , u(contrav)= u(nat)/cu
  !             v (covariant) = cv * v (naturel)   , v(contrav)= v(nat)/cv
  !
  !       on en tire :  u(covariant) = cu * cu * u(contravariant)
  !                     v(covariant) = cv * cv * v(contravariant)
  !
  !
  !     on a l'application (  x(X) , y(Y) )   avec - im/2 +1 <  X  < im/2
  !                                                          =     =
  !                                           et   - jm/2    <  Y  < jm/2
  !                                                          =     =
  !
  !      ...................................................
  !      ...................................................
  !      .  x  est la longitude du point  en radians       .
  !      .  y  est la  latitude du point  en radians       .
  !      .                                                 .
  !      .  on a :  cu(i,j) = rad * COS(y) * dx/dX         .
  !      .          cv( j ) = rad          * dy/dY         .
  !      .        aire(i,j) =  cu(i,j) * cv(j)             .
  !      .                                                 .
  !      . y, dx/dX, dy/dY calcules aux points concernes   .
  !      .                                                 .
  !      ...................................................
  !      ...................................................
  !
  !
  !
  !                                                           ,
  !    cv , bien que dependant de j uniquement,sera ici indice aussi en i
  !          pour un adressage plus facile en  ij  .
  !
  !
  !
  !  **************  aux points  u  et  v ,           *****************
  !      xprimu et xprimv sont respectivement les valeurs de  dx/dX
  !      yprimu et yprimv    .  .  .  .  .  .  .  .  .  .  .  dy/dY
  !      rlatu  et  rlatv    .  .  .  .  .  .  .  .  .  .  .la latitude
  !      cvu    et   cv      .  .  .  .  .  .  .  .  .  .  .    cv
  !
  !  **************  aux points u, v, scalaires, et z  ****************
  !      cu, cuv, cuscal, cuz sont respectiv. les valeurs de    cu
  !
  !
  !
  !         Exemple de distribution de variables sur la grille dans le
  !             domaine de travail ( X,Y ) .
  !     ................................................................
  !                  DX=DY= 1
  !
  !   
  !        +     represente  un  point scalaire ( p.exp  la pression )
  !        >     represente  la composante zonale du  vent
  !        V     represente  la composante meridienne du vent
  !        o     represente  la  vorticite
  !
  !       ----  , car aux poles , les comp.zonales covariantes sont nulles
  !
  !
  !
  !         i ->
  !
  !         1      2      3      4      5      6      7      8
  !  j
  !  v  1   + ---- + ---- + ---- + ---- + ---- + ---- + ---- + --
  !
  !         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
  !
  !     2   +   >  +   >  +   >  +   >  +   >  +   >  +   >  +  >
  !
  !         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
  !
  !     3   +   >  +   >  +   >  +   >  +   >  +   >  +   >  +  >
  !
  !         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
  !
  !     4   +   >  +   >  +   >  +   >  +   >  +   >  +   >  +  >
  !
  !         V   o  V   o  V   o  V   o  V   o  V   o  V   o  V  o
  !
  !     5   + ---- + ---- + ---- + ---- + ---- + ---- + ---- + --
  !
  !
  !      Ci-dessus,  on voit que le nombre de pts.en longitude est egal
  !                 a   IM = 8
  !      De meme ,   le nombre d'intervalles entre les 2 poles est egal
  !                 a   JM = 4
  !
  !      Les points scalaires ( + ) correspondent donc a des valeurs
  !       entieres  de  i ( 1 a IM )   et  de  j ( 1 a  JM +1 )   .
  !
  !      Les vents    U       ( > ) correspondent a des valeurs  semi-
  !       entieres  de i ( 1+ 0.5 a IM+ 0.5) et entieres de j ( 1 a JM+1)
  !
  !      Les vents    V       ( V ) correspondent a des valeurs entieres
  !       de     i ( 1 a  IM ) et semi-entieres de  j ( 1 +0.5  a JM +0.5)
  !
  !
  !
  WRITE(6,3) 
3 FORMAT( // 10x,' ....  INIGEOM  date du 01/06/98   ..... ', &
     //5x,'   Calcul des elongations cu et cv  comme sommes des 4 ' / &
     5x,' elong. cuij1, .. 4  , cvij1,.. 4  qui les entourent , aux &
     &'/ 5x,' memes endroits que les aires aireij1,...j4   . ' / )
  !
  !
  IF( nitergdiv/=2 ) THEN
      gamdi_gdiv = coefdis/ ( float(nitergdiv) -2. )
  ELSE
      gamdi_gdiv = 0.
  ENDIF
  IF( nitergrot/=2 ) THEN
      gamdi_grot = coefdis/ ( float(nitergrot) -2. )
  ELSE
      gamdi_grot = 0.
  ENDIF
  IF( niterh/=2 ) THEN
      gamdi_h = coefdis/ ( float(niterh) -2. )
  ELSE
      gamdi_h = 0.
  ENDIF

  WRITE(6,*) ' gamdi_gd ',gamdi_gdiv,gamdi_grot,gamdi_h,coefdis, &
     nitergdiv,nitergrot,niterh
  !
  pi    = 2.* ASIN(1.)
  !
  WRITE(6,990) 
  
  !     ----------------------------------------------------------------
  !
  IF( .NOT.fxyhypb )   THEN
      !
      !
      IF( ysinus )  THEN
          !
          WRITE(6,*) ' ***  Inigeom ,  Y = Sinus ( Latitude ) *** '
          !
          !   .... utilisation de f(x,y )  avec  y  =  sinus de la latitude  .....
          
          CALL  fxysinus (rlatu,yprimu,rlatv,yprimv,rlatu1,yprimu1, &
             rlatu2,yprimu2, &
             rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)
          
      ELSE
          !
          WRITE(6,*) '*** Inigeom ,  Y = Latitude  , der. sinusoid . ***'
          
          !  .... utilisation  de f(x,y) a tangente sinusoidale , y etant la latit. ...
          !
          
          pxo   = clon *pi /180.
          pyo   = 2.* clat* pi /180.
          !
          !  ....  determination de  transx ( pour le zoom ) par Newton-Raphson ...
          !
          itmax = 10
          eps   = .1e-7
          !
          xo1 = 0.
          DO iter = 1, itmax
            x1  = xo1
            f   = x1+ alphax *SIN(x1-pxo)
            df  = 1.+ alphax *COS(x1-pxo)
            x1  = x1 - f/df
            xdm = ABS( x1- xo1 )
            IF( xdm<=eps ) EXIT
            xo1 = x1
          ENDDO
          !
          transx = xo1
          
          itmay = 10
          eps   = .1e-7
          !
          yo1  = 0.
          DO iter = 1,itmay
            y1   = yo1
            f    = y1 + alphay* SIN(y1-pyo)
            df   = 1. + alphay* COS(y1-pyo)
            y1   = y1 -f/df
            ydm  = ABS(y1-yo1)
            IF(ydm<=eps) EXIT
            yo1  = y1
          ENDDO
          !
          transy = yo1
          
          CALL fxy ( rlatu,yprimu,rlatv,yprimv,rlatu1,yprimu1, &
             rlatu2,yprimu2, &
             rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025)
          
      ENDIF
      !
  ELSE
      !
      !   ....  Utilisation  de fxyhyper , f(x,y) a derivee tangente hyperbol.
      !   .....................................................................
      
      WRITE(6,*)'*** Inigeom , Y = Latitude  , der.tg. hyperbolique ***'
      
      CALL fxyhyper( clat, grossismy, dzoomy, tauy    , & 
         clon, grossismx, dzoomx, taux    , &
         rlatu,yprimu,rlatv, yprimv,rlatu1, yprimu1,rlatu2,yprimu2  , &
         rlonu,xprimu,rlonv,xprimv,rlonm025,xprimm025,rlonp025,xprimp025 )
      
      
  ENDIF
  !
  !  -------------------------------------------------------------------

  !
  rlatu(1)    =     ASIN(1.)
  rlatu(jjp1) =  - rlatu(1)
  !
  !
  !   ....  calcul  aux  poles  ....
  !
  yprimu(1)      = 0.
  yprimu(jjp1)   = 0.
  !
  !
  un4rad2 = 0.25 * rad * rad
  !
  !   --------------------------------------------------------------------
  !   --------------------------------------------------------------------
  !   -                                                                  -
  !   -  calcul  des aires ( aire,areau,airev, 1./aire, 1./airez  )      -
  !   -      et de   fext ,  force de coriolis  extensive  .             -
  !   -                                                                  -
  !   --------------------------------------------------------------------
  !   --------------------------------------------------------------------
  !
  !
  !
  !   A 1 point scalaire P (i,j) de la grille, reguliere en (X,Y) , sont
  !   affectees 4 aires entourant P , calculees respectivement aux points
  !            ( i + 1/4, j - 1/4 )    :    aireij1 (i,j)
  !            ( i + 1/4, j + 1/4 )    :    aireij2 (i,j)
  !            ( i - 1/4, j + 1/4 )    :    aireij3 (i,j)
  !            ( i - 1/4, j - 1/4 )    :    aireij4 (i,j)
  !
  !           ,
  !   Les cotes de chacun de ces 4 carres etant egaux a 1/2 suivant (X,Y).
  !   Chaque aire centree en 1 point scalaire P(i,j) est egale a la somme
  !   des 4 aires  aireij1,aireij2,aireij3,aireij4 qui sont affectees au
  !   point (i,j) .
  !   On definit en outre les coefficients  alpha comme etant egaux a
  !    (aireij / aire), c.a.d par exp.  alpha1(i,j)=aireij1(i,j)/aire(i,j)
  !
  !   De meme, toute aire centree en 1 point U est egale a la somme des
  !   4 aires aireij1,aireij2,aireij3,aireij4 entourant le point U .
  !    Idem pour  airev, airez .
  !
  !       On a ,pour chaque maille :    dX = dY = 1
  !
  !
  !                             . V
  !
  !                 aireij4 .        . aireij1
  !
  !                   U .       . P      . U
  !
  !                 aireij3 .        . aireij2
  !
  !                             . V
  !
  !
  !
  !
  !
  !   ....................................................................
  !
  !    Calcul des 4 aires elementaires aireij1,aireij2,aireij3,aireij4
  !   qui entourent chaque aire(i,j) , ainsi que les 4 elongations elemen
  !   taires cuij et les 4 elongat. cvij qui sont calculees aux memes 
  !     endroits  que les aireij   .    
  !
  !   ....................................................................
  !
  !     .......  do 35  :   boucle sur les  jjm + 1  latitudes   .....
  !
  !
  DO j = 1, jjp1
    !
    IF ( j== 1 )  THEN
        !
        yprm           = yprimu1(j)
        rlatm          = rlatu1(j)
        !
        coslatm        = COS( rlatm )
        radclatm       = 0.5* rad * coslatm
        !
        DO i = 1, iim
          xprp           = xprimp025( i )
          xprm           = xprimm025( i )
          areaij2( i,1 ) = un4rad2 * coslatm  * xprp * yprm
          areaij3( i,1 ) = un4rad2 * coslatm  * xprm * yprm
          cuij2  ( i,1 ) = radclatm * xprp
          cuij3  ( i,1 ) = radclatm * xprm
          cvij2  ( i,1 ) = 0.5* rad * yprm
          cvij3  ( i,1 ) = cvij2(i,1)
        ENDDO
        !
        DO  i = 1, iim
          areaij1( i,1 ) = 0.
          areaij4( i,1 ) = 0.
          cuij1  ( i,1 ) = 0.
          cuij4  ( i,1 ) = 0.
          cvij1  ( i,1 ) = 0.
          cvij4  ( i,1 ) = 0.
        ENDDO
        !
    END IF
    !
    IF ( j== jjp1 )  THEN
        yprp               = yprimu2(j-1)
        rlatp              = rlatu2 (j-1)
        !
        coslatp             = COS( rlatp )
        radclatp            = 0.5* rad * coslatp
        !
        DO i = 1,iim
          xprp              = xprimp025( i )
          xprm              = xprimm025( i )
          areaij1( i,jjp1 ) = un4rad2 * coslatp  * xprp * yprp
          areaij4( i,jjp1 ) = un4rad2 * coslatp  * xprm * yprp
          cuij1(i,jjp1)     = radclatp * xprp
          cuij4(i,jjp1)     = radclatp * xprm
          cvij1(i,jjp1)     = 0.5 * rad* yprp
          cvij4(i,jjp1)     = cvij1(i,jjp1)
        ENDDO
        !
        DO   i    = 1, iim
          areaij2( i,jjp1 ) = 0.
          areaij3( i,jjp1 ) = 0.
          cvij2  ( i,jjp1 ) = 0.
          cvij3  ( i,jjp1 ) = 0.
          cuij2  ( i,jjp1 ) = 0.
          cuij3  ( i,jjp1 ) = 0.
        ENDDO
        !
    END IF
    !
    
    IF ( j > 1 .AND. j < jjp1 )  THEN
        !
        rlatp    = rlatu2 ( j-1 )
        yprp     = yprimu2( j-1 )
        rlatm    = rlatu1 (  j  )
        yprm     = yprimu1(  j  )
        
        coslatm  = COS( rlatm )
        coslatp  = COS( rlatp )
        radclatp = 0.5* rad * coslatp
        radclatm = 0.5* rad * coslatm
        !
        DO  i = 1,iim
          xprp            = xprimp025( i )
          xprm            = xprimm025( i )
          
          ai14            = un4rad2 * coslatp * yprp
          ai23            = un4rad2 * coslatm * yprm
          areaij1 ( i,j ) = ai14 * xprp
          areaij2 ( i,j ) = ai23 * xprp
          areaij3 ( i,j ) = ai23 * xprm
          areaij4 ( i,j ) = ai14 * xprm
          cuij1   ( i,j ) = radclatp * xprp
          cuij2   ( i,j ) = radclatm * xprp
          cuij3   ( i,j ) = radclatm * xprm
          cuij4   ( i,j ) = radclatp * xprm
          cvij1   ( i,j ) = 0.5* rad * yprp
          cvij2   ( i,j ) = 0.5* rad * yprm
          cvij3   ( i,j ) = cvij2(i,j)
          cvij4   ( i,j ) = cvij1(i,j)
        ENDDO
        !
    END IF
    !
    !    ........       periodicite   ............
    !
    cvij1   (iip1,j) = cvij1   (1,j)
    cvij2   (iip1,j) = cvij2   (1,j)
    cvij3   (iip1,j) = cvij3   (1,j)
    cvij4   (iip1,j) = cvij4   (1,j)
    cuij1   (iip1,j) = cuij1   (1,j)
    cuij2   (iip1,j) = cuij2   (1,j)
    cuij3   (iip1,j) = cuij3   (1,j)
    cuij4   (iip1,j) = cuij4   (1,j)
    areaij1 (iip1,j) = areaij1 (1,j )
    areaij2 (iip1,j) = areaij2 (1,j )
    areaij3 (iip1,j) = areaij3 (1,j )
    areaij4 (iip1,j) = areaij4 (1,j )
    
  ENDDO
  !
  !    ..............................................................
  !
  DO j = 1, jjp1
    DO i = 1, iim
      area    ( i,j )  = areaij1(i,j) + areaij2(i,j) + areaij3(i,j) + &
         areaij4(i,j)
      alpha1  ( i,j )  = areaij1(i,j) / area(i,j)
      alpha2  ( i,j )  = areaij2(i,j) / area(i,j)
      alpha3  ( i,j )  = areaij3(i,j) / area(i,j)
      alpha4  ( i,j )  = areaij4(i,j) / area(i,j)
      alpha1p2( i,j )  = alpha1 (i,j) + alpha2 (i,j)
      alpha1p4( i,j )  = alpha1 (i,j) + alpha4 (i,j)
      alpha2p3( i,j )  = alpha2 (i,j) + alpha3 (i,j)
      alpha3p4( i,j )  = alpha3 (i,j) + alpha4 (i,j)
    ENDDO
    !
    !
    area    (iip1,j) = area    (1,j)
    alpha1  (iip1,j) = alpha1  (1,j)
    alpha2  (iip1,j) = alpha2  (1,j)
    alpha3  (iip1,j) = alpha3  (1,j)
    alpha4  (iip1,j) = alpha4  (1,j)
    alpha1p2(iip1,j) = alpha1p2(1,j)
    alpha1p4(iip1,j) = alpha1p4(1,j)
    alpha2p3(iip1,j) = alpha2p3(1,j)
    alpha3p4(iip1,j) = alpha3p4(1,j)
  ENDDO
  !
  
  DO j = 1,jjp1
    DO i = 1,iim
      areau       (i,j)= areaij1(i,j) + areaij2(i,j) + areaij4(i+1,j) + &
         areaij3(i+1,j)
      unsarea    ( i,j)= 1./ area(i,j)
      unsair_gam1( i,j)= unsarea(i,j)** ( - gamdi_gdiv )
      unsair_gam2( i,j)= unsarea(i,j)** ( - gamdi_h    )
      areasurg   ( i,j)= area(i,j)/ g
    ENDDO
    areau     (iip1,j)  = areau  (1,j)
    unsarea   (iip1,j)  = unsarea(1,j)
    unsair_gam1(iip1,j) = unsair_gam1(1,j)
    unsair_gam2(iip1,j) = unsair_gam2(1,j)
    areasurg   (iip1,j) = areasurg(1,j)
  ENDDO
!
!
  DO j = 1,jjm
    !
    DO i=1,iim
      areav     (i,j) = areaij2(i,j)+ areaij3(i,j)+ areaij1(i,j+1) + &
         areaij4(i,j+1)
    ENDDO
    DO i=1,iim
      areaz         = areaij2(i,j)+areaij1(i,j+1)+areaij3(i+1,j) + &
         areaij4(i+1,j+1)
      unsareaz(i,j) = 1./ areaz
      unsairz_gam(i,j)= unsareaz(i,j)** ( - gamdi_grot )
      fext    (i,j)   = areaz * SIN(rlatv(j))* 2.* omeg
    ENDDO
    areav     (iip1,j)  = areav(1,j)
    unsareaz  (iip1,j)  = unsareaz(1,j)
    fext      (iip1,j)  = fext(1,j)
    unsairz_gam(iip1,j) = unsairz_gam(1,j)
    !
  ENDDO
  !
  !
  !    .....      Calcul  des elongations cu,cv, cvu     .........
  !
  DO    j   = 1, jjm
    DO   i  = 1, iim
      cv(i,j) = 0.5 *( cvij2(i,j)+cvij3(i,j)+cvij1(i,j+1)+cvij4(i,j+1))
      cvu(i,j)= 0.5 *( cvij1(i,j)+cvij4(i,j)+cvij2(i,j)  +cvij3(i,j) )
      cuv(i,j)= 0.5 *( cuij2(i,j)+cuij3(i,j)+cuij1(i,j+1)+cuij4(i,j+1))
      unscv2(i,j) = 1./ ( cv(i,j)*cv(i,j) )
    ENDDO
    DO   i  = 1, iim
      cuvsurcv (i,j)    = areav(i,j)  * unscv2(i,j)
      cvsurcuv (i,j)    = 1./cuvsurcv(i,j)
      cuvscvgam1(i,j)   = cuvsurcv (i,j) ** ( - gamdi_gdiv )
      cuvscvgam2(i,j)   = cuvsurcv (i,j) ** ( - gamdi_h )
      cvscuvgam(i,j)    = cvsurcuv (i,j) ** ( - gamdi_grot )
    ENDDO
    cv       (iip1,j)  = cv       (1,j)
    cvu      (iip1,j)  = cvu      (1,j)
    unscv2   (iip1,j)  = unscv2   (1,j)
    cuv      (iip1,j)  = cuv      (1,j)
    cuvsurcv (iip1,j)  = cuvsurcv (1,j)
    cvsurcuv (iip1,j)  = cvsurcuv (1,j)
    cuvscvgam1(iip1,j) = cuvscvgam1(1,j)
    cuvscvgam2(iip1,j) = cuvscvgam2(1,j)
    cvscuvgam(iip1,j)  = cvscuvgam(1,j)
  ENDDO
  
  DO  j     = 2, jjm
    DO   i  = 1, iim
      cu(i,j) = 0.5*(cuij1(i,j)+cuij4(i+1,j)+cuij2(i,j)+cuij3(i+1,j))
      unscu2    (i,j)  = 1./ ( cu(i,j) * cu(i,j) )
      cvusurcu  (i,j)  =  areau(i,j) * unscu2(i,j)
      cusurcvu  (i,j)  = 1./ cvusurcu(i,j)
      cvuscugam1 (i,j) = cvusurcu(i,j) ** ( - gamdi_gdiv ) 
      cvuscugam2 (i,j) = cvusurcu(i,j) ** ( - gamdi_h    ) 
      cuscvugam (i,j)  = cusurcvu(i,j) ** ( - gamdi_grot )
    ENDDO
    cu       (iip1,j)  = cu(1,j)
    unscu2   (iip1,j)  = unscu2(1,j)
    cvusurcu (iip1,j)  = cvusurcu(1,j)
    cusurcvu (iip1,j)  = cusurcvu(1,j)
    cvuscugam1(iip1,j) = cvuscugam1(1,j)
    cvuscugam2(iip1,j) = cvuscugam2(1,j)
    cuscvugam (iip1,j) = cuscvugam(1,j)
  ENDDO
  
  !
  !   ....  calcul aux  poles  ....
  !
  DO    i      =  1, iip1
    cu    ( i, 1 )  =   0.
    unscu2( i, 1 )  =   0.
    cvu   ( i, 1 )  =   0.
    !
    cu    (i, jjp1) =   0.
    unscu2(i, jjp1) =   0.
    cvu   (i, jjp1) =   0.
  ENDDO
  !
  !    ..............................................................
  !
  DO j = 1, jjm
    DO i= 1, iim
      airvscu2  (i,j) = areav(i,j)/ ( cuv(i,j) * cuv(i,j) )
      aivscu2gam(i,j) = airvscu2(i,j)** ( - gamdi_grot )
    ENDDO
    airvscu2  (iip1,j)  = airvscu2(1,j)
    aivscu2gam(iip1,j)  = aivscu2gam(1,j)
  ENDDO
  
  DO j=2,jjm
    DO i=1,iim
      airuscv2   (i,j)    = areau(i,j)/ ( cvu(i,j) * cvu(i,j) )
      aiuscv2gam (i,j)    = airuscv2(i,j)** ( - gamdi_grot ) 
    ENDDO
    airuscv2  (iip1,j)  = airuscv2  (1,j)
    aiuscv2gam(iip1,j)  = aiuscv2gam(1,j)
  ENDDO

  !
  !   calcul des aires aux  poles :
  !   -----------------------------
  !
  apoln       = SSUM(iim,area(1,1),1)
  apols       = SSUM(iim,area(1,jjp1),1)
  unsapolnga1 = 1./ ( apoln ** ( - gamdi_gdiv ) )
  unsapolsga1 = 1./ ( apols ** ( - gamdi_gdiv ) )
  unsapolnga2 = 1./ ( apoln ** ( - gamdi_h    ) )
  unsapolsga2 = 1./ ( apols ** ( - gamdi_h    ) )
!
!-----------------------------------------------------------------------
!     gtitre='Coriolis version ancienne'
!     gfichier='fext1'
!     CALL writestd(fext,iip1*jjm)
!
!   changement F. Hourdin calcul conservatif pour fext
!   constang contient le produit a * cos ( latitude ) * omega
!
  DO i=1,iim
    constang(i,1) = 0.
  ENDDO
  DO j=1,jjm-1
    DO i=1,iim
      constang(i,j+1) = rad*omeg*cu(i,j+1)*COS(rlatu(j+1))
    ENDDO
  ENDDO
  DO i=1,iim
    constang(i,jjp1) = 0.
  ENDDO
  !
  !   periodicite en longitude
  !
  DO j=1,jjm
    fext(iip1,j)     = fext(1,j)
  ENDDO
  DO j=1,jjp1
    constang(iip1,j) = constang(1,j)
  ENDDO

  ! fin du changement

  !
  !-----------------------------------------------------------------------
  !
  WRITE(6,*) '   ***  Coordonnees de la grille  *** '
  WRITE(6,995)
  !
  WRITE(6,*) '   LONGITUDES  aux pts.   V  ( degres )  '
  WRITE(6,995)
  DO i=1,iip1
    rlonvv(i) = rlonv(i)*180./pi
  ENDDO
  WRITE(6,400) rlonvv
  !
  WRITE(6,995)
  WRITE(6,*) '   LATITUDES   aux pts.   V  ( degres )  '
  WRITE(6,995)
  DO i=1,jjm
    rlatuu(i)=rlatv(i)*180./pi
  ENDDO
  WRITE(6,400) (rlatuu(i),i=1,jjm)
  !
  DO i=1,iip1
    rlonvv(i)=rlonu(i)*180./pi
  ENDDO
  WRITE(6,995)
  WRITE(6,*) '   LONGITUDES  aux pts.   U  ( degres )  '
  WRITE(6,995)
  WRITE(6,400) rlonvv
  WRITE(6,995)
  
  WRITE(6,*) '   LATITUDES   aux pts.   U  ( degres )  '
  WRITE(6,995)
  DO i=1,jjp1
    rlatuu(i)=rlatu(i)*180./pi
  ENDDO
  WRITE(6,400) (rlatuu(i),i=1,jjp1)
  WRITE(6,995)
  !
444 FORMAT(f10.3,f6.0)
400 FORMAT(1x,8f8.2)
990 FORMAT(//)
995 FORMAT(/)
  !
  RETURN
END SUBROUTINE inigeom
