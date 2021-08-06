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

SUBROUTINE fyhyp ( yzoomdeg, grossism, dzoom,tau  , &
   rrlatu,yyprimu,rrlatv,yyprimv,rlatu2,yprimu2,rlatu1,yprimu1 , &
   champmin,champmax                                            ) 

!!    ...  Version du 01/04/2001 ....

  IMPLICIT NONE
!
!    ...   Auteur :  P. Le Van  ... 
!
!    .......    d'apres  formulations  de R. Sadourny  .......
!
!     Calcule les latitudes et derivees dans la grille du GCM pour une
!     fonction f(y) a tangente  hyperbolique  .
!
!     grossism etant le grossissement ( = 2 si 2 fois, = 3 si 3 fois , etc)
!     dzoom  etant  la distance totale de la zone du zoom ( en radians )
!     tau  la raideur de la transition de l'interieur a l'exterieur du zoom   
!
!
! N.B : Il vaut mieux avoir : grossism * dzoom  <  pi/2  (radians) ,en lati.
!      ********************************************************************
!
!
 include "dimensions.h"
 include "paramet.h"

  INTEGER,PARAMETER :: nmax=30000 , nmax2=2*nmax
  !
  !
  !     .......  arguments  d'entree    .......
  !
  REAL :: yzoomdeg, grossism,dzoom,tau 
  !         ( rentres  par  run.def )
  
  !     .......  arguments  de sortie   .......
  !
  REAL :: rrlatu(jjp1), yyprimu(jjp1),rrlatv(jjm), yyprimv(jjm), &
     rlatu1(jjm), yprimu1(jjm), rlatu2(jjm), yprimu2(jjm)
  
  !
  !     .....     champs  locaux    .....
  !
  
  REAL*8 :: ylat(jjp1), yprim(jjp1)
  REAL*8 :: yuv
  REAL*8,SAVE :: Yf(0:nmax2),Ytprim(0:nmax2),yt(0:nmax2)
  REAL*8 :: fhyp(0:nmax2),beta,fxm(0:nmax2)
  REAL*8 :: yypr(0:nmax2)
  REAL*8 :: yvrai(jjp1), yprimm(jjp1),ylatt(jjp1)
  REAL*8 :: pi,depi,pis2,epsilon,y0,pisjm
  REAL*8 :: yo1,yi,ylon2,ymoy,Yprimin,champmin,champmax
  REAL*8 :: yfi,Yf1,ffdy
  REAL*8 :: ypn
  REAL*8,SAVE :: deply,y00
  
  INTEGER :: i,j,it,ik,iter,jlat
  INTEGER,SAVE :: jpn
  integer :: jjpn
  REAL*8 :: a0,a1,a2,a3,yi2,heavyy0,heavyy0m
  REAL*8 :: fa(0:nmax2),fb(0:nmax2)
  REAL :: y0min,y0max
  
  REAL*8 :: heavyside
  EXTERNAL :: heavyside
  
  pi       = 2. * ASIN(1.)
  depi     = 2. * pi
  pis2     = pi/2.
  pisjm    = pi/ FLOAT(jjm)
  epsilon  = 1.e-3
  y0       =  yzoomdeg * pi/180. 
  
  IF( dzoom<1.)  THEN
      dzoom = dzoom * pi
  ELSEIF( dzoom< 12. ) THEN
      WRITE(6,*) ' Le param. dzoomy pour fyhyp est trop petit ! L aug &
         &menter et relancer ! '
      STOP 1
  ELSE
      dzoom = dzoom * pi/180.
  ENDIF
  
  WRITE(6,18)
  WRITE(6,*) ' yzoom( rad.),grossism,tau,dzoom (radians)'
  WRITE(6,24) y0,grossism,tau,dzoom
  
  DO i = 0, nmax2 
    yt(i) = - pis2  + FLOAT(i)* pi /nmax2
  ENDDO
  
  heavyy0m = heavyside( -y0 )
  heavyy0  = heavyside(  y0 )
  y0min    = 2.*y0*heavyy0m - pis2
  y0max    = 2.*y0*heavyy0  + pis2

  fa = 999.999
  fb = 999.999
  
  DO i = 0, nmax2 
    IF( yt(i)<y0 )  THEN
        fa (i) = tau*  (yt(i)-y0+dzoom/2. )
        fb(i) =   (yt(i)-2.*y0*heavyy0m +pis2) * ( y0 - yt(i) )
    ELSEIF ( yt(i)>y0 )  THEN
        fa(i) =   tau *(y0-yt(i)+dzoom/2. )
        fb(i) = (2.*y0*heavyy0 -yt(i)+pis2) * ( yt(i) - y0 ) 
    ENDIF
    
    IF( 200.* fb(i) < - fa(i) )   THEN
        fhyp ( i) = - 1.
    ELSEIF( 200. * fb(i) < fa(i) ) THEN
        fhyp ( i) =   1.
    ELSE  
        fhyp(i) =  TANH ( fa(i)/fb(i) )
    ENDIF
    
    IF( yt(i)==y0 )  fhyp(i) = 1.
    IF(yt(i)== y0min.OR.yt(i)== y0max ) fhyp(i) = -1.
    
  ENDDO

  !!  ....  Calcul  de  beta  ....
  !
  ffdy   = 0.
  
  DO i = 1, nmax2
    ymoy    = 0.5 * ( yt(i-1) + yt( i ) )
    IF( ymoy<y0 )  THEN
        fa(i)= tau * ( ymoy-y0+dzoom/2.) 
        fb(i) = (ymoy-2.*y0*heavyy0m +pis2) * ( y0 - ymoy )
    ELSEIF ( ymoy>y0 )  THEN
        fa(i)= tau * ( y0-ymoy+dzoom/2. ) 
        fb(i) = (2.*y0*heavyy0 -ymoy+pis2) * ( ymoy - y0 )
    ENDIF
    
    IF( 200.* fb(i) < - fa(i) )    THEN
        fxm ( i) = - 1.
    ELSEIF( 200. * fb(i) < fa(i) ) THEN
        fxm ( i) =   1.
    ELSE
        fxm(i) =  TANH ( fa(i)/fb(i) )
    ENDIF
    IF( ymoy==y0 )  fxm(i) = 1.
    IF (ymoy== y0min.OR.yt(i)== y0max ) fxm(i) = -1.
    ffdy = ffdy + fxm(i) * ( yt(i) - yt(i-1) )
    
  ENDDO
  
  beta  = ( grossism * ffdy - pi ) / ( ffdy - pi )
  
  IF( 2.*beta - grossism<= 0.)  THEN
      
      WRITE(6,*) ' **  Attention ! La valeur beta calculee dans la rou &
         &tine fyhyp est mauvaise ! '
      WRITE(6,*)'Modifier les valeurs de  grossismy ,tauy ou dzoomy', &
         ' et relancer ! ***  '
      CALL ABORT
      
  ENDIF
  !
  !   .....  calcul  de  Ytprim   .....
  !
  
  DO i = 0, nmax2
    Ytprim(i) = beta  + ( grossism - beta ) * fhyp(i)
  ENDDO
  
  !   .....  Calcul  de  Yf     ........
  
  Yf(0) = - pis2
  DO i = 1, nmax2
    yypr(i)    = beta + ( grossism - beta ) * fxm(i)
  ENDDO
  
  DO i=1,nmax2
    Yf(i)   = Yf(i-1) + yypr(i) * ( yt(i) - yt(i-1) )
  ENDDO
  
  !    ****************************************************************
  !
  !   .....   yuv  = 0.   si calcul des latitudes  aux pts.  U  .....
  !   .....   yuv  = 0.5  si calcul des latitudes  aux pts.  V  .....
  !
  WRITE(6,18)
  !
  DO ik = 1,4
    
    IF( ik==1 )  THEN
        yuv  = 0.
        jlat = jjm + 1
    ELSE IF ( ik==2 )  THEN
        yuv  = 0.5
        jlat = jjm 
    ELSE IF ( ik==3 )  THEN
        yuv  = 0.25
        jlat = jjm 
    ELSE IF ( ik==4 )  THEN
        yuv  = 0.75
        jlat = jjm 
    ENDIF
    !
    yo1   = 0.
    DO j =  1,jlat
      yo1   = 0.
      ylon2 =  - pis2 + pisjm * ( FLOAT(j)  + yuv  -1.)  
      yfi    = ylon2
      !
      DO it =  nmax2,0,-1
        IF( yfi>=Yf(it))  GO TO 350
      ENDDO
      it = 0
350   CONTINUE
      
      yi = yt(it)
      IF(it==nmax2)  THEN
          it       = nmax2 -1
          Yf(it+1) = pis2
      ENDIF
      !  .................................................................
      !  ....  Interpolation entre  yi(it) et yi(it+1)   pour avoir Y(yi)  
      !      .....           et   Y'(yi)                             .....
      !  .................................................................
      
      CALL coefpoly ( Yf(it),Yf(it+1),Ytprim(it), Ytprim(it+1), &
         yt(it),yt(it+1) ,   a0,a1,a2,a3   )      
      
      
      Yf1     = Yf(it)
      Yprimin = a1 + 2.* a2 * yi + 3.*a3 * yi *yi
      
      DO iter = 1,300
        yi = yi - ( Yf1 - yfi )/ Yprimin
        
        IF( ABS(yi-yo1)<=epsilon)  GO TO 550
        yo1      = yi
        yi2      = yi * yi
        Yf1      = a0 +  a1 * yi +     a2 * yi2  +     a3 * yi2 * yi
        Yprimin  =       a1      + 2.* a2 *  yi  + 3.* a3 * yi2
      ENDDO
      WRITE(6,*) ' Pas de solution ***** ',j,ylon2,iter
      STOP 2
550   CONTINUE
      !
      Yprimin   = a1  + 2.* a2 *  yi   + 3.* a3 * yi* yi
      yprim(j)  = pi / ( jjm * Yprimin )
      yvrai(j)  = yi 
      
    ENDDO
    
    DO j = 1, jlat -1
      IF( yvrai(j+1)< yvrai(j) )  THEN
          WRITE(6,*) ' PBS. avec  rlat(',j+1,') plus petit que rlat(',j, &
             ')'
          STOP 3
      ENDIF
    ENDDO
    
    WRITE(6,*) 'Reorganisation des latitudes pour avoir entre - pi/2' &
       ,' et  pi/2 '
    !
    IF( ik==1 )   THEN
        ypn = pis2 
        DO j = jlat,1,-1
          IF( yvrai(j)<= ypn ) GO TO 1502
        ENDDO
1502    CONTINUE
        
        jpn   = j
        y00   = yvrai(jpn)
        deply = pis2 -  y00
    ENDIF
    
    DO  j = 1, jjm +1 - jpn
      ylatt (j)  = -pis2 - y00  + yvrai(jpn+j-1)
      yprimm(j)  = yprim(jpn+j-1)
    ENDDO
    
    jjpn  = jpn
    IF( jlat== jjm ) jjpn = jpn -1
    
    DO j = 1,jjpn 
      ylatt (j + jjm+1 -jpn) = yvrai(j) + deply
      yprimm(j + jjm+1 -jpn) = yprim(j)
    ENDDO
    
    !      ***********   Fin de la reorganisation     *************
    !
    
    DO j = 1, jlat
      ylat(j) =  ylatt( jlat +1 -j )
      yprim(j) = yprimm( jlat +1 -j )
    ENDDO
    
    DO j = 1, jlat
      yvrai(j) = ylat(j)*180./pi
    ENDDO
    
    IF( ik==1 )  THEN
        DO j = 1, jlat
          rrlatu(j) =  ylat( j )
          yyprimu(j) = yprim( j )
        ENDDO
        
    ELSE IF ( ik== 2 )  THEN
        
        DO j = 1, jlat
          rrlatv(j) =  ylat( j )
          yyprimv(j) = yprim( j )
        ENDDO
        
    ELSE IF ( ik== 3 )  THEN
        DO j = 1, jlat
          rlatu2(j) =  ylat( j )
          yprimu2(j) = yprim( j )
        ENDDO
        
    ELSE IF ( ik== 4 )  THEN
        
        DO j = 1, jlat
          rlatu1(j) =  ylat( j )
          yprimu1(j) = yprim( j )
        ENDDO
        
    ENDIF
    
  ENDDO
  !
  WRITE(6,18)
  !
  !  .....     fin de la boucle  do 5000 .....
  
  DO j = 1, jjm
    ylat(j) = rrlatu(j) - rrlatu(j+1)
  ENDDO
  champmin =  1.e12
  champmax = -1.e12
  DO j = 1, jjm
    champmin = MIN( champmin, ylat(j) )
    champmax = MAX( champmax, ylat(j) )
  ENDDO
  champmin = champmin * 180./pi
  champmax = champmax * 180./pi
  
24 FORMAT(2x,'Parametres yzoom,gross,tau ,dzoom pour fyhyp ',4f8.3)
18 FORMAT(/)
68 FORMAT(1x,7f9.2)
  
  RETURN
END SUBROUTINE fyhyp
