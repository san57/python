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

SUBROUTINE vlsplt_ad_p(q,pente_max,masse,w,pbaru,pbarv,pdt,q_ad)
  !
  !     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
  !                TL = F. Chevallier 
  !
  !    ********************************************************************
  !     Shema  d'advection " pseudo amont " .
  !    ********************************************************************
  !     q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
  !
  !   pente_max facteur de limitation des pentes: 2 en general
  !                                               0 pour un schema amont
  !   pbaru,pbarv,w flux de masse en u ,v ,w
  !   pdt pas de temps
  !
  !   --------------------------------------------------------------------
  USE parallel
  USE mod_hallo
  USE Vampir
  USE Write_field_p
  IMPLICIT NONE
  !
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "logic.h"
  INCLUDE "comvert.h"
  INCLUDE "comconst.h"
  
  !
  !   Arguments:
  !   ----------
  REAL :: masse(ip1jmp1,llm),pente_max
  REAL :: pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm)
  REAL :: q(ip1jmp1,llm), q_ad(ip1jmp1,llm)
  REAL :: w(ip1jmp1,llm),pdt
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l
  !
  REAL :: zm(ip1jmp1,llm)
  REAL :: zm1(ip1jmp1,llm), zq1(ip1jmp1,llm)
  REAL :: zm2(ip1jmp1,llm), zq2(ip1jmp1,llm)
  REAL :: zm3(ip1jmp1,llm), zq3(ip1jmp1,llm)
  REAL :: zm4(ip1jmp1,llm), zq4(ip1jmp1,llm)
  REAL :: zm5(ip1jmp1,llm), zq5(ip1jmp1,llm)
  REAL :: mu(ip1jmp1,llm)
  REAL :: mv(ip1jm,llm)
  REAL :: mw(ip1jmp1,llm+1)
  REAL :: zq(ip1jmp1,llm)
  REAL :: zq_ad(ip1jmp1,llm)
  REAL :: zzpbar, zzw
  
  INTEGER :: ijb,ije
  TYPE(request),SAVE :: MyRequest1
  TYPE(request),SAVE :: MyRequest2
  
  
  CALL SetTag(MyRequest1,100)
  CALL SetTag(MyRequest2,101)
  
  zzpbar = 0.5 * pdt
  zzw    = pdt
  
  ijb=ij_begin
  ije=ij_end
  IF (north_pole) ijb=ijb+iip1
  IF (south_pole)  ije=ije-iip1
  
  DO l=1,llm
    DO ij = ijb,ije
      mu(ij,l)=pbaru(ij,l) * zzpbar
    ENDDO
  ENDDO
  
  
  ijb=ij_begin-iip1
  ije=ij_end
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  DO l=1,llm
    DO ij=ijb,ije
      mv(ij,l)=pbarv(ij,l) * zzpbar
    ENDDO
  ENDDO
  
  ijb=ij_begin
  ije=ij_end
  
  DO l=1,llm
    DO ij=ijb,ije
      mw(ij,l)=w(ij,l) * zzw
    ENDDO
  ENDDO
  
  DO ij=ijb,ije
    mw(ij,llm+1)=0.
  ENDDO
  
  ijb=ij_begin
  ije=ij_end
  zq(ijb:ije,:)=q(ijb:ije,:)
  zm(ijb:ije,:)=masse(ijb:ije,:)
  
  zm1(ijb:ije,:) = zm(ijb:ije,:)
  zq1(ijb:ije,:) = zq(ijb:ije,:)
  
  CALL VTe(VTadvection_ad)
  
  CALL VTb(VTvlx)
  CALL vlx_p(zq,pente_max,zm,mu,ij_begin,ij_end)
  CALL VTe(VTvlx)
  
  
  zm2(ijb:ije,:) = zm(ijb:ije,:)
  zq2(ijb:ije,:) = zq(ijb:ije,:)
  
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest1)
  CALL Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest1)
  CALL Register_Hallo(mv,ip1jm,llm,1,1,1,1,MyRequest1)
  CALL SendRequest(MyRequest1)
  CALL WaitRequest(MyRequest1)
  CALL VTe(VTHallo)
  
  CALL VTb(VTvly)
  CALL vly_p(zq,pente_max,zm,mv)
  CALL VTe(VTvly)
  
  zm3(ijb:ije,:) = zm(ijb:ije,:)
  zq3(ijb:ije,:) = zq(ijb:ije,:)
  CALL VTb(VTvlz)
  CALL vlz_p(zq,pente_max,zm,mw,ij_begin,ij_end)
  CALL VTe(VTvlz)
  
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest2)
  CALL SendRequest(MyRequest2)
  CALL WaitRequest(MyRequest2)
  CALL VTe(VTHallo)
  
  zm4(ijb:ije,:) = zm(ijb:ije,:)
  zq4(ijb:ije,:) = zq(ijb:ije,:)
  CALL VTb(VTvly)
  CALL vly_p(zq,pente_max,zm,mv)
  CALL VTe(VTvly)
  
  zm5(ijb:ije,:) = zm(ijb:ije,:)
  zq5(ijb:ije,:) = zq(ijb:ije,:)
  CALL VTb(VTvlx)
  CALL vlx_p(zq,pente_max,zm,mu,ij_begin,ij_end)
  CALL VTe(VTvlx)
  
  
  ijb=ij_begin
  ije=ij_end
  
  DO l=1,llm
    DO ij=ijb,ije
      q(ij,l)=zq(ij,l)
    ENDDO
  ENDDO
  
  DO l=1,llm
    DO ij=ijb,ije-iip1+1,iip1
      q(ij+iim,l)=q(ij,l)
    ENDDO
  ENDDO
  
  ! start AD
  ijb=ij_begin
  ije=ij_end
  
  zq_ad(ijb:ije,:) = 0.
  DO l=1,llm
    DO ij=ijb,ije-iip1+1,iip1
      q_ad(ij,l) = q_ad(ij,l) + q_ad(ij+iim,l)
      q_ad(ij+iim,l) = 0.
    ENDDO
  ENDDO
  
  DO l=1,llm
    DO ij=ijb,ije
      zq_ad(ij,l)=q_ad(ij,l)
      q_ad(ij,l)=0.
    ENDDO
  ENDDO
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq4,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(zq2,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(zm4,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(zm2,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(mv,ip1jm,llm,3,2,2,3,MyRequest2)
  CALL SendRequest(MyRequest2)
  CALL WaitRequest(MyRequest2)
  CALL VTe(VTHallo)
  

  
  CALL VTb(VTvlx)
  CALL vlx_ad_p(zq5,pente_max,zm5,mu,zq_ad,ij_begin,ij_end)
  CALL VTe(VTvlx)
  
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq_ad,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL SendRequest(MyRequest2)
  CALL WaitRequest(MyRequest2)
  CALL VTe(VTHallo)
  
  CALL VTb(VTvly)
  CALL vly_ad_p(zq4,pente_max,zm4,mv,zq_ad)
  CALL VTe(VTvly)

  CALL VTb(VTvlz)
  CALL vlz_ad_p(zq3,pente_max,zm3,mw,zq_ad,ij_begin,ij_end)
  CALL VTe(VTvlz)
  
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq_ad,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL SendRequest(MyRequest2)
  CALL WaitRequest(MyRequest2)
  CALL VTe(VTHallo)
  
  CALL VTb(VTvly)
  CALL vly_ad_p(zq2,pente_max,zm2,mv,zq_ad)
  CALL VTe(VTvly)

  CALL VTb(VTvlx)
  CALL vlx_ad_p(zq1,pente_max,zm1,mu,zq_ad,ij_begin,ij_end)
  CALL VTe(VTvlx)
  
  q_ad(ijb:ije,:) = zq_ad(ijb:ije,:)
  
  RETURN
END SUBROUTINE vlsplt_ad_p

SUBROUTINE vlx_ad_p(q,pente_max,masse,u_m,q_ad,ijb_x,ije_x)
  
  !     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
  !
  !    ********************************************************************
  !     Shema  d'advection " pseudo amont " .
  !    ********************************************************************
  !     nq,iq,q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
  !
  !
  !   --------------------------------------------------------------------
  USE Parallel
  IMPLICIT NONE
  !
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "logic.h"
  INCLUDE "comvert.h"
  INCLUDE "comconst.h"
  !
  !
  !   Arguments:
  !   ----------
  REAL :: masse(ip1jmp1,llm),pente_max
  REAL :: u_m( ip1jmp1,llm )
  REAL :: q(ip1jmp1,llm)
  REAL :: q_ad(ip1jmp1,llm)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l,j,i,iju,ijq,indu(ip1jmp1,llm),niju(llm)
  INTEGER :: n0,iadvplus(ip1jmp1,llm),nl(llm)
  !
  REAL :: new_m(ip1jmp1,llm),zu_m,zdum(ip1jmp1,llm)
  REAL :: dxq(ip1jmp1,llm),dxqu(ip1jmp1,llm)
  REAL :: dxq1(ip1jmp1,llm)
  REAL :: adxqu(ip1jmp1,llm),dxqmax(ip1jmp1,llm)
  REAL :: u_mq(ip1jmp1,llm)
  
  REAL :: dxq_ad(ip1jmp1,llm),dxqu_ad(ip1jmp1)
  REAL :: adxqu_ad(ip1jmp1),dxqmax_ad(ip1jmp1,llm)
  REAL :: u_mq_ad(ip1jmp1,llm)
  REAL :: a_ad, a1(ip1jmp1,llm), a2(ip1jmp1,llm), a3(ip1jmp1,llm)
  INTEGER :: ijb,ije,ijb_x,ije_x
  !   calcul de la pente a droite et a gauche de la maille
  
  IF (pente_max<=-1.e-5) THEN
      PRINT *,'pas de AD'
      STOP
  ENDIF

  
  ijb=ijb_x
  ije=ije_x
  
  IF (north_pole.AND.ijb==1) ijb=ijb+iip1
  IF (south_pole.AND.ije==ip1jmp1)  ije=ije-iip1
  
  !   calcul des pentes avec limitation, Van Leer scheme I:
  !   -----------------------------------------------------
  
  !   calcul de la pente aux points u
  DO l = 1, llm
    DO ij=ijb,ije-1
      dxqu(ij,l)=q(ij+1,l)-q(ij,l)
    ENDDO
    DO ij=ijb+iip1-1,ije,iip1
      dxqu(ij,l)=dxqu(ij-iim,l)
    ENDDO
    
    DO ij=ijb,ije
      IF (dxqu(ij,l) > 0.) THEN
          adxqu(ij,l)=dxqu(ij,l)
      ELSE
          adxqu(ij,l)=-dxqu(ij,l)
      ENDIF
    ENDDO
    
    !   calcul de la pente maximum dans la maille en valeur absolue
    
    DO ij=ijb+1,ije
      IF (adxqu(ij-1,l) < adxqu(ij,l)) THEN
          dxqmax(ij,l)=pente_max*adxqu(ij-1,l)
      ELSE
          dxqmax(ij,l)=pente_max*adxqu(ij,l)
      ENDIF
      
    ENDDO
    
    DO ij=ijb+iip1-1,ije,iip1
      dxqmax(ij-iim,l)=dxqmax(ij,l)
    ENDDO
    
    DO ij=ijb+1,ije
      IF(dxqu(ij-1,l)*dxqu(ij,l)>0) THEN
          dxq(ij,l)=dxqu(ij-1,l)+dxqu(ij,l)
      ELSE
          !   extremum local
          dxq(ij,l)=0.
      ENDIF
      dxq(ij,l)=0.5*dxq(ij,l)
      a1(ij,l) = dxq(ij,l)
      a2(ij,l) = a1(ij,l)
      IF (a1(ij,l) < 0.) THEN
          a2(ij,l) = -a1(ij,l)
      ENDIF
      a3(ij,l) = a2(ij,l)
      IF (a2(ij,l) > dxqmax(ij,l)) THEN
          a3(ij,l) = dxqmax(ij,l)
      ENDIF
      dxq1(ij,l) = dxq(ij,l)
      IF ((dxq1(ij,l)*a3(ij,l)) < 0.) THEN
          dxq(ij,l) = -a3(ij,l)
      ELSE
          dxq(ij,l) = a3(ij,l)
      ENDIF
    ENDDO
    
  ENDDO ! l=1,llm
  
  
  !   bouclage de la pente en iip1:
  !   -----------------------------
  
  DO l=1,llm
    DO ij=ijb+iip1-1,ije,iip1
      dxq(ij-iim,l)=dxq(ij,l)
    ENDDO
    DO ij=ijb,ije
      iadvplus(ij,l)=0
    ENDDO
    
  ENDDO
  
  !   calcul des flux a gauche et a droite
  
  !   on cumule le flux correspondant a toutes les mailles dont la masse
  !   au travers de la paroi pENDant le pas de temps.
  
  DO l=1,llm
    DO ij=ijb,ije-1
      IF (u_m(ij,l)>0.) THEN
          zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l)
          u_mq(ij,l)=u_m(ij,l)*(q(ij,l)+0.5*zdum(ij,l)*dxq(ij,l))
      ELSE
          zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l)
          u_mq(ij,l)=u_m(ij,l)*(q(ij+1,l)-0.5*zdum(ij,l)*dxq(ij+1,l))
      ENDIF
    ENDDO
  ENDDO
  !   detection des points ou on advecte plus que la masse de la
  !   maille
  DO l=1,llm
    DO ij=ijb,ije-1
      IF(zdum(ij,l)<0) THEN
          iadvplus(ij,l)=1
          u_mq(ij,l)=0.
      ENDIF
    ENDDO
  ENDDO
  DO l=1,llm
    DO ij=ijb+iip1-1,ije,iip1
      iadvplus(ij,l)=iadvplus(ij-iim,l)
    ENDDO
  ENDDO
  
  
  !   traitement special pour le cas ou on advecte en longitude plus que le
  !   contenu de la maille.
  !   cette partie est mal vectorisee.
  
  !  calcul du nombre de maille sur lequel on advecte plus que la maille.
  
  n0=0
  DO l=1,llm
    nl(l)=0
    DO ij=ijb,ije
      nl(l)=nl(l)+iadvplus(ij,l)
    ENDDO
    n0=n0+nl(l)
  ENDDO
  
  DO l=1,llm
    IF(nl(l)>0) THEN
        iju=0
        !   indicage des mailles concernees par le traitement special
        DO ij=ijb,ije
          IF(iadvplus(ij,l)==1.AND.MOD(ij,iip1)/=0) THEN
              iju=iju+1
              indu(iju,l)=ij
          ENDIF
        ENDDO
        niju(l)=iju
        
        !  traitement des mailles
        DO iju=1,niju(l)
          ij=indu(iju,l)
          j=(ij-1)/iip1+1
          zu_m=u_m(ij,l)
          u_mq(ij,l)=0.
          IF(zu_m>0.) THEN
              ijq=ij
              i=ijq-(j-1)*iip1
              !   accumulation pour les mailles completements advectees
              DO WHILE(zu_m>masse(ijq,l))
                u_mq(ij,l)=u_mq(ij,l)+q(ijq,l)*masse(ijq,l)
                zu_m=zu_m-masse(ijq,l)
                i=MOD(i-2+iim,iim)+1
                ijq=(j-1)*iip1+i
              ENDDO
              !   ajout de la maille non completement advectee
              u_mq(ij,l)=u_mq(ij,l)+zu_m* &
                 (q(ijq,l)+0.5*(1.-zu_m/masse(ijq,l))*dxq(ijq,l))
          ELSE
              ijq=ij+1
              i=ijq-(j-1)*iip1
              !   accumulation pour les mailles completements advectees
              DO WHILE(-zu_m>masse(ijq,l))
                u_mq(ij,l)=u_mq(ij,l)-q(ijq,l)*masse(ijq,l)
                zu_m=zu_m+masse(ijq,l)
                i=MOD(i,iim)+1
                ijq=(j-1)*iip1+i
              ENDDO
              !   ajout de la maille non completement advectee
              u_mq(ij,l)=u_mq(ij,l)+zu_m*(q(ijq,l)- &
                 0.5*(1.+zu_m/masse(ijq,l))*dxq(ijq,l))
          ENDIF
        ENDDO
    ENDIF
  ENDDO
  
  !   bouclage en latitude
  DO l=1,llm
    DO ij=ijb+iip1-1,ije,iip1
      u_mq(ij,l)=u_mq(ij-iim,l)
    ENDDO
  ENDDO
  
  !   calcul des tENDances
  DO l=1,llm
    DO ij=ijb+1,ije
      new_m(ij,l)=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
      q(ij,l)=(q(ij,l)*masse(ij,l)+ &
         u_mq(ij-1,l)-u_mq(ij,l)) &
         /new_m(ij,l)
    ENDDO
    DO ij=ijb+iip1-1,ije,iip1
      q(ij-iim,l)=q(ij,l)
    ENDDO
  ENDDO
  
  !Start of AD
  u_mq_ad(ij_begin:ij_end,:) = 0.
  dxq_ad(ij_begin:ij_end,:) = 0.
  !   calcul des tENDances
  DO l=1,llm
    DO ij=ije,ijb+iip1-1,-iip1
      q_ad(ij,l) = q_ad(ij,l) + q_ad(ij-iim,l)
      q_ad(ij-iim,l) = 0.
    ENDDO
    DO ij=ije,ijb+1,-1
      u_mq_ad(ij-1,l) = u_mq_ad(ij-1,l)+q_ad(ij,l)/new_m(ij,l)
      u_mq_ad(ij,l) = u_mq_ad(ij,l)-q_ad(ij,l)/new_m(ij,l)
      q_ad(ij,l) = q_ad(ij,l)*masse(ij,l)/new_m(ij,l)
    ENDDO
  ENDDO
  
  !   bouclage en latitude
  DO l=1,llm
    DO ij=ije,ijb+iip1-1,-iip1
      u_mq_ad(ij-iim,l) = u_mq_ad(ij-iim,l) + u_mq_ad(ij,l)
      u_mq_ad(ij,l) = 0.
    ENDDO
  ENDDO
  
  DO l=llm,1,-1
    IF(nl(l)>0) THEN
            
        !  traitement des mailles
        DO iju=1,niju(l)
          ij=indu(iju,l)
          j=(ij-1)/iip1+1
          zu_m=u_m(ij,l)
          IF(zu_m > 0.) THEN
              ijq=ij
              i=ijq-(j-1)*iip1
              !   accumulation pour les mailles completements advectees
              DO WHILE(zu_m > masse(ijq,l))
                q_ad(ijq,l) = q_ad(ijq,l) + u_mq_ad(ij,l)*masse(ijq,l)
                zu_m=zu_m-masse(ijq,l)
                i=MOD(i-2+iim,iim)+1
                ijq=(j-1)*iip1+i
              ENDDO
              !   ajout de la maille non completement advectee
              ! AD permute avec ligne precedentes -> sans impact
              q_ad(ijq,l) = q_ad(ijq,l) + zu_m * u_mq_ad(ij,l)
              dxq_ad(ijq,l) = dxq_ad(ijq,l) &
                 + zu_m*0.5*(1.-zu_m/masse(ijq,l))*u_mq_ad(ij,l)
          ELSE
              ijq=ij+1
              i=ijq-(j-1)*iip1
              !   accumulation pour les mailles completements advectees
              DO WHILE(-zu_m>masse(ijq,l))
                q_ad(ijq,l) = q_ad(ijq,l) - u_mq_ad(ij,l)*masse(ijq,l)
                zu_m=zu_m+masse(ijq,l)
                i=MOD(i,iim)+1
                ijq=(j-1)*iip1+i
              ENDDO
              !   ajout de la maille non completement advectee
              ! AD permute avec ligne precedentes -> sans impact
              q_ad(ijq,l) = q_ad(ijq,l) + zu_m * u_mq_ad(ij,l)
              dxq_ad(ijq,l) = dxq_ad(ijq,l)  &
                 - zu_m*0.5*(1.+zu_m/masse(ijq,l))*u_mq_ad(ij,l)
          ENDIF
          u_mq_ad(ij,l) = 0.
        ENDDO
    ENDIF
  ENDDO
  
  !   detection des points ou on advecte plus que la masse de la
  !   maille
  DO l=1,llm
    DO ij=ijb,ije-1
      IF(zdum(ij,l)<0) u_mq_ad(ij,l)=0.
    ENDDO
  ENDDO
  
  !   on cumule le flux correspondant a toutes les mailles dont la masse
  !   au travers de la paroi pENDant le pas de temps.
  
  DO l=1,llm
    DO ij=ijb,ije-1
      IF (u_m(ij,l) > 0.) THEN
          q_ad(ij,l) = q_ad(ij,l) + u_m(ij,l)* u_mq_ad(ij,l)
          dxq_ad(ij,l) = dxq_ad(ij,l) + u_m(ij,l)*0.5*zdum(ij,l)*u_mq_ad(ij,l)
      ELSE
          q_ad(ij+1,l) = q_ad(ij+1,l) + u_m(ij,l)* u_mq_ad(ij,l)
          dxq_ad(ij+1,l) = dxq_ad(ij+1,l) &
             - 0.5*u_m(ij,l)*zdum(ij,l)*u_mq_ad(ij,l)
      ENDIF
      u_mq_ad(ij,l) = 0.
    ENDDO
  ENDDO
  
  !   bouclage de la pente en iip1:
  DO l=1,llm
    DO ij=ije,ijb+iip1-1,-iip1
      dxq_ad(ij,l) = dxq_ad(ij,l) + dxq_ad(ij-iim,l)
      dxq_ad(ij-iim,l) = 0.
    ENDDO
  ENDDO
  
  !   calcul de la pente aux points u
  dxqmax_ad(ij_begin:ij_end,:) = 0.
  DO l = 1, llm
    adxqu_ad(ij_begin:ij_end) = 0.
    dxqu_ad(ij_begin:ij_end) = 0.
    DO ij=ijb+1,ije
      a_ad = dxq_ad(ij,l) 
      dxq_ad(ij,l) = 0.
      IF ((dxq1(ij,l)*a3(ij,l)) < 0.) THEN
          a_ad = - a_ad
      ENDIF
      IF (a2(ij,l) > dxqmax(ij,l)) THEN
          dxqmax_ad(ij,l) = dxqmax_ad(ij,l) + a_ad
          a_ad = 0.
      ENDIF
      IF (a1(ij,l) < 0.) THEN
          a_ad = - a_ad
      ENDIF
      dxq_ad(ij,l) = dxq_ad(ij,l) + a_ad
      a_ad = 0.
      dxq_ad(ij,l) = 0.5*dxq_ad(ij,l)
      IF(dxqu(ij-1,l)*dxqu(ij,l)>0) THEN
          dxqu_ad(ij) = dxqu_ad(ij) + dxq_ad(ij,l)
          dxqu_ad(ij-1) = dxqu_ad(ij-1) + dxq_ad(ij,l)
          dxq_ad(ij,l)=0.
      ELSE
          ! extremum local
          dxq_ad(ij,l)=0.
      ENDIF
    ENDDO
    
    DO ij=ije,ijb+iip1-1,-iip1
      dxqmax_ad(ij,l) = dxqmax_ad(ij,l) + dxqmax_ad(ij-iim,l)
      dxqmax_ad(ij-iim,l) = 0.
    ENDDO
    
    !   calcul de la pente maximum dans la maille en valeur absolue
    DO ij=ijb+1,ije
      IF (adxqu(ij-1,l) < adxqu(ij,l)) THEN
          adxqu_ad(ij-1) = adxqu_ad(ij-1) + pente_max*dxqmax_ad(ij,l)
      ELSE
          adxqu_ad(ij) = adxqu_ad(ij) + pente_max*dxqmax_ad(ij,l)
      ENDIF
      dxqmax_ad(ij,l) = 0.
    ENDDO
    
    DO ij=ijb,ije
      IF (dxqu(ij,l) > 0.) THEN
          dxqu_ad(ij) = dxqu_ad(ij) + adxqu_ad(ij)
      ELSE
          dxqu_ad(ij) = dxqu_ad(ij) - adxqu_ad(ij)
      ENDIF
      adxqu_ad(ij) = 0.
    ENDDO
    
    DO ij=ije,ijb+iip1-1,-iip1
      dxqu_ad(ij-iim) = dxqu_ad(ij-iim) + dxqu_ad(ij)
      dxqu_ad(ij) = 0.
    ENDDO
    
    DO ij=ijb,ije-1
      q_ad(ij+1,l) = q_ad(ij+1,l) + dxqu_ad(ij)
      q_ad(ij,l) = q_ad(ij,l) - dxqu_ad(ij)
      dxqu_ad(ij) = 0.
    ENDDO
    
  ENDDO
  
  ! last bit of forward !
  masse(ij_begin:ij_end,:) = new_m(ij_begin:ij_end,:)
  
  RETURN
END SUBROUTINE vlx_ad_p

SUBROUTINE vly_ad_p(q,pente_max,masse,masse_adv_v,q_ad)
  !
  !     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
  !
  !    ********************************************************************
  !     Shema  d'advection " pseudo amont " .
  !    ********************************************************************
  !     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
  !     dq 	       sont des arguments de sortie pour le s-pg ....
  !
  !
  !   --------------------------------------------------------------------
  USE parallel
  USE WRITE_field_p
  IMPLICIT NONE
  !
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "logic.h"
  INCLUDE "comvert.h"
  INCLUDE "comconst.h"
  INCLUDE "comgeom.h"
  !
  !
  !   Arguments:
  !   ----------
  REAL :: masse(ip1jmp1,llm),pente_max
  REAL :: masse_adv_v( ip1jm,llm)
  REAL :: q(ip1jmp1,llm)
  REAL :: q_ad(ip1jmp1,llm)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: i,ij,l
  !
  REAL :: areaj2,areajjm,areascb(iim),areasch(iim)
  REAL :: areascb_ad(iim),areasch_ad(iim)
  REAL :: dyq(ip1jmp1,llm),dyqv(ip1jm,llm)
  REAL :: dyq_ad(ip1jmp1,llm),dyqv_ad(ip1jm)
  REAL :: adyqv(ip1jm,llm),dyqmax(ip1jmp1,llm)
  REAL :: adyqv_ad(ip1jm),dyqmax_ad(ip1jmp1)
  REAL :: qbyv(ip1jm,llm)
  REAL :: qbyv_ad(ip1jm,llm)
  
  REAL :: qpns,qpsn,dyn1,dys1,dyn2,dys2,newmasse(ip1jmp1,llm)
  REAL :: qpns_ad, qpsn_ad, dyn1_ad,dys1_ad, dyn2_ad, dys2_ad
  LOGICAL,SAVE :: first
  
  REAL :: convpn,convps,convmpn,convmps
  REAL :: convpn_ad,convps_ad
  REAL :: massepn,masseps,qpn,qps
  REAL :: qpn_ad,qps_ad
  REAL,SAVE :: sinlon(iip1),sinlondlon(iip1)
  REAL,SAVE :: coslon(iip1),coslondlon(iip1)
  SAVE :: areaj2,areajjm
  REAL :: a, a_ad, a1(ip1jmp1,llm), a2(ip1jmp1,llm), a3(ip1jmp1,llm)
  !
  !
  REAL ::     SSUM
  EXTERNAL :: SSUM
  
  DATA first/.TRUE./
  INTEGER ijb,ije
  
  
  IF(first) THEN
      !        PRINT*,'Shema  Amont nouveau  appele dans  Vanleer   '
      first=.FALSE.
      DO i=2,iip1
        coslon(i)=COS(rlonv(i))
        sinlon(i)=SIN(rlonv(i))
        coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
        sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
      ENDDO
      coslon(1)=coslon(iip1)
      coslondlon(1)=coslondlon(iip1)
      sinlon(1)=sinlon(iip1)
      sinlondlon(1)=sinlondlon(iip1)
      areaj2 = SSUM( iim, area(iip2), 1 )
      areajjm= SSUM( iim, area(ip1jm -iim), 1 )
  ENDIF
  
  !
  a1(:,:)=0.
  a2(:,:)=0.
  a3(:,:)=0.
  
  DO l = 1, llm
    !
    !   --------------------------------
    !      CALCUL EN LATITUDE
    !   --------------------------------
    
    !   On commence par calculer la valeur du traceur moyenne sur le premier cercle
    !   de latitude autour du pole (qpns pour le pole nord et qpsn pour
    !    le pole nord) qui sera utilisee pour evaluer les pentes au pole.
    
    IF (north_pole) THEN
        DO i = 1, iim
          areascb(i) = area(i+ iip1) * q(i+ iip1,l)
        ENDDO
        qpns   = SSUM( iim,  areascb ,1 ) / areaj2
    ENDIF
    
    IF (south_pole) THEN
        DO i = 1, iim
          areasch(i) = area(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
        ENDDO
        qpsn   = SSUM( iim,  areasch ,1 ) / areajjm
    ENDIF
    
    
    !   calcul des pentes aux points v
    ijb=ij_begin-2*iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-iip1
    dyqv(ijb:ije,l)=0.
    
    DO ij=ijb,ije
      dyqv(ij,l)=q(ij,l)-q(ij+iip1,l)
      IF (dyqv(ij,l) > 0.) THEN
          adyqv(ij,l)=dyqv(ij,l)
      ELSE
          adyqv(ij,l)=-dyqv(ij,l)
      ENDIF
    ENDDO
    
    !   calcul des pentes aux points scalaires
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyq(ij,l)=.5*(dyqv(ij-iip1,l)+dyqv(ij,l))
      IF (adyqv(ij-iip1,l) > adyqv(ij,l)) THEN
          dyqmax(ij,l)=adyqv(ij,l)
      ELSE
          dyqmax(ij,l)=adyqv(ij-iip1,l)
      ENDIF
      dyqmax(ij,l)=pente_max*dyqmax(ij,l)
    ENDDO
    
    IF (north_pole) THEN
        DO ij=1,iip1
          dyq(ij,l)=qpns-q(ij+iip1,l)
        ENDDO
        
        dyn1=0.
        dyn2=0.
        DO ij=1,iim
          dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
          dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
        ENDDO
        DO ij=1,iip1
          dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
          dyq(ij,l)=0.
        ENDDO
    ENDIF
    
    IF (south_pole) THEN
        
        DO ij=1,iip1
          dyq(ip1jm+ij,l)=q(ip1jm+ij-iip1,l)-qpsn
        ENDDO
        
        dys1=0.
        dys2=0.
        
        DO ij=1,iim
          dys1=dys1+sinlondlon(ij)*dyq(ip1jm+ij,l)
          dys2=dys2+coslondlon(ij)*dyq(ip1jm+ij,l)
        ENDDO
        
        DO ij=1,iip1
          dyq(ip1jm+ij,l)=dys1*sinlon(ij)+dys2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
          dyq(ip1jm+ij,l)=0.
        ENDDO
    ENDIF
    
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end
    
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    !   calcul des pentes limitees
    
    DO ij=ijb,ije
      IF(dyqv(ij,l)*dyqv(ij-iip1,l)>0.) THEN
          a = dyq(ij,l)
          a1(ij,l) = a
          IF (a1(ij,l) < 0.) THEN
              a = -a
          ENDIF
          a2(ij,l) = a
          IF (a2(ij,l) > dyqmax(ij,l)) THEN
              a = dyqmax(ij,l)
          ENDIF
          a3(ij,l) = a
          IF ((dyq(ij,l)*a3(ij,l)) < 0.) THEN
              a = -a
          ENDIF
          dyq(ij,l)=a
      ELSE
          dyq(ij,l)=0.
      ENDIF
    ENDDO
    
  ENDDO
  
  
  ijb=ij_begin-iip1
  ije=ij_end
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  DO l=1,llm
    DO ij=ijb,ije
      IF(masse_adv_v(ij,l)>0) THEN
          qbyv(ij,l)=q(ij+iip1,l)+dyq(ij+iip1,l)* &
             0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
      ELSE
          qbyv(ij,l)=q(ij,l)-dyq(ij,l)* &
             0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
      ENDIF
      qbyv(ij,l)=masse_adv_v(ij,l)*qbyv(ij,l)
    ENDDO
  ENDDO
  
  
  ijb=ij_begin-2*iip1
  ije=ij_end+2*iip1
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end
  
  newmasse(ijb:ije,:) = masse(ijb:ije,:)
  
  ijb=ij_begin-2*iip1
  ije=ij_end+iip1
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  qbyv_ad(ijb:ije,:)=0.
  
  
  DO l=1,llm
    
    ijb=ij_begin-2*iip1
    ije=ij_end+2*iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      newmasse(ij,l)=masse(ij,l) &
         +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)
    ENDDO
    
    
    ijb=ij_begin
    ije=ij_end
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      q(ij,l)=(q(ij,l)*masse(ij,l)+qbyv(ij,l)-qbyv(ij-iip1,l)) &
         /newmasse(ij,l)
    ENDDO
    
    
    IF (north_pole) THEN
        convpn=SSUM(iim,qbyv(1,l),1)
        convmpn=ssum(iim,masse_adv_v(1,l),1)
        massepn=ssum(iim,masse(1,l),1)
        qpn=0.
        DO ij=1,iim
          qpn=qpn+masse(ij,l)*q(ij,l)
        ENDDO
        qpn=(qpn+convpn)/(massepn+convmpn)
        DO ij=1,iip1
          q(ij,l)=qpn
        ENDDO
        
        !Start of AD
        
        qpn_ad=0.
        DO ij=1,iip1
          qpn_ad = qpn_ad +q_ad(ij,l)
          q_ad(ij,l) = 0.
        ENDDO
        convpn_ad = qpn_ad /(massepn+convmpn)
        qpn_ad = qpn_ad /(massepn+convmpn)
        DO ij=1,iim
          q_ad(ij,l) = q_ad(ij,l) + newmasse(ij,l)*qpn_ad
        ENDDO
        qpn_ad=0.
        
        DO ij= 1, iim
          qbyv_ad(ij,l) = qbyv_ad(ij,l) + convpn_ad
        END DO
        convpn_ad = 0.
    ENDIF
    
    IF (south_pole) THEN
        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
        masseps=ssum(iim, masse(ip1jm+1,l),1)
        qps=0.
        DO ij = ip1jm+1,ip1jmp1-1
          qps=qps+masse(ij,l)*q(ij,l)
        ENDDO
        qps=(qps+convps)/(masseps+convmps)
        DO ij=ip1jm+1,ip1jmp1
          q(ij,l)=qps
        ENDDO
        
        !Start of AD
        qps_ad = 0.
        DO ij=ip1jm+1,ip1jmp1
          qps_ad = qps_ad +q_ad(ij,l)
          q_ad(ij,l) = 0.
        ENDDO
        convps_ad = qps_ad /(masseps+convmps)
        qps_ad = qps_ad /(masseps+convmps)
        DO ij = ip1jm+1,ip1jmp1-1
          q_ad(ij,l) = q_ad(ij,l) + newmasse(ij,l)*qps_ad
        ENDDO
        qps_ad = 0.
        
        DO ij = ip1jm-iim, ip1jm-1
          qbyv_ad(ij,l) = qbyv_ad(ij,l) - convps_ad
        ENDDO
        convps_ad = 0.
    ENDIF
    
    
    
    ijb=ij_begin-2*iip1
    ije=ij_end+iip1
    
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      qbyv_ad(ij,l) = qbyv_ad(ij,l) + q_ad(ij,l)/newmasse(ij,l)
    ENDDO
    
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-2*iip1
    
    DO ij=ijb,ije
      qbyv_ad(ij,l) = qbyv_ad(ij,l) - q_ad(ij+iip1,l)/newmasse(ij+iip1,l)
    ENDDO
    
    
    ijb=ij_begin
    ije=ij_end
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      q_ad(ij,l)=q_ad(ij,l)*masse(ij,l)/newmasse(ij,l)
    ENDDO
    
    
  END DO
  
  
  ijb=ij_begin-2*iip1
  ije=ij_end+iip1
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  DO l=1,llm
    DO ij=ijb,ije
      qbyv_ad(ij,l)=masse_adv_v(ij,l)*qbyv_ad(ij,l)
    ENDDO
  ENDDO
  
  
  ijb=ij_begin
  ije=ij_end
  
  IF (north_pole) ijb=ij_begin+iip1
  IF (south_pole)  ije=ij_end
  
  DO l=1,llm
    DO ij=ijb,ije
      IF(masse_adv_v(ij-iip1,l) > 0.) THEN
          q_ad(ij,l) = q_ad(ij,l) + qbyv_ad(ij-iip1,l)
      ENDIF
    ENDDO
  ENDDO
  
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  DO l=1,llm
    DO ij=ijb,ije
      IF(.NOT.(masse_adv_v(ij,l) > 0.)) THEN
          q_ad(ij,l) = q_ad(ij,l) + qbyv_ad(ij,l)
      ENDIF
    ENDDO
  ENDDO
  
  
  ijb=ij_begin-iip1
  ije=ij_end+iip1
  
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end
  dyq_ad(ijb:ije,:) = 0.
  
  IF (north_pole) ijb=ij_begin+iip1
  IF (south_pole)  ije=ij_end
  
  DO l=1,llm
    DO ij=ijb,ije
      IF(masse_adv_v(ij-iip1,l) > 0.) THEN 
          dyq_ad(ij,l) = dyq_ad(ij,l) + qbyv_ad(ij-iip1,l)* &
             0.5*(1.-masse_adv_v(ij-iip1,l)/masse(ij,l))
      ENDIF
    ENDDO
  ENDDO
  
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  DO l=1,llm
    DO ij=ijb,ije
      IF(.NOT.(masse_adv_v(ij,l) > 0.)) THEN
          dyq_ad(ij,l) = dyq_ad(ij,l) - qbyv_ad(ij,l)* &
             0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
      ENDIF
    ENDDO
  ENDDO
  
  
  
  DO l = 1, llm
    
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    dyqmax_ad(ijb:ije) = 0.
    
    DO ij=ijb,ije
      IF(dyqv(ij,l)*dyqv(ij-iip1,l).GT.0.) THEN
          a_ad = dyq_ad(ij,l)
          dyq_ad(ij,l) = 0.
          IF ((dyq(ij,l)*a3(ij,l)) < 0.) THEN
              a_ad = -a_ad
          ENDIF
          IF (a2(ij,l) > dyqmax(ij,l)) THEN
              dyqmax_ad(ij) = dyqmax_ad(ij) + a_ad
              a_ad = 0.
          ENDIF
          IF (a1(ij,l) < 0.) THEN
              a_ad = -a_ad
          ENDIF
          dyq_ad(ij,l) = dyq_ad(ij,l) + a_ad
      ELSE
          dyq_ad(ij,l)=0.
          
      ENDIF
    ENDDO
    
    IF (north_pole) THEN
        
        DO ij=1,iip1
          dyq_ad(ij,l)=0.
        ENDDO
        
        dyn1_ad = 0.
        dyn2_ad = 0.
        DO ij=1,iip1
          dyn1_ad = dyn1_ad + dyq_ad(ij,l)*sinlon(ij)
          dyn2_ad = dyn2_ad + dyq_ad(ij,l)*coslon(ij)
        ENDDO
        
        DO ij=1,iim
          dyq_ad(ij,l) = dyq_ad(ij,l) + dyn1_ad*sinlondlon(ij)
          dyq_ad(ij,l) = dyq_ad(ij,l) + dyn2_ad*coslondlon(ij)
        ENDDO
        
        dyn1_ad = 0.
        dyn2_ad = 0.
        qpns_ad = 0.
        DO ij=1,iip1
          qpns_ad = qpns_ad + dyq_ad(ij,l)
          q_ad(ij+iip1,l) = q_ad(ij+iip1,l) - dyq_ad(ij,l)
        ENDDO
        
    ENDIF
    
    
    IF (south_pole) THEN
        
        DO ij=1,iip1
          dyq_ad(ip1jm+ij,l)=0.
        ENDDO
        
        dys1_ad = 0.
        dys2_ad = 0.
        DO ij=1,iip1
          dys1_ad = dys1_ad + dyq_ad(ip1jm+ij,l)*sinlon(ij)
          dys2_ad = dys2_ad + dyq_ad(ip1jm+ij,l)*coslon(ij)
        ENDDO
        
        DO ij=1,iim
          dyq_ad(ip1jm+ij,l) = dyq_ad(ip1jm+ij,l) + dys1_ad*sinlondlon(ij)
          dyq_ad(ip1jm+ij,l) = dyq_ad(ip1jm+ij,l) + dys2_ad*coslondlon(ij)
        ENDDO
        
        dys1_ad = 0.
        dys2_ad = 0.
        qpsn_ad = 0.
        
        DO ij=1,iip1
          q_ad(ip1jm+ij-iip1,l) = q_ad(ip1jm+ij-iip1,l) + dyq_ad(ip1jm+ij,l)
          qpsn_ad = qpsn_ad - dyq_ad(ip1jm+ij,l)
        ENDDO
        
    ENDIF
    
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole) ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyqmax_ad(ij) = pente_max*dyqmax_ad(ij)
    ENDDO
    
    
    ijb=ij_begin-iip1
    ije=ij_end
    
    IF (north_pole) ijb=ij_begin
    IF (south_pole) ije=ij_end-iip1
    
    adyqv_ad(ijb:ije) = 0.
    dyqv_ad(ijb:ije) = 0.
    
    
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      IF (adyqv(ij-iip1,l) > adyqv(ij,l)) THEN
          adyqv_ad(ij) = adyqv_ad(ij) + dyqmax_ad(ij)
      ENDIF
      dyqv_ad(ij) = dyqv_ad(ij) + .5*dyq_ad(ij,l)
    ENDDO
    
    
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-2*iip1
    
    DO ij=ijb,ije
      IF (.NOT.(adyqv(ij,l) > adyqv(ij+iip1,l)) ) THEN
          adyqv_ad(ij) = adyqv_ad(ij) + dyqmax_ad(ij+iip1)
      ENDIF
      dyqv_ad(ij) = dyqv_ad(ij) + .5*dyq_ad(ij+iip1,l)
      dyq_ad(ij+iip1,l)=0.
    ENDDO
    
    ijb=ij_begin-iip1
    ije=ij_end
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      IF (dyqv(ij,l) > 0.) THEN
          dyqv_ad(ij) = dyqv_ad(ij) + adyqv_ad(ij)
      ELSE
          dyqv_ad(ij) = dyqv_ad(ij) - adyqv_ad(ij)
      ENDIF
    ENDDO
    
    
    
    ijb=ij_begin
    ije=ij_end
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end
    
    DO ij=ijb,ije
      q_ad(ij,l) = q_ad(ij,l) - dyqv_ad(ij-iip1)
    ENDDO
    
    
    ijb=ij_begin
    ije=ij_end
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-iip1
        
    DO ij=ijb,ije
      q_ad(ij,l) = q_ad(ij,l) + dyqv_ad(ij)
    ENDDO
    
    
    
    IF (north_pole) THEN
        
        areascb_ad(:) = 0.
        qpns_ad = qpns_ad/areaj2
        DO i = 1, iim
          areascb_ad(i) = areascb_ad(i) + qpns_ad
        ENDDO
        qpns_ad = 0.
        
        DO i = 1, iim
          q_ad(i+ iip1,l) = q_ad(i+ iip1,l) + area(i+ iip1) * areascb_ad(i)
        ENDDO
        
    ENDIF
    
    
    IF (south_pole) THEN
        
        areasch_ad(:) = 0.
        qpsn_ad = qpsn_ad/areajjm
        DO i = 1, iim
          areasch_ad(i) = areasch_ad(i) + qpsn_ad
        ENDDO
        qpsn_ad = 0.
        
        DO i = 1, iim
          q_ad(i+ip1jm-iip1,l) = q_ad(i+ip1jm-iip1,l) &
             + area(i+ip1jm-iip1) * areasch_ad(i)
        ENDDO
        
    ENDIF
    
  ENDDO
  
  
  !last bit of forward
  ijb=ij_begin
  ije=ij_end
  IF (north_pole) ijb=ij_begin+iip1
  IF (south_pole)  ije=ij_end-iip1
  
  DO l=1,llm
    DO ij=ijb,ije
      masse(ij,l) = newmasse(ij,l)
    ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE vly_ad_p
      


SUBROUTINE vly_ad_p_bis(q,pente_max,masse,masse_adv_v,q_ad)
  !
  !     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
  !
  !    ********************************************************************
  !     Shema  d'advection " pseudo amont " .
  !    ********************************************************************
  !     q,masse_adv_v,w sont des arguments d'entree  pour le s-pg ....
  !     dq 	       sont des arguments de sortie pour le s-pg ....
  !
  !
  !   --------------------------------------------------------------------
  USE parallel
  IMPLICIT NONE
  !
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "logic.h"
  INCLUDE "comvert.h"
  INCLUDE "comconst.h"
  INCLUDE "comgeom.h"
  !
  !
  !   Arguments:
  !   ----------
  REAL :: masse(ip1jmp1,llm),pente_max
  REAL :: masse_adv_v( ip1jm,llm)
  REAL :: q(ip1jmp1,llm)
  REAL :: q_ad(ip1jmp1,llm)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: i,ij,l
  !
  REAL :: areaj2,areajjm,areascb(iim),areasch(iim)
  REAL :: areascb_ad(iim),areasch_ad(iim)
  REAL :: dyq(ip1jmp1,llm),dyqv(ip1jm,llm)
  REAL :: dyq_ad(ip1jmp1,llm),dyqv_ad(ip1jm)
  REAL :: adyqv(ip1jm,llm),dyqmax(ip1jmp1,llm)
  REAL :: adyqv_ad(ip1jm),dyqmax_ad(ip1jmp1)
  REAL :: qbyv(ip1jm,llm)
  REAL :: qbyv_ad(ip1jm,llm)
  
  REAL :: qpns,qpsn,dyn1,dys1,dyn2,dys2,newmasse(ip1jmp1,llm)
  REAL :: qpns_ad, qpsn_ad, dyn1_ad,dys1_ad, dyn2_ad, dys2_ad
  LOGICAL,SAVE ::first
  
  REAL :: convpn,convps,convmpn,convmps
  REAL :: convpn_ad,convps_ad
  REAL :: massepn,masseps,qpn,qps
  REAL :: qpn_ad,qps_ad
  REAL,SAVE :: sinlon(iip1),sinlondlon(iip1)
  REAL,SAVE :: coslon(iip1),coslondlon(iip1)
  SAVE :: areaj2,areajjm
  REAL :: a, a_ad, a1(ip1jmp1,llm), a2(ip1jmp1,llm), a3(ip1jmp1,llm)
  !
  !
  REAL ::     SSUM
  EXTERNAL :: SSUM
  
  DATA first/.TRUE./
  INTEGER :: ijb,ije
  
  
  IF(first) THEN
      !        PRINT*,'Shema  Amont nouveau  appele dans  Vanleer   '
      first=.FALSE.
      DO i=2,iip1
        coslon(i)=COS(rlonv(i))
        sinlon(i)=SIN(rlonv(i))
        coslondlon(i)=coslon(i)*(rlonu(i)-rlonu(i-1))/pi
        sinlondlon(i)=sinlon(i)*(rlonu(i)-rlonu(i-1))/pi
      ENDDO
      coslon(1)=coslon(iip1)
      coslondlon(1)=coslondlon(iip1)
      sinlon(1)=sinlon(iip1)
      sinlondlon(1)=sinlondlon(iip1)
      areaj2 = SSUM( iim, area(iip2), 1 )
      areajjm= SSUM( iim, area(ip1jm -iim), 1 )
  ENDIF
  
  DO l = 1, llm
    !
    !   --------------------------------
    !      CALCUL EN LATITUDE
    !   --------------------------------
    
    !   On commence par calculer la valeur du traceur moyenne sur le premier cercle
    !   de latitude autour du pole (qpns pour le pole nord et qpsn pour
    !    le pole nord) qui sera utilisee pour evaluer les pentes au pole.
    
    IF (north_pole) THEN
        DO i = 1, iim
          areascb(i) = area(i+ iip1) * q(i+ iip1,l)
        ENDDO
        qpns   = SSUM( iim,  areascb ,1 ) / areaj2
    ENDIF
    
    IF (south_pole) THEN
        DO i = 1, iim
          areasch(i) = area(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
        ENDDO
        qpsn   = SSUM( iim,  areasch ,1 ) / areajjm
    ENDIF
    
    
    !   calcul des pentes aux points v
    ijb=ij_begin-2*iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyqv(ij,l)=q(ij,l)-q(ij+iip1,l)
      IF (dyqv(ij,l) > 0.) THEN
          adyqv(ij,l)=dyqv(ij,l)
      ELSE
          adyqv(ij,l)=-dyqv(ij,l)
      ENDIF
    ENDDO
    
    !   calcul des pentes aux points scalaires
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyq(ij,l)=.5*(dyqv(ij-iip1,l)+dyqv(ij,l))
      IF (adyqv(ij-iip1,l) > adyqv(ij,l)) THEN
          dyqmax(ij,l)=adyqv(ij,l)
      ELSE
          dyqmax(ij,l)=adyqv(ij-iip1,l)
      ENDIF
      dyqmax(ij,l)=pente_max*dyqmax(ij,l)
    ENDDO
    
    IF (north_pole) THEN
        DO ij=1,iip1
          dyq(ij,l)=qpns-q(ij+iip1,l)
        ENDDO
        
        dyn1=0.
        dyn2=0.
        DO ij=1,iim
          dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
          dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
        ENDDO
        DO ij=1,iip1
          dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
          dyq(ij,l)=0.
        ENDDO
    ENDIF
    
    IF (south_pole) THEN
        
        DO ij=1,iip1
          dyq(ip1jm+ij,l)=q(ip1jm+ij-iip1,l)-qpsn
        ENDDO
        
        dys1=0.
        dys2=0.
        
        DO ij=1,iim
          dys1=dys1+sinlondlon(ij)*dyq(ip1jm+ij,l)
          dys2=dys2+coslondlon(ij)*dyq(ip1jm+ij,l)
        ENDDO
        
        DO ij=1,iip1
          dyq(ip1jm+ij,l)=dys1*sinlon(ij)+dys2*coslon(ij)
        ENDDO
        
        DO ij=1,iip1
          dyq(ip1jm+ij,l)=0.
        ENDDO
    ENDIF
    
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    !   calcul des pentes limitees
    
    DO ij=ijb,ije
      IF(dyqv(ij,l)*dyqv(ij-iip1,l)>0.) THEN
          a = dyq(ij,l)
          a1(ij,l) = a
          IF (a1(ij,l) < 0.) THEN
              a = -a
          ENDIF
          a2(ij,l) = a
          IF (a2(ij,l) > dyqmax(ij,l)) THEN
              a = dyqmax(ij,l)
          ENDIF
          a3(ij,l) = a
          IF ((dyq(ij,l)*a3(ij,l)) < 0.) THEN
              a = -a
          ENDIF
          dyq(ij,l)=a
      ELSE
          dyq(ij,l)=0.
      ENDIF
    ENDDO
    
  ENDDO
  
  
  ijb=ij_begin-iip1
  ije=ij_end
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  DO l=1,llm
    DO ij=ijb,ije
      IF(masse_adv_v(ij,l)>0) THEN
          qbyv(ij,l)=q(ij+iip1,l)+dyq(ij+iip1,l)* &
             0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
      ELSE
          qbyv(ij,l)=q(ij,l)-dyq(ij,l)* &
             0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
      ENDIF
      qbyv(ij,l)=masse_adv_v(ij,l)*qbyv(ij,l)
    ENDDO
  ENDDO
  
  ijb=ij_begin-iip1
  ije=ij_end+2*iip1
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end
  
  newmasse(ijb:ije,:) = masse(ijb:ije,:)
  
  DO l=1,llm
    
    ijb=ij_begin-iip1
    ije=ij_end+2*iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      newmasse(ij,l)=masse(ij,l) &
         +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)
    ENDDO
    
    ijb=ij_begin
    ije=ij_end
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      q(ij,l)=(q(ij,l)*masse(ij,l)+qbyv(ij,l)-qbyv(ij-iip1,l)) &
         /newmasse(ij,l)
    ENDDO
    
    IF (north_pole) THEN
        convpn=SSUM(iim,qbyv(1,l),1)
        convmpn=ssum(iim,masse_adv_v(1,l),1)
        massepn=ssum(iim,masse(1,l),1)
        qpn=0.
        DO ij=1,iim
          qpn=qpn+masse(ij,l)*q(ij,l)
        ENDDO
        qpn=(qpn+convpn)/(massepn+convmpn)
        DO ij=1,iip1
          q(ij,l)=qpn
        ENDDO
    ENDIF
    
    
    IF (south_pole) THEN
        
        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
        masseps=ssum(iim, masse(ip1jm+1,l),1)
        qps=0.
        DO ij = ip1jm+1,ip1jmp1-1
          qps=qps+masse(ij,l)*q(ij,l)
        ENDDO
        qps=(qps+convps)/(masseps+convmps)
        DO ij=ip1jm+1,ip1jmp1
          q(ij,l)=qps
        ENDDO
    ENDIF
    
    !Start of AD
    ijb=ij_begin-2*iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-iip1
    
    qbyv_ad(ijb:ije,l) = 0.
    
    IF (south_pole) THEN
        qps_ad = 0.
        DO ij=ip1jm+1,ip1jmp1
          qps_ad = qps_ad +q_ad(ij,l)
          q_ad(ij,l) = 0.
        ENDDO
        convps_ad = qps_ad /(masseps+convmps)
        qps_ad = qps_ad /(masseps+convmps)
        DO ij = ip1jm+1,ip1jmp1-1
          q_ad(ij,l) = q_ad(ij,l) + newmasse(ij,l)*qps_ad
        ENDDO
        qps_ad = 0.
        
        DO ij = ip1jm-iim, ip1jm-1
          qbyv_ad(ij,l) = qbyv_ad(ij,l) - convps_ad
        ENDDO
        convps_ad = 0.
    ENDIF
    
    IF (north_pole) THEN
        qpn_ad=0.
        DO ij=1,iip1
          qpn_ad = qpn_ad +q_ad(ij,l)
          q_ad(ij,l) = 0.
        ENDDO
        convpn_ad = qpn_ad /(massepn+convmpn)
        qpn_ad = qpn_ad /(massepn+convmpn)
        DO ij=1,iim
          q_ad(ij,l) = q_ad(ij,l) + newmasse(ij,l)*qpn_ad
        ENDDO
        qpn_ad=0.
        
        DO ij= 1, iim
          qbyv_ad(ij,l) = qbyv_ad(ij,l) + convpn_ad
        ENDDO
        convpn_ad = 0.
        
    ENDIF
    
    ijb=ij_begin-2*iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      qbyv_ad(ij,l) = qbyv_ad(ij,l) + q_ad(ij,l)/newmasse(ij,l)
    ENDDO
    
    
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-2*iip1
    DO ij=ijb,ije    
      qbyv_ad(ij,l) = qbyv_ad(ij,l) - q_ad(ij+iip1,l)/newmasse(ij+iip1,l)
    ENDDO
    
    ijb=ij_begin-iip1
    ije=ij_end+2*iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      q_ad(ij,l)=q_ad(ij,l)*masse(ij,l)/newmasse(ij,l)
    ENDDO
    
    
  ENDDO
  
  ijb=ij_begin-2*iip1
  ije=ij_end+iip1
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  dyq_ad(ijb:ije,:) = 0.
  
  DO l=1,llm
    DO ij=ijb,ije
      qbyv_ad(ij,l)=masse_adv_v(ij,l)*qbyv_ad(ij,l)
      IF(masse_adv_v(ij,l) > 0.) THEN
          q_ad(ij+iip1,l) = q_ad(ij+iip1,l) + qbyv_ad(ij,l)
          dyq_ad(ij+iip1,l) = dyq_ad(ij+iip1,l) + qbyv_ad(ij,l)* &
             0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
      ELSE
          q_ad(ij,l) = q_ad(ij,l) + qbyv_ad(ij,l)
          dyq_ad(ij,l) = dyq_ad(ij,l) - qbyv_ad(ij,l)* &
             0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
      ENDIF
    ENDDO
  ENDDO
  
  ijb=ij_begin-iip1
  ije=ij_end
  IF (north_pole) ijb=ij_begin+iip1
  IF (south_pole)  ije=ij_end-iip1
  
  
  DO l = 1, llm
    
    DO ij=ijb,ije
      dyqmax_ad(ij) = 0.
      IF(dyqv(ij,l)*dyqv(ij-iip1,l)>0.) THEN
          a_ad = dyq_ad(ij,l)
          dyq_ad(ij,l) = 0.
          IF ((dyq(ij,l)*a3(ij,l)) < 0.) THEN
              a_ad = -a_ad
          ENDIF
          IF (a2(ij,l) > dyqmax(ij,l)) THEN
              dyqmax_ad(ij) = dyqmax_ad(ij) + a_ad
              a_ad = 0.
          ENDIF
          IF (a1(ij,l) < 0.) THEN
              a_ad = -a_ad
          ENDIF
          dyq_ad(ij,l) = dyq_ad(ij,l) + a_ad
      ELSE
          dyq_ad(ij,l)=0.
      ENDIF
    ENDDO
    
    
    IF (north_pole) THEN
        
        DO ij=1,iip1
          dyq_ad(ij,l)=0.
        ENDDO
        
        dyn1_ad = 0.
        dyn2_ad = 0.
        DO ij=1,iip1
          dyn1_ad = dyn1_ad + dyq_ad(ij,l)*sinlon(ij)
          dyn2_ad = dyn2_ad + dyq_ad(ij,l)*coslon(ij)
        ENDDO
        DO ij=1,iim
          dyq_ad(ij,l) = dyq_ad(ij,l) + dyn1_ad*sinlondlon(ij)
          dyq_ad(ij,l) = dyq_ad(ij,l) + dyn2_ad*coslondlon(ij)
        ENDDO
        
        dyn1_ad = 0.
        dyn2_ad = 0.
        qpns_ad = 0.
        
        DO ij=1,iip1
          qpns_ad = qpns_ad + dyq_ad(ij,l)
          q_ad(ij+iip1,l) = q_ad(ij+iip1,l) - dyq_ad(ij,l)
          dyq_ad(ij,l) = 0.
        ENDDO
        
    ENDIF
    
    
    IF (south_pole) THEN
        
        DO ij=1,iip1
          dyq_ad(ij,l)=0.
        ENDDO
        
        dys1_ad = 0.
        dys2_ad = 0.
        DO ij=1,iip1
          dys1_ad = dys1_ad + dyq_ad(ip1jm+ij,l)*sinlon(ij)
          dys2_ad = dys2_ad + dyq_ad(ip1jm+ij,l)*coslon(ij)
        ENDDO
        
        DO ij=1,iim
          dyq_ad(ip1jm+ij,l) = dyq_ad(ip1jm+ij,l) + dys1_ad*sinlondlon(ij)
          dyq_ad(ip1jm+ij,l) = dyq_ad(ip1jm+ij,l) + dys2_ad*coslondlon(ij)
        ENDDO
        dys1_ad = 0.
        dys2_ad = 0.
        
        qpsn_ad = 0.
        DO ij=1,iip1
          q_ad(ip1jm+ij-iip1,l) = q_ad(ip1jm+ij-iip1,l) + dyq_ad(ip1jm+ij,l)
          qpsn_ad = qpsn_ad - dyq_ad(ip1jm+ij,l)
          dyq_ad(ip1jm+ij,l) = 0.
        ENDDO
        
    ENDIF
    
    
    
    ijb=ij_begin-2*iip1
    ije=ij_end
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-iip1
    adyqv_ad(ijb:ije) = 0.
    
    ijb=ij_begin-iip1
    ije=ij_end
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyqmax_ad(ij) = pente_max*dyqmax_ad(ij)
      IF (adyqv(ij-iip1,l) > adyqv(ij,l)) THEN
          adyqv_ad(ij) = adyqv_ad(ij) + dyqmax_ad(ij)
      ELSE
          adyqv_ad(ij-iip1) = adyqv_ad(ij-iip1) + dyqmax_ad(ij)
      ENDIF
    ENDDO
    
    
    ijb=ij_begin
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end
    
    DO ij=ijb,ije 
      dyqv_ad(ij-iip1) = .5*dyq_ad(ij,l)
    ENDDO
    
    ijb=ij_begin-iip1
    ije=ij_end
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole) ije=ije-iip1
    
    DO ij=ijb,ije   
      dyqv_ad(ij) = dyqv_ad(ij) + .5*dyq_ad(ij,l)
    ENDDO
    
    ijb=ij_begin-iip1
    ije=ij_end
    IF (north_pole) ijb=ij_begin
    
    DO ij=ijb,ije
      IF (dyqv(ij,l) > 0.) THEN
          dyqv_ad(ij) = dyqv_ad(ij) + adyqv_ad(ij)
          adyqv_ad(ij)=0.
      ELSE
          dyqv_ad(ij) = dyqv_ad(ij) - adyqv_ad(ij)
          adyqv_ad(ij)=0.
      ENDIF
    ENDDO
    
    
    ijb=ij_begin
    ije=ij_end
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ije-iip1
    DO ij=ijb,ije 
      q_ad(ij,l) = (q_ad(ij,l) + dyqv_ad(ij)) -dyqv_ad(ij-iip1)
    ENDDO
    
    IF (north_pole) THEN
        DO ij=1,iip1
          q_ad(ij,l) = q_ad(ij,l) + dyqv_ad(ij)
        ENDDO
    ENDIF
    
    IF (south_pole) THEN
        DO ij=ip1jm+1,ip1jmp1
          q_ad(ij,l) = q_ad(ij,l) - dyqv_ad(ij)
        ENDDO
    ENDIF
    
    areascb_ad(:) = 0.
    areasch_ad(:) = 0.
    qpns_ad = qpns_ad/areaj2
    qpsn_ad = qpsn_ad/areajjm
    DO i = 1, iim
      areascb_ad(i) = areascb_ad(i) + qpns_ad
      areasch_ad(i) = areasch_ad(i) + qpsn_ad
    ENDDO
    qpns_ad = 0.
    qpsn_ad = 0.
    
    DO i = 1, iim
      q_ad(i+ iip1,l) = q_ad(i+ iip1,l) + area(i+ iip1) * areascb_ad(i)
      q_ad(i+ ip1jm- iip1,l) = q_ad(i+ ip1jm- iip1,l) &
         + area(i+ ip1jm- iip1) * areasch_ad(i)
      areascb_ad(i) = 0.
      areasch_ad(i) = 0.
    ENDDO
    
  ENDDO
    
  
  !last bit of forward
  DO l=1,llm
    DO ij=iip2,ip1jm
      masse(ij,l) = newmasse(ij,l)
    ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE vly_ad_p_bis

SUBROUTINE vlz_ad_p(q,pente_max,masse,w,q_ad,ijb_x,ije_x)
  
  !
  !     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
  !
  !    ********************************************************************
  !     Shema  d'advection " pseudo amont " .
  !    ********************************************************************
  !    q,pbaru,pbarv,w sont des arguments d'entree  pour le s-pg ....
  !     dq 	       sont des arguments de sortie pour le s-pg ....
  !
  !
  !   --------------------------------------------------------------------
  USE Parallel
  IMPLICIT NONE
  !
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "logic.h"
  INCLUDE "comvert.h"
  INCLUDE "comconst.h"
  !
  !
  !   Arguments:
  !   ----------
  REAL :: masse(ip1jmp1,llm),pente_max
  REAL :: q(ip1jmp1,llm)
  REAL :: q_ad(ip1jmp1,llm)
  REAL :: w(ip1jmp1,llm+1)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l
  !
  REAL :: wq(ip1jmp1,llm+1),newmasse(ip1jmp1,llm)
  REAL :: wq_ad(ip1jmp1,llm+1)
  
  REAL :: dzq(ip1jmp1,llm),dzqw(ip1jmp1,llm),adzqw(ip1jmp1,llm),dzqmax(ip1jmp1,llm)
  REAL :: dzq_ad(ip1jmp1,llm),dzqw_ad(ip1jmp1,llm),adzqw_ad(ip1jmp1,llm)
  REAL :: dzqmax_ad
  REAL :: sigw
  REAL :: a, a1(ip1jmp1,llm), a2(ip1jmp1,llm), a3(ip1jmp1,llm), a_ad
  
  REAL ::      SSUM
  INTEGER :: ijb,ije,ijb_x,ije_x
  EXTERNAL :: SSUM, convflu
  EXTERNAL :: filtreg
  
  !    On oriente tout dans le sens de la pression c'est a dire dans le
  !    sens de W
  ijb=ijb_x
  ije=ije_x
  
  DO l=2,llm
    DO ij=ijb,ije
      dzqw(ij,l)=q(ij,l-1)-q(ij,l)
      IF (dzqw(ij,l) > 0.) THEN
          adzqw(ij,l)=dzqw(ij,l)
      ELSE
          adzqw(ij,l)=-dzqw(ij,l)
      ENDIF
    ENDDO
  ENDDO
  
  DO l=2,llm-1
    DO ij=ijb,ije
      IF(dzqw(ij,l)*dzqw(ij,l+1)>0.) THEN
          dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
      ELSE
          dzq(ij,l)=0.
      ENDIF
      IF (adzqw(ij,l) > adzqw(ij,l+1)) THEN
          dzqmax(ij,l)=pente_max*adzqw(ij,l+1)
      ELSE
          dzqmax(ij,l)=pente_max*adzqw(ij,l)
      ENDIF
      a = dzq(ij,l)
      a1(ij,l) = a
      IF (a1(ij,l) < 0.) THEN
          a = -a
      ENDIF
      a2(ij,l) = a
      IF (a2(ij,l) > dzqmax(ij,l)) THEN
          a = dzqmax(ij,l)
      ENDIF
      a3(ij,l) = a
      IF ((dzq(ij,l)*a3(ij,l)) < 0.) THEN
          a = -a
      ENDIF
      dzq(ij,l)=a
    ENDDO
  ENDDO
  
  DO ij=ijb,ije
    dzq(ij,1)=0.
    dzq(ij,llm)=0.
  ENDDO
  
  ! ---------------------------------------------------------------
  !   .... calcul des termes d'advection verticale  .......
  ! ---------------------------------------------------------------
  
  ! calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dq pour calculer dq
  
  DO l = 1,llm-1
    DO  ij = ijb,ije
      IF(w(ij,l+1)>0.) THEN
          sigw=w(ij,l+1)/masse(ij,l+1)
          wq(ij,l+1)=w(ij,l+1)*(q(ij,l+1)+0.5*(1.-sigw)*dzq(ij,l+1))
      ELSE
          sigw=w(ij,l+1)/masse(ij,l)
          wq(ij,l+1)=w(ij,l+1)*(q(ij,l)-0.5*(1.+sigw)*dzq(ij,l))
      ENDIF
    ENDDO
  ENDDO
  
  DO ij=ijb,ije
    wq(ij,llm+1)=0.
    wq(ij,1)=0.
  ENDDO
  
  DO l=1,llm
    DO ij=ijb,ije
      newmasse(ij,l)=masse(ij,l)+w(ij,l+1)-w(ij,l)
      q(ij,l)=(q(ij,l)*masse(ij,l)+wq(ij,l+1)-wq(ij,l)) &
         /newmasse(ij,l)
    ENDDO
  ENDDO
  
  !Start of AD
  
  wq_ad(ijb:ije,:) = 0.
  DO l=1,llm
    DO ij=ijb,ije
      wq_ad(ij,l) = wq_ad(ij,l) - q_ad(ij,l)/newmasse(ij,l)
      wq_ad(ij,l+1) = wq_ad(ij,l+1) + q_ad(ij,l)/newmasse(ij,l)
      q_ad(ij,l) = q_ad(ij,l) * masse(ij,l) / newmasse(ij,l)
    ENDDO
  ENDDO
  
  DO ij=ijb,ije
    wq_ad(ij,llm+1)=0.
    wq_ad(ij,1)=0.
  ENDDO
  
  dzq_ad(ijb:ije,:) = 0.
  DO l = 1,llm-1
    DO  ij = ijb,ije
      IF(w(ij,l+1)>0.) THEN
          sigw=w(ij,l+1)/masse(ij,l+1)
          q_ad(ij,l+1) = q_ad(ij,l+1) + wq_ad(ij,l+1) * w(ij,l+1)
          dzq_ad(ij,l+1) = dzq_ad(ij,l+1) + w(ij,l+1)*0.5*(1.-sigw) * wq_ad(ij,l+1)
          wq_ad(ij,l+1) = 0.
      ELSE
          sigw=w(ij,l+1)/masse(ij,l)
          q_ad(ij,l) = q_ad(ij,l) + wq_ad(ij,l+1) * w(ij,l+1)
          dzq_ad(ij,l) = dzq_ad(ij,l) - w(ij,l+1)*0.5*(1.+sigw) * wq_ad(ij,l+1)
          wq_ad(ij,l+1) = 0.
      ENDIF
    ENDDO
  ENDDO

  DO ij=ijb,ije
    dzq_ad(ij,1)=0.
    dzq_ad(ij,llm)=0.
  ENDDO
  
  dzqw_ad(ijb:ije,:) = 0.
  adzqw_ad(ijb:ije,:) = 0.
  DO l=2,llm-1
    DO ij=ijb,ije
      a_ad = dzq_ad(ij,l)
      dzq_ad(ij,l) = 0.
      IF ((dzq(ij,l)*a3(ij,l)) < 0.) THEN
          a_ad = - a_ad
      ENDIF
      dzqmax_ad = 0.
      IF (a2(ij,l) > dzqmax(ij,l)) THEN
          dzqmax_ad = dzqmax_ad + a_ad
          a_ad = 0.
      ENDIF
      IF (a1(ij,l) < 0.) THEN
          a_ad = - a_ad
      ENDIF
      dzq_ad(ij,l) = dzq_ad(ij,l) + a_ad
      a_ad = 0.
      IF (adzqw(ij,l) > adzqw(ij,l+1)) THEN
          adzqw_ad(ij,l+1)=adzqw_ad(ij,l+1) + pente_max*dzqmax_ad
          dzqmax_ad = 0.
      ELSE
          adzqw_ad(ij,l)=adzqw_ad(ij,l) + pente_max*dzqmax_ad
          dzqmax_ad = 0.
      ENDIF
      IF(dzqw(ij,l)*dzqw(ij,l+1) > 0.) THEN
          dzqw_ad(ij,l) = dzqw_ad(ij,l) + 0.5*dzq_ad(ij,l)
          dzqw_ad(ij,l+1) = dzqw_ad(ij,l+1) + 0.5*dzq_ad(ij,l)
          dzq_ad(ij,l) = 0.
      ELSE
          dzq_ad(ij,l) = 0.
      ENDIF
    ENDDO
  ENDDO
  
  DO l=2,llm
    DO ij=ijb,ije
      IF (dzqw(ij,l) > 0.) THEN
          dzqw_ad(ij,l)= dzqw_ad(ij,l) + adzqw_ad(ij,l)
          adzqw_ad(ij,l) = 0.
      ELSE
          dzqw_ad(ij,l)= dzqw_ad(ij,l) - adzqw_ad(ij,l)
          adzqw_ad(ij,l) = 0.
      ENDIF
      q_ad(ij,l-1) = q_ad(ij,l-1) + dzqw_ad(ij,l)
      q_ad(ij,l) = q_ad(ij,l) - dzqw_ad(ij,l)
      dzqw_ad(ij,l) = 0.
    ENDDO
  ENDDO
  
  ! last bit of forward !
  masse(ijb:ije,:) = newmasse(ijb:ije,:)
  
  RETURN
END SUBROUTINE vlz_ad_p
