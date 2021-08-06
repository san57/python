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

SUBROUTINE vlsplt_p(q,pente_max,masse,w,pbaru,pbarv,pdt)
  !
  !     Auteurs:   P.Le Van, F.Hourdin, F.Forget 
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
  USE write_field_p
  USE times
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
  REAL :: q(ip1jmp1,llm)
  REAL :: w(ip1jmp1,llm),pdt
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l
  !
  REAL :: zm(ip1jmp1,llm)
  REAL :: mu(ip1jmp1,llm)
  REAL :: mv(ip1jm,llm)
  REAL :: mw(ip1jmp1,llm+1)
  REAL :: zq(ip1jmp1,llm)
  REAL :: zzpbar, zzw
  
  INTEGER :: ijb,ije
  TYPE(request),SAVE :: MyRequest1
  TYPE(request),SAVE :: MyRequest2
  
  
  CALL resume_timer(timer_vanleer)
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
  
  CALL suspend_timer(timer_vanleer)
  CALL VTb(VTvlx)
  CALL vlx_p(zq,pente_max,zm,mu,ij_begin,ij_end)
  CALL VTe(VTvlx)
  CALL resume_timer(timer_vanleer)
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest1)
  CALL Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest1)
  CALL Register_Hallo(mv,ip1jm,llm,1,1,1,1,MyRequest1)
  CALL SendRequest(MyRequest1)
  CALL WaitRequest(MyRequest1)
  CALL VTe(VTHallo)
  
  
  CALL suspend_timer(timer_vanleer)
  
  CALL VTb(VTvly)
  CALL vly_p(zq,pente_max,zm,mv)
  CALL VTe(VTvly)
  
  CALL VTb(VTvlz)
  CALL vlz_p(zq,pente_max,zm,mw,ij_begin,ij_end)
  CALL VTe(VTvlz)
  CALL resume_timer(timer_vanleer)
  
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest2)
  CALL SendRequest(MyRequest2)
  CALL WaitRequest(MyRequest2)
  CALL VTe(VTHallo)
  
  CALL suspend_timer(timer_vanleer)
  CALL VTb(VTvly)
  CALL vly_p(zq,pente_max,zm,mv)
  CALL VTe(VTvly)
  
  CALL VTb(VTvlx)
  CALL vlx_p(zq,pente_max,zm,mu,ij_begin,ij_end)
  CALL VTe(VTvlx)
  
  CALL resume_timer(timer_vanleer)
  
  
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
  
  
  CALL suspend_timer(timer_vanleer)
  
  RETURN
END SUBROUTINE vlsplt_p


SUBROUTINE vlx_p(q,pente_max,masse,u_m,ijb_x,ije_x)
  
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
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
  INTEGER :: n0,iadvplus(ip1jmp1,llm),nl(llm)
  !
  REAL :: new_m,zu_m,zdum(ip1jmp1,llm)
  REAL :: dxq(ip1jmp1,llm),dxqu(ip1jmp1)
  REAL :: zz(ip1jmp1)
  REAL :: adxqu(ip1jmp1),dxqmax(ip1jmp1,llm)
  REAL :: u_mq(ip1jmp1,llm)
  
  REAL ::     SSUM
  EXTERNAL :: SSUM
  
  INTEGER :: ijb,ije,ijb_x,ije_x
  
  !   calcul de la pente a droite et a gauche de la maille
  
  ijb=ijb_x
  ije=ije_x
  
  IF (north_pole.AND.ijb==1) ijb=ijb+iip1
  IF (south_pole.AND.ije==ip1jmp1)  ije=ije-iip1
  
  IF (pente_max>-1.e-5) THEN
      
      !   calcul des pentes avec limitation, Van Leer scheme I:
      !   -----------------------------------------------------
      
      !   calcul de la pente aux points u
      !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)         
      DO l = 1, llm
        
        DO ij=ijb,ije-1
          dxqu(ij)=q(ij+1,l)-q(ij,l)
          !              IF(u_m(ij,l).LT.0.) STOP'limx n admet pas les U<0'
          !              sigu(ij)=u_m(ij,l)/masse(ij,l)
        ENDDO
        DO ij=ijb+iip1-1,ije,iip1
          dxqu(ij)=dxqu(ij-iim)
          !              sigu(ij)=sigu(ij-iim)
        ENDDO
        
        DO ij=ijb,ije
          adxqu(ij)=ABS(dxqu(ij))
        ENDDO
        
        !   calcul de la pente maximum dans la maille en valeur absolue
        
        DO ij=ijb+1,ije
          dxqmax(ij,l)=pente_max* &
             MIN(adxqu(ij-1),adxqu(ij))
          ! limitation subtile
          !  MIN(adxqu(ij-1)/sigu(ij-1),adxqu(ij)/(1.-sigu(ij)))
          
          
        ENDDO
        
        DO ij=ijb+iip1-1,ije,iip1
          dxqmax(ij-iim,l)=dxqmax(ij,l)
        ENDDO
        
        DO ij=ijb+1,ije
#ifdef CRAY
          dxq(ij,l)=cvmgp(dxqu(ij-1)+dxqu(ij),0.,dxqu(ij-1)*dxqu(ij))
#else
          IF(dxqu(ij-1)*dxqu(ij)>0) THEN
              dxq(ij,l)=dxqu(ij-1)+dxqu(ij)
          ELSE
              !   extremum local
              dxq(ij,l)=0.
          ENDIF
#endif
          dxq(ij,l)=0.5*dxq(ij,l)
          dxq(ij,l)=SIGN(MIN(ABS(dxq(ij,l)),dxqmax(ij,l)),dxq(ij,l))
        ENDDO
        
      ENDDO ! l=1,llm
      !$OMP END DO NOWAIT
      
  ELSE ! (pente_max.lt.-1.e-5)
      
      !   Pentes produits:
      !   ----------------
      !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llm
        DO ij=ijb,ije-1
          dxqu(ij)=q(ij+1,l)-q(ij,l)
        ENDDO
        DO ij=ijb+iip1-1,ije,iip1
          dxqu(ij)=dxqu(ij-iim)
        ENDDO
        
        DO ij=ijb+1,ije
          zz(ij)=dxqu(ij-1)*dxqu(ij)
          zz(ij)=zz(ij)+zz(ij)
          IF(zz(ij)>0) THEN
              dxq(ij,l)=zz(ij)/(dxqu(ij-1)+dxqu(ij))
          ELSE
              !   extremum local
              dxq(ij,l)=0.
          ENDIF
        ENDDO
        
      ENDDO
      !$OMP END DO NOWAIT
  ENDIF ! (pente_max.lt.-1.e-5)
  
  !   bouclage de la pente en iip1:
  !   -----------------------------
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb+iip1-1,ije,iip1
      dxq(ij-iim,l)=dxq(ij,l)
    ENDDO
    DO ij=ijb,ije
      iadvplus(ij,l)=0
    ENDDO
    
  ENDDO
  !$OMP END DO NOWAIT
  
  !   calcul des flux a gauche et a droite
  
#ifdef CRAY
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb,ije-1
      zdum(ij,l)=cvmgp(1.-u_m(ij,l)/masse(ij,l), &
         1.+u_m(ij,l)/masse(ij+1,l), &
         u_m(ij,l))
      zdum(ij,l)=0.5*zdum(ij,l)
      u_mq(ij,l)=cvmgp( &
         q(ij,l)+zdum(ij,l)*dxq(ij,l), &
         q(ij+1,l)-zdum(ij,l)*dxq(ij+1,l), &
         u_m(ij,l))
      u_mq(ij,l)=u_m(ij,l)*u_mq(ij,l)
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
#else
  !   on cumule le flux correspondant a toutes les mailles dont la masse
  !   au travers de la paroi pENDant le pas de temps.
  !	PRINT*,'Cumule ....'
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
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
  !$OMP END DO NOWAIT
#endif
  !   detection des points ou on advecte plus que la masse de la
  !   maille
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb,ije-1
      IF(zdum(ij,l)<0) THEN
          iadvplus(ij,l)=1
          u_mq(ij,l)=0.
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb+iip1-1,ije,iip1
      iadvplus(ij,l)=iadvplus(ij-iim,l)
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  
  !   traitement special pour le cas ou on advecte en longitude plus que le
  !   contenu de la maille.
  !   cette partie est mal vectorisee.
  
  !  calcul du nombre de maille sur lequel on advecte plus que la maille.
  
  n0=0
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    nl(l)=0
    DO ij=ijb,ije
      nl(l)=nl(l)+iadvplus(ij,l)
    ENDDO
    n0=n0+nl(l)
  ENDDO
  !$OMP END DO NOWAIT
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    IF(nl(l)>0) THEN
        iju=0
        !   indicage des mailles concernees par le traitement special
        DO ij=ijb,ije
          IF(iadvplus(ij,l)==1.AND.MOD(ij,iip1)/=0) THEN
              iju=iju+1
              indu(iju)=ij
          ENDIF
        ENDDO
        niju=iju
        
        !  traitement des mailles
        DO iju=1,niju
          ij=indu(iju)
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
  !$OMP END DO NOWAIT
  
  
  !   bouclage en latitude
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb+iip1-1,ije,iip1
      u_mq(ij,l)=u_mq(ij-iim,l)
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  !   calcul des tENDances
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb+1,ije
      new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
      q(ij,l)=(q(ij,l)*masse(ij,l)+ &
         u_mq(ij-1,l)-u_mq(ij,l)) &
         /new_m
      masse(ij,l)=new_m
    ENDDO
    !   ModIF Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
    DO ij=ijb+iip1-1,ije,iip1
      q(ij-iim,l)=q(ij,l)
      masse(ij-iim,l)=masse(ij,l)
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  RETURN
END SUBROUTINE vlx_p


SUBROUTINE vly_p(q,pente_max,masse,masse_adv_v)
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
  !
  !      Local 
  !   ---------
  !
  INTEGER :: i,ij,l
  !
  REAL :: areaj2,areajjm,areascb(iim),areasch(iim)
  REAL :: dyq(ip1jmp1,llm),dyqv(ip1jm)
  REAL :: adyqv(ip1jm),dyqmax(ip1jmp1)
  REAL :: qbyv(ip1jm,llm)
  
  REAL :: qpns,qpsn,dyn1,dys1,dyn2,dys2,newmasse
  LOGICAL :: first
  SAVE :: first
  !$OMP THREADPRIVATE(first)
  
  REAL :: convpn,convps,convmpn,convmps
  REAL :: massepn,masseps,qpn,qps
  REAL :: sinlon(iip1),sinlondlon(iip1)
  REAL :: coslon(iip1),coslondlon(iip1)
  SAVE :: sinlon,coslon,sinlondlon,coslondlon
  !$OMP THREADPRIVATE(sinlon,coslon,sinlondlon,coslondlon)
  SAVE :: areaj2,areajjm
  !$OMP THREADPRIVATE(areaj2,areajjm)
  !
  !
  REAL ::     SSUM
  EXTERNAL :: SSUM
  
  DATA first/.TRUE./
  INTEGER :: ijb,ije
  
  IF(first) THEN
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
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
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
      dyqv(ij)=q(ij,l)-q(ij+iip1,l)
      adyqv(ij)=ABS(dyqv(ij))
    ENDDO
    
    !   calcul des pentes aux points scalareas
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyq(ij,l)=.5*(dyqv(ij-iip1)+dyqv(ij))
      dyqmax(ij)=MIN(adyqv(ij-iip1),adyqv(ij))
      dyqmax(ij)=pente_max*dyqmax(ij)
    ENDDO
    
    !   calcul des pentes aux poles
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
    
    !   filtrage de la derivee
    
    !   calcul des pentes limites aux poles
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      IF(dyqv(ij)*dyqv(ij-iip1)>0.) THEN
          dyq(ij,l)=SIGN(MIN(ABS(dyq(ij,l)),dyqmax(ij)),dyq(ij,l))
      ELSE
          dyq(ij,l)=0.
      ENDIF
    ENDDO
    
  ENDDO
  !$OMP END DO NOWAIT
  
  ijb=ij_begin-iip1
  ije=ij_end
  IF (north_pole) ijb=ij_begin
  IF (south_pole)  ije=ij_end-iip1
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
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
  !$OMP END DO NOWAIT
  
  ijb=ij_begin
  ije=ij_end
  IF (north_pole) ijb=ij_begin+iip1
  IF (south_pole)  ije=ij_end-iip1
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
  DO l=1,llm
    DO ij=ijb,ije
      newmasse=masse(ij,l) &
         +masse_adv_v(ij,l)-masse_adv_v(ij-iip1,l)
      
      q(ij,l)=(q(ij,l)*masse(ij,l)+qbyv(ij,l)-qbyv(ij-iip1,l)) &
         /newmasse
      masse(ij,l)=newmasse
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
  ENDDO
  !$OMP END DO NOWAIT
  
  RETURN
END SUBROUTINE vly_p



SUBROUTINE vlz_p(q,pente_max,masse,w,ijb_x,ije_x)
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
  REAL :: w(ip1jmp1,llm+1)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l
  !
  REAL,SAVE :: wq(ip1jmp1,llm+1)
  REAL :: newmasse
  
  REAL,SAVE :: dzq(ip1jmp1,llm),dzqw(ip1jmp1,llm),adzqw(ip1jmp1,llm)
  REAL :: dzqmax
  REAL :: sigw
  
  REAL ::     SSUM
  EXTERNAL :: SSUM
  
  INTEGER :: ijb,ije,ijb_x,ije_x
  !    On oriente tout dans le sens de la pression c'est a dire dans le
  !    sens de W
  
  ijb=ijb_x
  ije=ije_x
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
  DO l=2,llm
    DO ij=ijb,ije
      dzqw(ij,l)=q(ij,l-1)-q(ij,l)
      adzqw(ij,l)=ABS(dzqw(ij,l))
    ENDDO
  ENDDO
  !$OMP END DO
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=2,llm-1
    DO ij=ijb,ije
#ifdef CRAY
      dzq(ij,l)=0.5* &
         cvmgp(dzqw(ij,l)+dzqw(ij,l+1),0.,dzqw(ij,l)*dzqw(ij,l+1))
#else
      IF(dzqw(ij,l)*dzqw(ij,l+1)>0.) THEN
          dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
      ELSE
          dzq(ij,l)=0.
      ENDIF
#endif
      dzqmax=pente_max*MIN(adzqw(ij,l),adzqw(ij,l+1))
      dzq(ij,l)=SIGN(MIN(ABS(dzq(ij,l)),dzqmax),dzq(ij,l))
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  !$OMP MASTER
  DO ij=ijb,ije
    dzq(ij,1)=0.
    dzq(ij,llm)=0.
  ENDDO
  !$OMP END MASTER
  !$OMP BARRIER
  ! ---------------------------------------------------------------
  !   .... calcul des termes d'advection verticale  .......
  ! ---------------------------------------------------------------
  
  ! calcul de  - d( q   * w )/ d(sigma)    qu'on ajoute a  dq pour calculer dq
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
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
  !$OMP END DO NOWAIT
  
  !$OMP MASTER
  DO ij=ijb,ije
    wq(ij,llm+1)=0.
    wq(ij,1)=0.
  ENDDO
  !$OMP END MASTER
  !$OMP BARRIER
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb,ije
      newmasse=masse(ij,l)+w(ij,l+1)-w(ij,l)
      q(ij,l)=(q(ij,l)*masse(ij,l)+wq(ij,l+1)-wq(ij,l)) &
         /newmasse
      masse(ij,l)=newmasse
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  
  RETURN
END SUBROUTINE vlz_p
