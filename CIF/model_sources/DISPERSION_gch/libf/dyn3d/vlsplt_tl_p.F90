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

SUBROUTINE vlsplt_tl_p(q,pente_max,masse,w,pbaru,pbarv,pdt,q_tl)
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
  USE write_field_p
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
  REAL :: q(ip1jmp1,llm), q_tl(ip1jmp1,llm)
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
  REAL :: zq_tl(ip1jmp1,llm)
  REAL :: zzpbar, zzw
  
  INTEGER :: ijb,ije
  TYPE(request) :: MyRequest1
  TYPE(request) :: MyRequest2
  
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
  zq_tl(ijb:ije,:)=q_tl(ijb:ije,:)
  zm(ijb:ije,:)=masse(ijb:ije,:)
  
  CALL VTb(VTvlx)
  CALL vlx_tl_p(zq,pente_max,zm,mu,zq_tl,ij_begin,ij_end)
  CALL VTe(VTvlx)
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest1)
  CALL Register_Hallo(zq_tl,ip1jmp1,llm,2,2,2,2,MyRequest1)
  CALL Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest1)
  CALL Register_Hallo(mv,ip1jm,llm,1,1,1,1,MyRequest1)
  CALL SendRequest(MyRequest1)
  CALL WaitRequest(MyRequest1)
  CALL VTe(VTHallo)
  
  CALL VTb(VTvly)
  CALL vly_tl_p(zq,pente_max,zm,mv,zq_tl)
  CALL VTe(VTvly)
  
  CALL VTb(VTvlz)
  CALL vlz_tl_p(zq,pente_max,zm,mw,zq_tl,ij_begin,ij_end)
  CALL VTe(VTvlz)
  
  CALL VTb(VTHallo)
  CALL Register_Hallo(zq,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(zq_tl,ip1jmp1,llm,2,2,2,2,MyRequest2)
  CALL Register_Hallo(zm,ip1jmp1,llm,1,1,1,1,MyRequest2)
  CALL SendRequest(MyRequest2)
  CALL WaitRequest(MyRequest2)
  CALL VTe(VTHallo)
  
  CALL VTb(VTvly)
  CALL vly_tl_p(zq,pente_max,zm,mv,zq_tl)
  CALL VTe(VTvly)
  
  CALL VTb(VTvlx)
  CALL vlx_tl_p(zq,pente_max,zm,mu,zq_tl,ij_begin,ij_end)
  CALL VTe(VTvlx)
  
  ijb=ij_begin
  ije=ij_end
  
  DO l=1,llm
    DO ij=ijb,ije
      q_tl(ij,l)=zq_tl(ij,l)           
      q(ij,l)=zq(ij,l)
    ENDDO
  ENDDO
  
  
  DO l=1,llm
    DO ij=ijb,ije-iip1+1,iip1
      q_tl(ij+iim,l)=q_tl(ij,l)
      q(ij+iim,l)=q(ij,l)
    ENDDO
  ENDDO
  
  
  RETURN
END SUBROUTINE vlsplt_tl_p


SUBROUTINE vlx_tl_p(q,pente_max,masse,u_m,q_tl,ijb_x,ije_x)
  
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
  REAL :: q_tl(ip1jmp1,llm)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l,j,i,iju,ijq,indu(ip1jmp1),niju
  INTEGER :: n0,iadvplus(ip1jmp1,llm),nl(llm)
  !
  REAL :: new_m,zu_m,zdum(ip1jmp1,llm)
  REAL :: dxq(ip1jmp1,llm),dxqu(ip1jmp1)
  REAL :: adxqu(ip1jmp1),dxqmax(ip1jmp1,llm)
  REAL :: u_mq(ip1jmp1,llm)
  
  REAL :: dxq_tl(ip1jmp1,llm),dxqu_tl(ip1jmp1)
  REAL :: adxqu_tl(ip1jmp1),dxqmax_tl(ip1jmp1,llm)
  REAL :: u_mq_tl(ip1jmp1,llm)
  REAL :: a, a_tl
  
  REAL ::     SSUM
  EXTERNAL :: SSUM
  
  INTEGER :: ijb,ije,ijb_x,ije_x
  
  !   calcul de la pente a droite et a gauche de la maille
  
  IF (pente_max<=-1.e-5) THEN
      PRINT *,'pas de TL'
      STOP
  ENDIF
  ijb=ijb_x
  ije=ije_x
  
  IF (north_pole.AND.ijb==1) ijb=ijb+iip1
  IF (south_pole.AND.ije==ip1jmp1)  ije=ije-iip1
  
  
  !   calcul des pentes avec limitation, Van Leer scheme I:
  !   -----------------------------------------------------
  
  !   calcul de la pente aux points u
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)         
  DO l = 1, llm
    
    DO ij=ijb,ije-1
      dxqu_tl(ij)=q_tl(ij+1,l)-q_tl(ij,l)
      dxqu(ij)=q(ij+1,l)-q(ij,l)
    ENDDO
    DO ij=ijb+iip1-1,ije,iip1
      dxqu_tl(ij)=dxqu_tl(ij-iim)
      dxqu(ij)=dxqu(ij-iim)
    ENDDO
    
    DO ij=ijb,ije
      IF (dxqu(ij) > 0.) THEN
          adxqu_tl(ij)=dxqu_tl(ij)
          adxqu(ij)=dxqu(ij)
      ELSE
          adxqu_tl(ij)=-dxqu_tl(ij)
          adxqu(ij)=-dxqu(ij)
      ENDIF
    ENDDO
    
    !   calcul de la pente maximum dans la maille en valeur absolue
    
    DO ij=ijb+1,ije
      IF (adxqu(ij-1) < adxqu(ij)) THEN
          dxqmax_tl(ij,l)=pente_max*adxqu_tl(ij-1)
          dxqmax(ij,l)=pente_max*adxqu(ij-1)
      ELSE
          dxqmax_tl(ij,l)=pente_max*adxqu_tl(ij)
          dxqmax(ij,l)=pente_max*adxqu(ij)
      ENDIF
      
    ENDDO
    
    DO ij=ijb+iip1-1,ije,iip1
      dxqmax_tl(ij-iim,l)=dxqmax_tl(ij,l)
      dxqmax(ij-iim,l)=dxqmax(ij,l)
    ENDDO
    
    DO ij=ijb+1,ije
      IF(dxqu(ij-1)*dxqu(ij)>0) THEN
          dxq_tl(ij,l)=dxqu_tl(ij-1)+dxqu_tl(ij)
          dxq(ij,l)=dxqu(ij-1)+dxqu(ij)
      ELSE
          !   extremum local
          dxq_tl(ij,l)=0.
          dxq(ij,l)=0.
      ENDIF
      dxq_tl(ij,l)=0.5*dxq_tl(ij,l)
      dxq(ij,l)=0.5*dxq(ij,l)
      a_tl = dxq_tl(ij,l)
      a = dxq(ij,l)
      IF (a < 0.) THEN
          a_tl = -a_tl
          a = -a
      ENDIF
      IF (a > dxqmax(ij,l)) THEN
          a_tl = dxqmax_tl(ij,l)
          a = dxqmax(ij,l)
      ENDIF
      IF ((dxq(ij,l)*a) < 0.) THEN
          a_tl = -a_tl
          a = -a
      ENDIF
      dxq_tl(ij,l)=a_tl
      dxq(ij,l)=a
    ENDDO
    
  ENDDO ! l=1,llm
  
  !   bouclage de la pente en iip1:
  !   -----------------------------
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb+iip1-1,ije,iip1
      dxq_tl(ij-iim,l)=dxq_tl(ij,l)
      dxq(ij-iim,l)=dxq(ij,l)
    ENDDO
    DO ij=ijb,ije
      iadvplus(ij,l)=0
    ENDDO
    
  ENDDO
  !$OMP END DO NOWAIT

  !   calcul des flux a gauche et a droite

  !   on cumule le flux correspondant a toutes les mailles dont la masse
  !   au travers de la paroi pENDant le pas de temps.
  
  DO l=1,llm
    DO ij=ijb,ije-1
      IF (u_m(ij,l)>0.) THEN
          zdum(ij,l)=1.-u_m(ij,l)/masse(ij,l)
          u_mq_tl(ij,l)=u_m(ij,l)*q_tl(ij,l) + &
             u_m(ij,l)*0.5*zdum(ij,l)*dxq_tl(ij,l)
          u_mq(ij,l)=u_m(ij,l)*(q(ij,l)+0.5*zdum(ij,l)*dxq(ij,l))
      ELSE
          zdum(ij,l)=1.+u_m(ij,l)/masse(ij+1,l)
          u_mq_tl(ij,l)=u_m(ij,l)*q_tl(ij+1,l)- &
             u_m(ij,l)*0.5*zdum(ij,l)*dxq_tl(ij+1,l)
          u_mq(ij,l)=u_m(ij,l)*(q(ij+1,l)-0.5*zdum(ij,l)*dxq(ij+1,l))
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  !   detection des points ou on advecte plus que la masse de la
  !   maille
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb,ije-1
      IF(zdum(ij,l)<0) THEN
          iadvplus(ij,l)=1
          u_mq_tl(ij,l)=0.
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
          u_mq_tl(ij,l)=0.
          u_mq(ij,l)=0.
          IF(zu_m>0.) THEN
              ijq=ij
              i=ijq-(j-1)*iip1
              !   accumulation pour les mailles completements advectees
              DO WHILE(zu_m>masse(ijq,l))
                u_mq_tl(ij,l)=u_mq_tl(ij,l)+q_tl(ijq,l)*masse(ijq,l)
                u_mq(ij,l)=u_mq(ij,l)+q(ijq,l)*masse(ijq,l)
                zu_m=zu_m-masse(ijq,l)
                i=MOD(i-2+iim,iim)+1
                ijq=(j-1)*iip1+i
              ENDDO
              !   ajout de la maille non completement advectee
              u_mq_tl(ij,l)=u_mq_tl(ij,l)+zu_m* &
                 (q_tl(ijq,l)+0.5*(1.-zu_m/masse(ijq,l))*dxq_tl(ijq,l))
              u_mq(ij,l)=u_mq(ij,l)+zu_m* &
                 (q(ijq,l)+0.5*(1.-zu_m/masse(ijq,l))*dxq(ijq,l))
          ELSE
              ijq=ij+1
              i=ijq-(j-1)*iip1
              !   accumulation pour les mailles completements advectees
              DO WHILE(-zu_m>masse(ijq,l))
                u_mq_tl(ij,l)=u_mq_tl(ij,l)-q_tl(ijq,l)*masse(ijq,l)
                u_mq(ij,l)=u_mq(ij,l)-q(ijq,l)*masse(ijq,l)
                zu_m=zu_m+masse(ijq,l)
                i=MOD(i,iim)+1
                ijq=(j-1)*iip1+i
              ENDDO
              !   ajout de la maille non completement advectee
              u_mq_tl(ij,l)=u_mq_tl(ij,l)+zu_m*(q_tl(ijq,l)- &
                 0.5*(1.+zu_m/masse(ijq,l))*dxq_tl(ijq,l))
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
      u_mq_tl(ij,l)=u_mq_tl(ij-iim,l)
      u_mq(ij,l)=u_mq(ij-iim,l)
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  !   calcul des tENDances
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb+1,ije
      new_m=masse(ij,l)+u_m(ij-1,l)-u_m(ij,l)
      q_tl(ij,l)=(q_tl(ij,l)*masse(ij,l)+ &
         u_mq_tl(ij-1,l)-u_mq_tl(ij,l)) &
         /new_m
      q(ij,l)=(q(ij,l)*masse(ij,l)+ &
         u_mq(ij-1,l)-u_mq(ij,l)) &
         /new_m
      masse(ij,l)=new_m
    ENDDO
    !   ModIF Fred 22 03 96 correction d'un bug (les scopy ci-dessous)
    DO ij=ijb+iip1-1,ije,iip1
      q_tl(ij-iim,l)=q_tl(ij,l)
      q(ij-iim,l)=q(ij,l)
      masse(ij-iim,l)=masse(ij,l)
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  
  RETURN
END SUBROUTINE vlx_tl_p
SUBROUTINE vly_tl_p(q,pente_max,masse,masse_adv_v,q_tl)
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
  REAL :: q_tl(ip1jmp1,llm)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: i,ij,l
  !
  REAL :: areaj2,areajjm,areascb(iim),areasch(iim)
  REAL :: areascb_tl(iim),areasch_tl(iim)
  REAL :: dyq(ip1jmp1,llm),dyqv(ip1jm)
  REAL :: dyq_tl(ip1jmp1,llm),dyqv_tl(ip1jm)
  REAL :: adyqv(ip1jm),dyqmax(ip1jmp1)
  REAL :: adyqv_tl(ip1jm),dyqmax_tl(ip1jmp1)
  REAL :: qbyv(ip1jm,llm)
  REAL :: qbyv_tl(ip1jm,llm)
  
  REAL :: qpns,qpsn,dyn1,dys1,dyn2,dys2,newmasse
  REAL :: qpns_tl, qpsn_tl, dyn1_tl,dys1_tl, dyn2_tl, dys2_tl
  LOGICAL :: first
  SAVE :: first
  !$OMP THREADPRIVATE(first)
  
  REAL :: convpn,convps,convmpn,convmps
  REAL :: convpn_tl,convps_tl
  REAL :: massepn,masseps,qpn,qps
  REAL :: qpn_tl,qps_tl
  REAL :: sinlon(iip1),sinlondlon(iip1)
  REAL :: coslon(iip1),coslondlon(iip1)
  SAVE :: sinlon,coslon,sinlondlon,coslondlon
  !$OMP THREADPRIVATE(sinlon,coslon,sinlondlon,coslondlon)
  SAVE :: areaj2,areajjm
  !$OMP THREADPRIVATE(areaj2,areajjm)
  REAL :: a, a_tl
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
    !
    !   On commence par calculer la valeur du traceur moyenne sur le premier cercle
    !   de latitude autour du pole (qpns pour le pole nord et qpsn pour
    !    le pole nord) qui sera utilisee pour evaluer les pentes au pole.
    
    IF (north_pole) THEN
        DO i = 1, iim
          areascb_tl(i) = area(i+ iip1) * q_tl(i+ iip1,l)
          areascb(i) = area(i+ iip1) * q(i+ iip1,l)
        ENDDO
        qpns_tl   = SSUM( iim,  areascb_tl ,1 ) / areaj2
        qpns   = SSUM( iim,  areascb ,1 ) / areaj2
    ENDIF
    
    IF (south_pole) THEN
        DO i = 1, iim
          areasch_tl(i) = area(i+ ip1jm- iip1) * q_tl(i+ ip1jm- iip1,l)
          areasch(i) = area(i+ ip1jm- iip1) * q(i+ ip1jm- iip1,l)
        ENDDO
        qpsn_tl   = SSUM( iim,  areasch_tl ,1 ) / areajjm
        qpsn   = SSUM( iim,  areasch ,1 ) / areajjm
    ENDIF
    
    
    
    !   calcul des pentes aux points v
    
    ijb=ij_begin-2*iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyqv_tl(ij)=q_tl(ij,l)-q_tl(ij+iip1,l)
      dyqv(ij)=q(ij,l)-q(ij+iip1,l)
      IF (dyqv(ij) > 0.) THEN
          adyqv_tl(ij)=dyqv_tl(ij)
          adyqv(ij)=dyqv(ij)
      ELSE
          adyqv_tl(ij)=-dyqv_tl(ij)
          adyqv(ij)=-dyqv(ij)
      ENDIF
    ENDDO
    
    !   calcul des pentes aux points scalaires
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1
    
    DO ij=ijb,ije
      dyq_tl(ij,l)=.5*(dyqv_tl(ij-iip1)+dyqv_tl(ij))
      dyq(ij,l)=.5*(dyqv(ij-iip1)+dyqv(ij))
      IF (adyqv(ij-iip1) > adyqv(ij)) THEN
          dyqmax_tl(ij)=adyqv_tl(ij)
          dyqmax(ij)=adyqv(ij)
      ELSE
          dyqmax_tl(ij)=adyqv_tl(ij-iip1)
          dyqmax(ij)=adyqv(ij-iip1)
      ENDIF
      dyqmax_tl(ij)=pente_max*dyqmax_tl(ij)
      dyqmax(ij)=pente_max*dyqmax(ij)
    ENDDO
    
    !   calcul des pentes aux poles
    IF (north_pole) THEN
        DO ij=1,iip1
          dyq_tl(ij,l)=qpns_tl-q_tl(ij+iip1,l)
          dyq(ij,l)=qpns-q(ij+iip1,l)
        ENDDO
    ENDIF
    
    IF (south_pole) THEN
        DO ij=1,iip1
          dyq_tl(ip1jm+ij,l)=q_tl(ip1jm+ij-iip1,l)-qpsn_tl
          dyq(ip1jm+ij,l)=q(ip1jm+ij-iip1,l)-qpsn
        ENDDO
    ENDIF
    
    !   filtrage de la derivee
    dyn1_tl=0.
    dyn1=0.
    dys1_tl=0.
    dys1=0.
    dyn2_tl=0.
    dyn2=0.
    dys2_tl=0.
    dys2=0.
    DO ij=1,iim
      dyn1_tl=dyn1_tl+sinlondlon(ij)*dyq_tl(ij,l)
      dyn1=dyn1+sinlondlon(ij)*dyq(ij,l)
      dys1_tl=dys1_tl+sinlondlon(ij)*dyq_tl(ip1jm+ij,l)
      dys1=dys1+sinlondlon(ij)*dyq(ip1jm+ij,l)
      dyn2_tl=dyn2_tl+coslondlon(ij)*dyq_tl(ij,l)
      dyn2=dyn2+coslondlon(ij)*dyq(ij,l)
      dys2_tl=dys2_tl+coslondlon(ij)*dyq_tl(ip1jm+ij,l)
      dys2=dys2+coslondlon(ij)*dyq(ip1jm+ij,l)
    ENDDO
    DO ij=1,iip1
      dyq_tl(ij,l)=dyn1_tl*sinlon(ij)+dyn2_tl*coslon(ij)
      dyq(ij,l)=dyn1*sinlon(ij)+dyn2*coslon(ij)
      dyq_tl(ip1jm+ij,l)=dys1_tl*sinlon(ij)+dys2_tl*coslon(ij)
      dyq(ip1jm+ij,l)=dys1*sinlon(ij)+dys2*coslon(ij)
    ENDDO
    
    !   calcul des pentes limites aux poles
    
    IF (north_pole) THEN
        DO ij=1,iip1
          dyq_tl(ij,l)=0.
          dyq(ij,l)=0.
        END DO
    ENDIF
    IF (south_pole) THEN
        DO ij=1,iip1
          dyq_tl(ip1jm+ij,l)=0.
          dyq(ip1jm+ij,l)=0.
        END DO
    ENDIF
    
    
    !   calcul des pentes limitees
    ijb=ij_begin-iip1
    ije=ij_end+iip1
    IF (north_pole) ijb=ij_begin+iip1
    IF (south_pole)  ije=ij_end-iip1

    DO ij=ijb,ije
      IF(dyqv(ij)*dyqv(ij-iip1)>0.) THEN
          a_tl = dyq_tl(ij,l)
          a = dyq(ij,l)
          IF (a < 0.) THEN
              a_tl = -a_tl
              a = -a
          ENDIF
          IF (a > dyqmax(ij)) THEN
              a_tl = dyqmax_tl(ij)
              a = dyqmax(ij)
          ENDIF
          IF ((dyq(ij,l)*a) < 0.) THEN
              a_tl = -a_tl
              a = -a
          ENDIF
          dyq_tl(ij,l)=a_tl
          dyq(ij,l)=a
      ELSE
          dyq_tl(ij,l)=0.
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
          qbyv_tl(ij,l)=q_tl(ij+iip1,l)+dyq_tl(ij+iip1,l)* &
             0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
          qbyv(ij,l)=q(ij+iip1,l)+dyq(ij+iip1,l)* &
             0.5*(1.-masse_adv_v(ij,l)/masse(ij+iip1,l))
      ELSE
          qbyv_tl(ij,l)=q_tl(ij,l)-dyq_tl(ij,l)* &
             0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
          qbyv(ij,l)=q(ij,l)-dyq(ij,l)* &
             0.5*(1.+masse_adv_v(ij,l)/masse(ij,l))
      ENDIF
      qbyv_tl(ij,l)=masse_adv_v(ij,l)*qbyv_tl(ij,l)
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
      q_tl(ij,l)=(q_tl(ij,l)*masse(ij,l)+qbyv_tl(ij,l)-qbyv_tl(ij-iip1,l)) &
         /newmasse
      q(ij,l)=(q(ij,l)*masse(ij,l)+qbyv(ij,l)-qbyv(ij-iip1,l)) &
         /newmasse
      masse(ij,l)=newmasse
    ENDDO
    
    IF (north_pole) THEN
        convpn_tl= SSUM(iim,qbyv_tl(1,l),1)
        convpn=SSUM(iim,qbyv(1,l),1)
        convmpn=ssum(iim,masse_adv_v(1,l),1)
        massepn=ssum(iim,masse(1,l),1)
        qpn_tl=0.
        qpn=0.
        DO ij=1,iim
          qpn_tl=qpn_tl+masse(ij,l)*q_tl(ij,l)
          qpn=qpn+masse(ij,l)*q(ij,l)
        ENDDO
        qpn_tl=(qpn_tl+convpn_tl)/(massepn+convmpn)
        qpn=(qpn+convpn)/(massepn+convmpn)
        DO ij=1,iip1
          q_tl(ij,l)=qpn_tl
          q(ij,l)=qpn
        ENDDO
    ENDIF
    
    IF (south_pole) THEN
        
        convps_tl=-SSUM(iim,qbyv_tl(ip1jm-iim,l),1)
        convps=-SSUM(iim,qbyv(ip1jm-iim,l),1)
        convmps=-ssum(iim,masse_adv_v(ip1jm-iim,l),1)
        masseps=ssum(iim, masse(ip1jm+1,l),1)
        qps_tl=0.
        qps=0.
        DO ij = ip1jm+1,ip1jmp1-1
          qps_tl=qps_tl+masse(ij,l)*q_tl(ij,l)
          qps=qps+masse(ij,l)*q(ij,l)
        ENDDO
        qps_tl=(qps_tl+convps_tl)/(masseps+convmps)
        qps=(qps+convps)/(masseps+convmps)
        DO ij=ip1jm+1,ip1jmp1
          q_tl(ij,l)=qps_tl
          q(ij,l)=qps
        ENDDO
    ENDIF
  ENDDO
  !$OMP END DO NOWAIT
  
  RETURN
END SUBROUTINE vly_tl_p
SUBROUTINE vlz_tl_p(q,pente_max,masse,w,q_tl,ijb_x,ije_x)
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
  REAL :: q_tl(ip1jmp1,llm)
  REAL :: w(ip1jmp1,llm+1)
  !
  !      Local 
  !   ---------
  !
  INTEGER :: ij,l
  !
  REAL,SAVE :: wq(ip1jmp1,llm+1)
  REAL :: newmasse
  REAL,SAVE :: wq_tl(ip1jmp1,llm+1)
  
  REAL,SAVE :: dzq(ip1jmp1,llm),dzqw(ip1jmp1,llm),adzqw(ip1jmp1,llm)
  REAL,SAVE :: dzq_tl(ip1jmp1,llm),dzqw_tl(ip1jmp1,llm),adzqw_tl(ip1jmp1,llm)
  REAL :: dzqmax,dzqmax_tl
  REAL :: sigw
  REAL :: a, a_tl
  
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
      dzqw_tl(ij,l)=q_tl(ij,l-1)-q_tl(ij,l)
      dzqw(ij,l)=q(ij,l-1)-q(ij,l)
      IF (dzqw(ij,l) > 0.) THEN
          adzqw_tl(ij,l)=dzqw_tl(ij,l)
          adzqw(ij,l)=dzqw(ij,l)
      ELSE
          adzqw_tl(ij,l)=-dzqw_tl(ij,l)
          adzqw(ij,l)=-dzqw(ij,l)
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=2,llm-1
    DO ij=ijb,ije
      IF(dzqw(ij,l)*dzqw(ij,l+1)>0.) THEN
          dzq_tl(ij,l)=0.5*(dzqw_tl(ij,l)+dzqw_tl(ij,l+1))
          dzq(ij,l)=0.5*(dzqw(ij,l)+dzqw(ij,l+1))
      ELSE
          dzq_tl(ij,l)=0.
          dzq(ij,l)=0.
      ENDIF
      IF (adzqw(ij,l) > adzqw(ij,l+1)) THEN
          dzqmax_tl=pente_max*adzqw_tl(ij,l+1)
          dzqmax=pente_max*adzqw(ij,l+1)
      ELSE
          dzqmax_tl=pente_max*adzqw_tl(ij,l)
          dzqmax=pente_max*adzqw(ij,l)
      ENDIF
      a_tl = dzq_tl(ij,l)
      a = dzq(ij,l)
      IF (a < 0.) THEN
          a_tl = -a_tl
          a = -a
      ENDIF
      IF (a > dzqmax) THEN
          a_tl = dzqmax_tl
          a = dzqmax
      ENDIF
      IF ((dzq(ij,l)*a) < 0.) THEN
          a_tl = -a_tl
          a = -a
      ENDIF
      dzq_tl(ij,l)=a_tl
      dzq(ij,l)=a
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  !$OMP MASTER
  DO ij=ijb,ije
    dzq_tl(ij,1)=0.
    dzq(ij,1)=0.
    dzq_tl(ij,llm)=0.
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
          wq_tl(ij,l+1)=w(ij,l+1)*(q_tl(ij,l+1)+0.5*(1.-sigw)*dzq_tl(ij,l+1))
          wq(ij,l+1)=w(ij,l+1)*(q(ij,l+1)+0.5*(1.-sigw)*dzq(ij,l+1))
      ELSE
          sigw=w(ij,l+1)/masse(ij,l)
          wq_tl(ij,l+1)=w(ij,l+1)*(q_tl(ij,l)-0.5*(1.+sigw)*dzq_tl(ij,l))
          wq(ij,l+1)=w(ij,l+1)*(q(ij,l)-0.5*(1.+sigw)*dzq(ij,l))
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  !$OMP MASTER
  DO ij=ijb,ije
    wq_tl(ij,llm+1)=0.
    wq(ij,llm+1)=0.
    wq_tl(ij,1)=0.
    wq(ij,1)=0.
  ENDDO
  !$OMP END MASTER
  !$OMP BARRIER
  
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
    DO ij=ijb,ije
      newmasse=masse(ij,l)+w(ij,l+1)-w(ij,l)
      q_tl(ij,l)=(q_tl(ij,l)*masse(ij,l)+wq_tl(ij,l+1)-wq_tl(ij,l)) &
         /newmasse
      q(ij,l)=(q(ij,l)*masse(ij,l)+wq(ij,l+1)-wq(ij,l)) &
         /newmasse
      masse(ij,l)=newmasse
    ENDDO
  ENDDO
  !$OMP END DO NOWAIT
  
  
  RETURN
END SUBROUTINE vlz_tl_p
