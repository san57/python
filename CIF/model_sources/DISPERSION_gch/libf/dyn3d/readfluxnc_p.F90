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

SUBROUTINE readfluxnc_p(irec,massem,pbarun,pbarvn,wn,tetan,phin, &
   nrec,areafi,phisfi, &
   t,zmfu, zmfd, zde_u,zen_d, coefkz, &
   dake, mpke, phike, updke, dndke, wghtke,&
   entr_therm,fm_therm,&
   yu1,yv1,ftsol,pctsrf,phis)
      
  USE dimphy
  USE vampir
  IMPLICIT NONE
  
 include "dimensions.h"
 include "paramet.h"
 include "logic.h"
 include "comvert.h"
 include "comconst.h"
 include "comgeom2.h"
 include "tracstoke.h"
  
  INTEGER :: irec,nrec
  
  INTEGER :: ig,l
  INTEGER,parameter :: nbsrf=4 ! nombre de sous-fractions pour une maille
  
  ! Convection
  REAL :: zmfd(klon,llm),zen_d(klon,llm)
  REAL :: zmfu(klon,llm),zde_u(klon,llm)
  REAL :: mfd(klon,llm),en_d(klon,llm)
  REAL :: mfu(klon,llm),de_u(klon,llm)
  !Convection Kerry Emanuel
  REAL :: dake(klon,llm),mpke(klon,llm)
  REAL :: phike(klon,llm,llm)
  REAL :: updke(klon,llm),dndke(klon,llm)
  REAL :: wghtke(klon,llm)

  !Thermiques
  REAL :: entr_therm(klon,llm)
  REAL :: fm_therm(klon,llm)

  REAL :: t(klon,llm)
  
  REAL :: aready(iip1,jjp1)
  REAL :: coefkz(klon,llm)
  REAL :: yu1(klon), yv1(klon)
  REAL :: ftsol(klon,nbsrf),pctsrf(klon,nbsrf)
  
  REAL :: massefi(klon,llm)
  
  !    Flux masse
  REAL :: massem(iip1,jjp1,llm),tetan(iip1,jjp1,llm)
  REAL :: pbarun(iip1,jjp1,llm),pbarvn(iip1,jjm,llm)
  REAL :: wn(iip1,jjp1,llm),phin(iip1,jjp1,llm)
  REAL :: phis(iip1,jjp1)
  
  REAL :: areafi(klon),phisfi(klon)
  
  REAL :: zdtvr, ziadvtrac, ziadvtrac2
  REAL :: zcontrole(klon),zmass,zflux
  
  INTEGER :: zklon,zklev,zrec
  INTEGER,SAVE :: zim,zjm,zlm
  
  REAL :: zpi
  
  zpi=2.*ASIN(1.)
      
      
!==================================================================
!   Si le numero du record est 0 alors: INITIALISATION
!==================================================================
!
  
  IF(irec==0) THEN
      
      !==================================================================
      !   ouverture des fichiers
      !==================================================================
      
      CALL read_dstoke_p(0,zdtvr,ziadvtrac,ziadvtrac2)

      CALL read_fstoke0_p(0, zrec,zim,zjm,zlm, aready,phis, massem,pbarun,pbarvn,wn,tetan,phin)
      ! 
      IF(physic)THEN
          CALL VTb(VTread_pstoke0)
          CALL read_pstoke0(0, &
             zrec,zklon,zklev,areafi,phisfi, &
             t,mfu,mfd,de_u,en_d,coefkz, &
             dake,mpke,phike,updke,dndke,wghtke,&    !variables schema KE
             entr_therm,fm_therm,&                   !variables thermiques
             yu1,yv1,ftsol,pctsrf)
          CALL VTe(VTread_pstoke0)
      ENDIF
      
      nrec=zrec
      istdyn=ziadvtrac
      istphy=ziadvtrac2
      dtvr  =zdtvr
          
      !==================================================================
      !   Fin des initialisations
  ELSE                      ! irec=0
      !==================================================================
          
          
      !-----------------------------------------------------------------------
      !   Lecture des fichiers fluxmass et  physique:
      !   -----------------------------------------------------
      
      CALL read_fstoke0_p(irec, zrec,zim,zjm,zlm, aready,phis, massem,pbarun,pbarvn,wn,tetan,phin)
      
      IF(physic)THEN
          
          CALL VTb(VTread_pstoke0)
          CALL read_pstoke0(irec, &
             zrec,zklon,zklev,areafi,phisfi, &
             t,mfu,mfd,de_u,en_d,coefkz, &
             dake,mpke,phike,updke,dndke,wghtke,&    !variables schema KE
             entr_therm,fm_therm,&                   !variables thermiques
             yu1,yv1,ftsol,pctsrf)
          CALL VTe(VTread_pstoke0)
          
          DO l=1,llm
            DO ig=1,klon
              zmfu(ig,l)=mfu(ig,l)
              zmfd(ig,l)=mfd(ig,l)
              zde_u(ig,l)=de_u(ig,l)
              zen_d(ig,l)=en_d(ig,l)
            ENDDO
          ENDDO
          !-----------------------------------------------------------------------
          !   PETIT CONTROLE SUR LES FLUX CONVECTIFS...
          !-----------------------------------------------------------------------
              
          CALL gr_dyn_fi_p(llm,iip1,jjp1,klon,massem,massefi)
          
          DO ig=1,klon
            zcontrole(ig)=1.
          ENDDO
          DO l=2,llm
            DO ig=1,klon
              zmass=MAX(massefi(ig,l),massefi(ig,l-1))/areafi(ig)
              zflux=MAX(ABS(zmfu(ig,l)),ABS(zmfd(ig,l)))*dtphys
              IF(zflux>0.9*zmass) THEN
                  zcontrole(ig)=MIN(zcontrole(ig),0.9*zmass/zflux)
              ENDIF
            ENDDO
          ENDDO

!        do ig=1,klon
!           if(zcontrole(ig).lt.0.99999) then
!pb               print*,'ATTENTION !!! on reduit les flux de masse '
!pb               print*,'convectifs au point ig=',ig
!           endif
!        enddo
          
          DO l=1,llm
            DO ig=1,klon
              zmfu(ig,l)=zmfu(ig,l)*zcontrole(ig)
              zmfd(ig,l)=zmfd(ig,l)*zcontrole(ig)
              zde_u(ig,l)=zde_u(ig,l)*zcontrole(ig)
              zen_d(ig,l)=zen_d(ig,l)*zcontrole(ig)
            ENDDO
          ENDDO
      ENDIF
      
      
  ENDIF                     ! irec=0
  
  RETURN
END SUBROUTINE readfluxnc_p
