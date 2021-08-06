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

SUBROUTINE timeloop_ad(iday_step, itau, nsplit, nsplit_dyn, &
           nsplit_phy, dtdyn, p, paprs, pplay, q, q_ad, qfi, &
           eflux, qstotraj, qfistotraj, efluxstotraj,&
           w, coefkz, nsec, wfunc_int, &
           periodflux, nbounds, bounddays, sflux_ad, masse, &
           pbaru, pbarv, t, zmfu, zmfd, zde_u, zen_d, pente_max, nobs, tabobs, &
           dake,mpke,phike,updke,dndke,wghtke,&
           entr_therm,fm_therm, &
           nbsat, nlevavgmax,  infos,nbformul, ps, &
           prescrstotraj, prodstotraj, &
           refprescr, refprod, temp, pmid, refjrates, depvel_loc, &
           nblossmax, nbprodmax, iday, month,  &
           pscale_ad, scale_ad)


  USE parallel
  USE mod_hallo
  USE vampir
  USE dimphy
  USE mod_interface_dyn_phys
  USE mod_phys_lmdz_para, ONLY : klon_mpi_para_end,scatter_phy=>scatter
  USE times
  USE Write_Field
  USE Write_field_p
  USE SPECIES_NAME
  USE CONSTANTS
  USE VARIABLES
  IMPLICIT NONE
  
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "comconst.h"
  INCLUDE "comvert.h"
  INCLUDE "comgeom2.h"
  INCLUDE "logic.h"
  INCLUDE "temps.h"
  INCLUDE "control.h"
  INCLUDE "tracstoke.h"
  INCLUDE "ajout.h"
  
  INTEGER :: iii, iq, j, i, l, ii, jj
  INTEGER :: itau, iday_step, itausplit, itauobs
  INTEGER :: isplit,nsplit_dyn,nsplit,nsplit_phy
  REAL :: dtdyn
  REAL :: masse(iip1,jjp1,llm)
  REAL :: pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm) !flux de masse
  REAL, DIMENSION(iip1,jjp1,llm)     :: w
  REAL, DIMENSION(iip1,jjp1,llmp1)   :: p
  REAL :: t(klon,llm)
  REAL :: zmfd(klon,llm),zen_d(klon,llm)
  REAL :: zmfu(klon,llm),zde_u(klon,llm)
  REAL :: dake(klon,llm),mpke(klon,llm)
  REAL :: phike(klon,llm,llm)
  REAL :: updke(klon,llm),dndke(klon,llm),wghtke(klon,llm)
  REAL :: entr_therm(klon,llm),fm_therm(klon,llm+1)   !fm a une dim de plus
  REAL, DIMENSION(klon,llm)   ::  coefkz
  REAL, DIMENSION(llm,nsec,iip1,jjp1)  :: wfunc_int
  INTEGER :: nsec, days_done, jjb,jje, nobs
  INTEGER :: isec
  REAL, DIMENSION(iip1,jjp1,llm,iqmax)   :: q, q_ad
  INTEGER  :: i1, i2, periodflux
  INTEGER  :: nbounds
  REAL     :: t1, t2, coefffc
  INTEGER, DIMENSION(nday+1)  :: bounddays
  REAL, DIMENSION(iip1,jjp1,nbounds,iqmax)   :: sflux_ad
  REAL :: eflux(klon,iqmax)
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: eflux_ad
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: source_ad
  REAL :: pente_max
  REAL :: paprs(klon,llm+1),pplay(klon,llm)
  REAL :: qstotraj(iip1,jjp1,llm,iqmax,nsplit*nsplit_dyn)
  REAL :: qfistotraj(klon,llm,iqmax,nsplit*nsplit_phy)
  REAL :: efluxstotraj(klon,iqmax,nsplit*nsplit_phy)
  REAL :: qfi(klon,llm,iqmax), qfi_ad(klon,llm,iqmax)
  REAL :: mr_ad, ppmv_ad
  REAL, DIMENSION(llm) :: p_mid,dp,dpp
  REAL, DIMENSION(llm) :: veci, veci_ad
  REAL, ALLOCATABLE, DIMENSION(:)  :: vecf, vecf_ad
  REAL :: tabobs(nobs,10)
!gch
  INTEGER ::  nbsat, nlevavgmax(nbsat), nbformul(nbsat)
  REAL :: nlevavglu
  INTEGER ::  iq1, nlev, nlevavg,nblinfo, nbs, chosenlev
  INTEGER, DIMENSION(nbsat) :: infos
  REAL, ALLOCATABLE, DIMENSION(:) :: qintavg, qintavg_ad
  REAL, ALLOCATABLE, DIMENSION(:) :: amat, deltap, qprim
  REAL, ALLOCATABLE, DIMENSION(:) :: tmp_ad
  REAL, ALLOCATABLE, DIMENSION(:)  :: pavg,pavg_mid, qa0, ak
  REAL, ALLOCATABLE, DIMENSION(:)  ::  qa0lu,dpavg
  REAL, DIMENSION(klon,llm,iprescrmax) :: locprescr, locprescr_ad
  REAL, DIMENSION(klon,llm,iprodmax)   :: locprod, locprod_ad
  REAL, DIMENSION(klon,llm,nsplit*nsplit_phy,iprescrmax):: prescrstotraj
  REAL, DIMENSION(klon,llm,nsplit*nsplit_phy,iprodmax):: prodstotraj
  REAL :: mr_log, mr_log_ad
  REAL :: convpres, mair, conv, frac1, frac2
  REAL :: ps(iip1,jjp1)
  REAL, DIMENSION(klon,llm) :: temp, pmid
  REAL, DIMENSION(klon,llm,iprescrmax) :: refprescr, refprescr_ad
  REAL, DIMENSION(klon,llm,iprodmax)   :: refprod, refprod_ad
  REAL, DIMENSION(klon,llm,ijratesmax) :: refjrates
  REAL, DIMENSION(klon,idepmax)        :: depvel_loc
  REAL, DIMENSION(klon,iprodmax)       :: phypscale_ad
  REAL, DIMENSION(klon,iprescrmax)     :: physcale_ad
  REAL, DIMENSION(iip1,jjp1,iprodmax)  :: intpscale_ad
  REAL, DIMENSION(iip1,jjp1,iprescrmax):: intscale_ad
  REAL, DIMENSION(iip1,jjp1,nbounds,iprodmax):: pscale_ad
  REAL, DIMENSION(iip1,jjp1,nbounds,iprescrmax):: scale_ad
  INTEGER :: nblossmax, nbprodmax, iday, month, month1, month2, dec
  
!
  ALLOCATE( eflux_ad(klon,iqmax) )
  ALLOCATE( source_ad(iip1,jjp1,iqmax) )

  DO isplit=nsplit,1,-1
    ! backward for itau, forward for isplit...
    itausplit=(itau-1)*nsplit+nsplit-isplit+1 ! sub time step
    itauobs=nday * iday_step * nsplit - itausplit + 1
!    print*, 'itausplit', itausplit, itau, nsplit, isplit
    !     AD statement of comparison with obs
    IF ( any( tabobs(:,1) <= itauobs ) &
            .AND. any( tabobs(:,1) + tabobs(:,2) > itauobs )) THEN
      DO i = nobs,1,-1
        jj=tabobs(i,3)
        ii=tabobs(i,4)
        IF ( tabobs(i,1) <= itauobs &
                .AND. tabobs(i,1) + tabobs(i,2) > itauobs) THEN
!          PRINT*, 'iobs', i, 'itausplit',itausplit, 'itau', itau, 'itauobs', itauobs, tabobs(i,:)
          iq = abs ( int(tabobs(i,6) ) )
          IF (iq <= iqmax .and. iq > 0) THEN
            IF (tabobs(i,6) < 0 ) iq = abs (int( tabobs(i,6)/100. )) !nec. pr adv_mass
            ppmv_ad = tabobs(i,5)
            mr_ad = ppmv_ad / adv_mass(iq) *dry_mass *1.e6 !gch
  
            IF (tabobs(i,6) < 0 ) THEN ! satellite
              iq = abs (int( tabobs(i,6)/100. ))
              iq1 = int ( abs(tabobs(i,6))-iq*100)
              ALLOCATE(qprim(llm))
       
              ! lecture des pavg, ak et qa0 voulus
              ! decodage du numero de ligne voulu ds fic du mois
              nblinfo=int((int(tabobs(i,6))-tabobs(i,6))*10000000.)
              IF( nbformul(iq1)==1 .OR. nbformul(iq1)==2 ) THEN !  general formula (Warning diff 1-2 in mair and equation used)
                 READ(infos(iq1),rec=nblinfo) nlevavglu
                 nlevavg=int(nlevavglu)
                 ALLOCATE(pavg(nlevavg+1)) ! pressure levels incl. surface pressure
                 ALLOCATE(pavg_mid(nlevavg)) ! pressure at the middle of the layer
                 ALLOCATE(qa0lu(nlevavg)) ! prior profile
                 ALLOCATE(ak(nlevavg))  ! averaging kernels
                 read(infos(iq1),rec=nblinfo)nlevavglu,ak(:),qa0lu(:),pavg(:)
                 ! dpavg to compute correction of FC
                 ! pavg_mid = pressure of middle of layer, for intex
                 ! rk: surface = nlevavg+1
                 ! dpavg=level thickness
                 ALLOCATE(dpavg(nlevavg))
                 do nlev=1,nlevavg
                   dpavg(nlev)=pavg(nlev+1)-pavg(nlev)
                   pavg_mid(nlev)=(pavg(nlev+1)+pavg(nlev))/2.
                 enddo
                 !interpolate model concentrations on averaging kernel pressure grid
                 DO l = 1, llm
                   p_mid(llm+1-l) = (p(tabobs(i,4),tabobs(i,3),l)+p(tabobs(i,4),tabobs(i,3),l+1))/2.
                   dpp(llm+1-l) = p(tabobs(i,4),tabobs(i,3),l)-p(tabobs(i,4),tabobs(i,3),l+1)
                 ENDDO
                 nbs=60+mpi_rank
                 read(nbs,rec=i)qprim
                 ALLOCATE(vecf(nlevavg))
                 ALLOCATE(vecf_ad(nlevavg))
                 DO l = 1, llm
                   veci(llm+1-l) = qprim(l)
                 ENDDO
                 CALL intex (llm,nlevavg,p_mid,pavg_mid,veci,vecf)
                 ! correction by FC
                 coefffc = dot_product(veci,dpp) / dot_product(vecf,dpavg) &
                              / sum(dpp) * sum(dpavg)
                 vecf(:) = vecf(:) * coefffc
                 !PRINT*, 'coeff', coefffc,' obs numero',i ! veci, vecf
                 ! surface index 1 as in LMDz
                 ALLOCATE(qintavg(nlevavg))
                 ALLOCATE(qintavg_ad(nlevavg))
                 DO l = 1, nlevavg
                   qintavg(nlevavg+1-l) = vecf(l)
                 ENDDO
                 !AD model-equivalent of observation, in kg/kg
                 ! amat and qa0 so that surface index 1 + conversion
                 ALLOCATE(amat(nlevavg))
                 ALLOCATE(qa0(nlevavg))
                 ALLOCATE(deltap(nlevavg))
                 conv=adv_mass(iq)/dry_mass /1.e6 ! ppm -> vmr -> mmr
                 DO nlev=1,nlevavg
                   amat(nlev)=ak(nlevavg+1-nlev)
                   qa0(nlev)=qa0lu(nlevavg+1-nlev)*conv
                   deltap(nlev)=dpavg(nlevavg+1-nlev)
                 ENDDO
                 mair=0.
                 DO nlev=1,nlevavg
                   IF (nbformul(iq1)==1) mair=mair+deltap(nlev)
                   IF (nbformul(iq1)==2) mair=mair+amat(nlev)*deltap(nlev)
                 ENDDO
                 ! AD here ????????CONVERSION NEEDED????? WHERE WHEN WHY ???
                 mr_ad = mr_ad/mair
                 DO nlev=1,nlevavg
                   IF (nbformul(iq1)==1) qintavg_ad(nlev)=mr_ad*amat(nlev) * deltap(nlev)
                   IF (nbformul(iq1)==2) qintavg_ad(nlev)=mr_ad*amat(nlev) * deltap(nlev)
                 ENDDO
                 mr_ad=0.
  
                 !AD interpolate model concentrations on averaging kernel pressure grid
                 DO l = 1, nlevavg
                   vecf_ad(l) = qintavg_ad(nlevavg+1-l)
                   vecf_ad(l) = vecf_ad(l) * coefffc
                   qintavg_ad(nlevavg+1-l) = 0.
                 ENDDO
                 CALL intex_ad (llm,nlevavg,p_mid,pavg_mid,veci,vecf,veci_ad,vecf_ad)
                 DO l = 1, llm
                   q_ad(tabobs(i,4),tabobs(i,3),l,iq)=q_ad(tabobs(i,4),tabobs(i,3),l,iq) &
                                                      + veci_ad(llm+1-l)
                   veci_ad(llm+1-l) = 0.
                 ENDDO
  
                 DEALLOCATE(vecf)
                 DEALLOCATE(vecf_ad)
                 DEALLOCATE(qintavg)
                 DEALLOCATE(qintavg_ad)
                 DEALLOCATE(amat)
                 DEALLOCATE(deltap)
                 DEALLOCATE(qa0lu)
                 DEALLOCATE(qprim)
                 deallocate(ak)
                 deallocate(qa0)
                 deallocate(pavg)
                 deallocate(pavg_mid)
                 DEALLOCATE(dpavg)
              ENDIF ! satellite formula 1 and 2
            ELSE  ! surface
              l=nint((tabobs(i,6) - iq) *100)! extract level information
!              WRITE(*,*) 'EEEEEEEEEEE', trim(species(iq)), l, iq, q_ad(tabobs(i,3),tabobs(i,2),l,iq), mr_ad
!              WRITE(*,*)'OBS',i, tabobs(i,:)
!              WRITE(*,*)'OBS',i, mr_ad
!              WRITE(*,*)'OBS',i, q_ad(tabobs(i,3),tabobs(i,2),l,iq)
              q_ad(tabobs(i,4),tabobs(i,3),l,iq) = q_ad(tabobs(i,4),tabobs(i,3),l,iq) + mr_ad
!              WRITE(*,*)'OBS',i, q_ad(tabobs(i,3),tabobs(i,2),l,iq)
              
              
            ENDIF ! ! tabobs(i,5) < 0.
  
!            tabobs(i,5) = 0.
          END IF
        ENDIF  !tabobs(i,1) == itausplit
      ENDDO !i=nobs...
    ENDIF ! any( tabobs(:,1) == itausplit
    
    !
    !-----------------------------------------------------------------------
    !     PHYSICS
    !-----------------------------------------------------------------------
    !     ---------------------------------------------------------
    IF(physic) THEN

        ! Putting the tracer on the dynamics grid
        CALL gr_dyn_fi_p(llm*iqmax,iip1,jjp1,klon,q,qfi)
        
        jjb=jj_begin
        jje=jj_end
        IF (north_pole) jjb=jjb+1
        IF (south_pole) jje=jje-1
        
        DO j=jjb,jje
          q_ad(1,j,:,:)=q_ad(1,j,:,:)+q_ad(iip1,j,:,:)
        ENDDO
        CALL gr_fi_dyn_ad_p(llm*iqmax,iip1,jjp1,klon,q_ad,qfi_ad)
        
        !-----------------------------------------------------------------------c
        
        eflux_ad(:,:) = 0.
        physcale_ad(:,:)=0.
        phypscale_ad(:,:)=0.
        locprod_ad(:,:,:) = 0.
        locprescr_ad(:,:,:) = 0.

        DO iii=nsplit_phy,1,-1
          !     
          !     AD statement of physics
          !     
          
          eflux(1:klon,:) = efluxstotraj(1:klon,:,(isplit-1)*nsplit_phy+iii)
          qfi(1:klon,:,:) = qfistotraj(1:klon,:,:,(isplit-1)*nsplit_phy+iii)
          locprescr(:,:,:) = prescrstotraj(:,:,(isplit-1)*nsplit_phy+iii,:)
          locprod(:,:,:) = prodstotraj(:,:,(isplit-1)*nsplit_phy+iii,:)
          
          !                    CALL VTb(VTphysiq)
          CALL phytrac_ad(dtphys, t, paprs, pplay, &
                          zmfu, zmfd, zde_u, zen_d, coefkz, qfi, eflux, eflux_ad, qfi_ad, &
                          dake,phike,mpke,updke,dndke,wghtke,&
                          entr_therm,fm_therm, &
                          locprod,locprescr, temp, pmid,refjrates,depvel_loc, &
                          locprod_ad, locprescr_ad, nblossmax, nbprodmax)
          !                    CALL VTe(VTphysiq)
          
        ENDDO         ! nsplit_phy
        !
        ! Translating the tracer to the dynamics grid
        CALL gr_dyn_fi_ad_p(llm*iqmax,iip1,jjp1,klon,q_ad,qfi_ad)
        qfi_ad(1:klon,:,:)=0.
        
    ENDIF             ! if (physic)
    !
    !     AD statement of advection
    !     
    
    DO iii=nsplit_dyn,1,-1
      q(:,jj_begin:jj_end,:,:) = qstotraj(:,jj_begin:jj_end,:,:,(isplit-1)*nsplit_dyn+iii)
      DO iq=1,iqmax
        CALL VTb(VTadvection_ad)
        CALL vlsplt_ad_p(q(1,1,1,iq),pente_max,masse,w,pbaru,pbarv,dtdyn,q_ad(1,1,1,iq))
        CALL VTe(VTadvection_ad)
      END DO
    ENDDO
    ! AD statement of introduce source
    ! careful with time indices in AD!!!
    ! SACS
    i = iday-day_ini
    IF (i <= 15) THEN
      month1 = month-1
      month2 = month
      if (month1 < 1) month1 = 12
      frac1 = (15-i)/30.
      frac2 = 1. - frac1
    ELSE
      month1 = month
      month2 = month+1
      if (month2 > 12) month2 = 1
      frac2 = (i-15)/30.
      frac1 = 1. - frac2
    ENDIF
    dec=0 ! decalage pour les flx_ch2o qui sont globaux!!!
    IF (mpi_rank>0) dec=klon_mpi_para_end(mpi_rank-1)
!
    days_done = (itau-1) / iday_step
    j = MIN( (nday-days_done-1)/periodflux+1, nbounds )
    t1 = bounddays(j)*nsplit*iday_step+1.
    t2 = bounddays(j+1)*nsplit*iday_step+1.
    ! 
    ! SACS
    i1 = j
    i2 = j+1
    DO i = 1, klon
      DO l=1,llm
        phypscale_ad(i, :) = phypscale_ad(i, :) + locprod_ad(i,l,:)*refprod(i,l,:)
        physcale_ad(i, :)  = physcale_ad(i, :) + locprescr_ad(i,l,:)*refprescr(i,l,:)
      ENDDO
    ENDDO
    locprod_ad(:,:,:) = 0.
    locprescr_ad(:,:,:) = 0.
    
    IF(iprodmax > 0)call gr_dyn_fi_ad_p(iprodmax,iip1,jjp1,klon,intpscale_ad,phypscale_ad)
    IF(iprescrmax > 0)call gr_dyn_fi_ad_p(iprescrmax,iip1,jjp1,klon,intscale_ad,physcale_ad)
    physcale_ad(:,:)=0.
    phypscale_ad(:,:)=0.
    
    call gr_dyn_fi_ad_p(iqmax,iip1,jjp1,klon,source_ad,eflux_ad)
    eflux_ad(:,:)=0.
    
    DO i = 1, iip1
      IF(iprodmax > 0) THEN
        DO iq = 1, iprodmax
          pscale_ad(i,jj_begin:jj_end,i1,iq) = pscale_ad(i,jj_begin:jj_end,i1,iq) + &
                  intpscale_ad(i,jj_begin:jj_end,iq)
!                  (1-(nday*nsplit*iday_step+1-itausplit - t1)/(t2 - t1)) &
!                  *intpscale_ad(i,jj_begin:jj_end,iq)
!          pscale_ad(i,jj_begin:jj_end,i2,iq) = pscale_ad(i,jj_begin:jj_end,i2,iq) + &
!                 (nday*nsplit*iday_step+1-itausplit - t1)/(t2 - t1) &
!                  *intpscale_ad(i,jj_begin:jj_end,iq)
        END DO
      END IF
      
      IF(iprescrmax > 0) THEN
        DO iq = 1, iprescrmax
          scale_ad(i,jj_begin:jj_end,i1,iq) = scale_ad(i,jj_begin:jj_end,i1,iq) + &
                 intscale_ad(i,jj_begin:jj_end,iq)
!                 (1-(nday*nsplit*iday_step+1-itausplit - t1)/(t2 - t1)) &
!                  *intscale_ad(i,jj_begin:jj_end,iq)
!          scale_ad(i,jj_begin:jj_end,i2,iq) = scale_ad(i,jj_begin:jj_end,i2,iq) + &
!                 (nday*nsplit*iday_step+1-itausplit - t1)/(t2 - t1) &
!                  *intscale_ad(i,jj_begin:jj_end,iq)
        END DO
      END IF
    
      DO iq = 1, iqmax
        sflux_ad(i,jj_begin:jj_end,i1,iq) = sflux_ad(i,jj_begin:jj_end,i1,iq) + &
                source_ad(i,jj_begin:jj_end,iq)
!            (1-(nday*nsplit*iday_step+1-itausplit - t1)/(t2 - t1)) &
!             *source_ad(i,jj_begin:jj_end,iq)
!        sflux_ad(i,jj_begin:jj_end,i2,iq) = sflux_ad(i,jj_begin:jj_end,i2,iq) + &
!            (nday*nsplit*iday_step+1-itausplit - t1)/(t2 - t1) &
!             *source_ad(i,jj_begin:jj_end,iq)
      ENDDO
    ENDDO
!
    !     
  ENDDO               ! isplit=nsplit,1,-1

  DEALLOCATE( eflux_ad )
  DEALLOCATE( source_ad )

RETURN
END SUBROUTINE timeloop_ad
