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

SUBROUTINE phytrac (pdtphys, t_seri, paprs, pplay, areafi, &
               pmfu, pmfd, pde_u, pen_d, coefh, tr_seri, eflux, &
               entr_therm,fm_therm,da,phi,mp,upd,dnd,wght, &
               refprod,refprescr,temp, pmid,refjrates,depvel_loc, &
               sumd0,sump0,suml0,sume0,nbprodmax,nblossmax)
  USE dimphy
  USE SPECIES_NAME
  USE CONSTANTS
  IMPLICIT NONE
  
  !======================================================================
  ! Auteur(s) FH
  ! Objet: Moniteur general des tendances traceurs
  !======================================================================
  INCLUDE "YOMCST.h"
  INCLUDE "dimensions.h"
  INCLUDE "indicesol.h"
  INCLUDE "temps.h"
  INCLUDE "control.h"
  INCLUDE "paramet.h"
  INCLUDE "clesph0.h"
  !======================================================================
  
  ! Arguments:
  !
  !   EN ENTREE:
  !   ==========
  !
  !   divers:
  !   -------
  !
  ! integer :: iqmax ! nombre de traceurs auxquels on applique la physique
  real :: pdtphys  ! pas d'integration pour la physique (seconde)
  REAL :: t_seri(klon,llm) ! temperature
  REAL :: tr_seri(klon,llm,iqmax) ! traceur  
  REAL :: paprs(klon,llm+1)  ! pression pour chaque inter-couche (en Pa)
  REAL :: pplay(klon,llm)  ! pression pour le mileu de chaque couche (en Pa)
  !
  !   convection:
  !   -----------
  !
  REAL :: pmfu(klon,llm)  ! flux de masse dans le panache montant
  REAL :: pmfd(klon,llm)  ! flux de masse dans le panache descendant
  REAL :: pde_u(klon,llm) ! flux detraine dans le panache montant
  REAL :: pen_d(klon,llm) ! flux entraine dans le panache descendant
 
!Kerry Emanuel
  REAL :: da(klon,llm)  ! taux de detrainement de lascendance adiabatique
  REAL :: phi(klon,llm,llm)  ! flux de masse melange du a lenvironnement
  REAL :: mp(klon,llm) ! flux de masse dans la descente insaturee
  REAL :: upd(klon,llm) ! saturated updraft mass flux
  REAL :: dnd(klon,llm) ! saturated downdraft mass flux
  REAL :: wght(klon,llm) ! poids des couches pour la couche d alimentation
  !
  !   Thermiques:
  !   -----------
  !
  REAL :: fm_therm(klon,llm+1)
  REAL :: entr_therm(klon,llm)
  !
  !   Couche limite:
  !   --------------
  !
  REAL :: coefh(klon,llm) ! coeff melange CL
  
  !
  ! Sources and sinks:
  ! ------------------------------
  !
  REAL :: source(klon)   
  REAL :: eflux(klon,iqmax)

  !
  INTEGER :: k, it, i, l, iq
  REAL :: delp(klon,llm)
  REAL :: zmasse(klon,llm) !(thermals)

  !Chemistry
  REAL, DIMENSION(klon,llm,iprescrmax) :: refprescr
  REAL, DIMENSION(klon,llm,iprodmax) :: refprod
  REAL, DIMENSION(klon,llm,ijratesmax) :: refjrates
  REAL, DIMENSION(klon,llm) :: temp,pmid
  REAL, DIMENSION(klon,idepmax) :: depvel_loc
  REAL, DIMENSION(klon)   :: areafi
  REAL :: boltzman
  REAL :: d_chem(klon,klev,iqmax) ! Trend for tracers
  INTEGER :: nblossmax,nbprodmax
  REAL :: sumd0(idepmax,klon),sume0(iqmax,klon)
  REAL :: suml0(llm-1,iqmax,nblossmax,klon),sump0(llm-1,iqmax,nbprodmax,klon)
  !
  !
  ! Variables locales pour effectuer les appels en serie
  !----------------------------------------------------
  !
  REAL :: d_tr(klon,llm) ! tendances de traceurs 
  REAL :: d_tr_cv(klon,llm) ! tendance de traceurs  convection 
  REAL :: d_tr_th(klon,llm) ! tendance de traceurs thermique
  !
  !   Controls
  !-------------
  LOGICAL,SAVE :: pblayer,convection
  DATA pblayer,convection /.TRUE.,.TRUE./

 !
  !======================================================================
  
  boltzman=1.3806e-23 ! JK-1
  sumd0=-99.

!  WRITE(*,*) 'Phytrac in timeloop !'


  DO k = 1, llm
    delp(:,k) = paprs(:,k)-paprs(:,k+1)
  ENDDO
  zmasse(:,:) = delp(:,:) / rg

  IF (convection) THEN
    DO it=1, iqmax
      if (iflag_con .eq. 2) then
       CALL nflxtr(pdtphys, pmfu, pmfd, pde_u, pen_d, &
         paprs, tr_seri(1,1,it), d_tr_cv)
      else if ((iflag_con .eq. 3).or.(iflag_con .eq. 30)) then
       CALL cvltr(pdtphys, da, phi, mp, paprs, &
         tr_seri(:,:,it), upd, dnd, wght, d_tr_cv)
      end if
      tr_seri(:,:,it) = tr_seri(:,:,it) + d_tr_cv
      IF (iflag_con .ge. 3) then
        DO k = 1, llm
        DO i = 1, klon
          tr_seri(i,k,it)=MAX(tr_seri(i,k,it),0.)
          tr_seri(i,k,it)=MIN(tr_seri(i,k,it),1.e10)
        ENDDO
        ENDDO
      ENDIF
      IF (iflag_con .eq. 30) then
        CALL thermcell_dq(klon,llm,pdtphys, &
        fm_therm,entr_therm,zmasse,tr_seri(:,:,it),d_tr_th)
        d_tr_th(:,:)=pdtphys*d_tr_th(:,:)
        tr_seri(:,:,it) = tr_seri(:,:,it)+d_tr_th(:,:)
      ENDIF

    ENDDO !iq
  ENDIF ! convection
  
  IF (pblayer) THEN
      DO it=1, iqmax
        sume0(it,:)=0.
        source(:) = eflux(:,it)
        sume0(it,:)=sume0(it,:)+source(:)*areafi(:)*pdtphys
        
        ! deposition
        DO iq=1, idepmax
          IF (dep_species(iq)%name == species(it)%name) THEN
            ! Conversion kg/kg ---> kg/cm3 : velocity in cm/s
            source(:) = source(:) - depvel_loc(:,iq) * tr_seri(:,1,it) &
                * 1.e-6  * (dry_mass*1.66e-23) * pmid(:,1) / (temp(:,1) * boltzman)
            sumd0(iq,:) = depvel_loc(:,iq) * tr_seri(:,1,it) &
                * 1.e-6 * (dry_mass*1.66e-23) * pmid(:,1) / (temp(:,1) * boltzman)  &
                * pdtphys * areafi(:)
          END IF
        END DO
        
        CALL cltrac(pdtphys, coefh,t_seri, tr_seri(1,1,it), source, paprs, pplay, delp, d_tr )
        tr_seri(:,:,it) = tr_seri(:,:,it) + d_tr(:,:)
      ENDDO
  ENDIF ! couche limite

  IF (do_chemistry) THEN
    CALL comp_chemistry(tr_seri, refprod, refprescr, refjrates, temp, pmid, pdtphys, &
                        d_chem, sump0, suml0, nbprodmax, nblossmax)
    tr_seri(:,:,:) = tr_seri(:,:,:) + d_chem(:,:,:)
    
    !!!!!!!!! WARNING! REMOVING NEGATIVE CONCENTRATIONS !!!!
!    do i=1,klon
!      do l=1,llm
!        do iq=1,iqmax
!          if(tr_seri(i,l,iq) .le. 0.)tr_seri(i,l,iq)=1.e-20
!        enddo
!      enddo
!    enddo
  ENDIF ! chemistry

!  WRITE(*,*) 'END Phytrac in timeloop !'

  
  RETURN
END SUBROUTINE phytrac
!---------------------------------------------------------------------
!---------------------------------------------------------------------
SUBROUTINE phytrac_tl ( pdtphys, t_seri,paprs,pplay, &
                        pmfu, pmfd, pde_u, pen_d, &
                        coefh, tr_seri, eflux, eflux_tl, tr_seri_tl, &
                        da,phi,mp,upd,dnd,wght,entr_therm,fm_therm, &
                        refprod,refprescr,temp,pmid,refjrates,depvel_loc, &
                        refprod_tl,refprescr_tl)
  USE parallel
  USE dimphy
  USE Write_Field_phy
  USE SPECIES_NAME
  USE CONSTANTS
  IMPLICIT NONE
  
  
  !======================================================================
  ! Auteur(s) FH
  ! Objet: Moniteur general des tendances traceurs
  !======================================================================
  INCLUDE "YOMCST.h"
  INCLUDE "dimensions.h"
  INCLUDE "indicesol.h"
  INCLUDE "temps.h"
  INCLUDE "control.h"
  INCLUDE "paramet.h"
 include "clesph0.h"
  !======================================================================
  
  ! Arguments:
  !
  !   EN ENTREE:
  !   ==========
  !
  !   divers:
  !   -------
  !
  ! INTEGER :: iqmax ! nombre de traceurs auxquels on applique la physique
  REAL :: pdtphys  ! pas d'integration pour la physique (seconde)
  REAL :: t_seri(klon,llm) ! temperature
  REAL :: tr_seri(klon,llm,iqmax) ! traceur  
  REAL :: tr_seri_tl(klon,llm,iqmax) ! traceur  
  REAL :: paprs(klon,llm+1)  ! pression pour chaque inter-couche (en Pa)
  REAL :: pplay(klon,llm)  ! pression pour le mileu de chaque couche (en Pa)
  !
  !   convection:
  !   -----------
  !
  REAL :: pmfu(klon,llm)  ! flux de masse dans le panache montant
  REAL :: pmfd(klon,llm)  ! flux de masse dans le panache descendant
  REAL :: pde_u(klon,llm) ! flux detraine dans le panache montant
  REAL :: pen_d(klon,llm) ! flux entraine dans le panache descendant
 ! Kerry Emanuel
  REAL :: da(klon,llm)  ! taux de detrainement de lascendance adiabatique
  REAL :: phi(klon,llm,llm)  ! flux de masse melange du a lenvironnement
  REAL :: mp(klon,llm) ! flux de masse dans la descente insaturee
  REAL :: upd(klon,llm) ! saturated updraft mass flux
  REAL :: dnd(klon,llm) ! saturated downdraft mass flux
  REAL :: wght(klon,llm) !poids des couches pour la couche dalimentation
 ! Thermiques
  REAL :: entr_therm(klon,llm)
  REAL :: fm_therm(klon,llm+1)
  !
  !
  !   Couche limite:
  !   --------------
  !
  REAL :: coefh(klon,llm) ! coeff melange CL
  !
  ! Sources et puits des traceurs:
  ! ------------------------------
  !
  REAL :: source(klon)       
  REAL :: source_tl(klon)      
  REAL :: eflux(klon,iqmax)   
  REAL :: eflux_tl(klon,iqmax)   
  
  !Chemistry
  REAL, DIMENSION(klon,llm,iprescrmax) :: refprescr, refprescr_tl
  REAL, DIMENSION(klon,llm,iprodmax) :: refprod, refprod_tl
  REAL, DIMENSION(klon,llm,ijratesmax) :: refjrates
  REAL, DIMENSION(klon,llm) :: temp,pmid
  REAL, DIMENSION(klon,idepmax) :: depvel_loc
  REAL :: boltzman
  REAL :: d_chem(klon,klev,iqmax) ! tendance de traceurs  chimie 
  REAL :: d_chem_tl(klon,klev,iqmax) ! tendance de traceurs  chimie
!
  
  !======================================================================
  !
  INTEGER :: k, it,i,l,iq
  REAL :: delp(klon,llm)
  REAL :: zmasse(klon,llm)
  !
  ! Variables locales pour effectuer les appels en serie
  !----------------------------------------------------
  !
  REAL :: d_tr_cl(klon,llm) ! tendances de traceurs 
  REAL :: d_tr_cv(klon,llm) ! tendance de traceurs  convection 
  REAL :: d_tr_cl_tl(klon,llm) ! tendances de traceurs 
  REAL :: d_tr_cv_tl(klon,llm) ! tendance de traceurs  convection 
  REAL :: d_tr_th_tl(klon,llm)
  REAL :: d_tr_th(klon,llm)
  !
  !   Controles
  !-------------
  LOGICAL,SAVE :: pblayer,convection
  DATA pblayer,convection /.TRUE.,.TRUE./

  !
  !======================================================================

  boltzman=1.3806e-23 ! JK-1
  DO k = 1, llm
    delp(:,k) = paprs(:,k)-paprs(:,k+1)
  ENDDO
  zmasse(:,:) = delp(:,:) / rg

  IF (convection) THEN
    DO it=1, iqmax
      if (iflag_con .eq. 2) then
        CALL nflxtr_tl(pdtphys, pmfu, pmfd, pde_u, pen_d, &
           paprs, tr_seri(1,1,it), d_tr_cv, &
           tr_seri_tl(1,1,it), d_tr_cv_tl)
      else if ((iflag_con .eq. 3).or.(iflag_con .eq. 30)) then
        CALL cvltr_tl(pdtphys,da,phi,mp,paprs,&
         tr_seri(:,:,it),tr_seri_tl(:,:,it),upd,dnd,&
         wght,d_tr_cv,d_tr_cv_tl)
      endif
      !
      tr_seri_tl(:,:,it) = tr_seri_tl(:,:,it) + d_tr_cv_tl(:,:)
      tr_seri(:,:,it) = tr_seri(:,:,it) + d_tr_cv(:,:)
      DO k = 1,llm
        DO i = 1,klon
          if (tr_seri(i,k,it) .lt. 0.) tr_seri_tl(i,k,it) = 0.
          if (tr_seri(i,k,it) .gt. 1.e10) tr_seri_tl(i,k,it) = 0.
          tr_seri(i,k,it) = MAX(tr_seri(i,k,it),0.)
          tr_seri(i,k,it) = MIN(tr_seri(i,k,it),1.e10)
        END DO
      END DO
      if (iflag_con .eq. 30) then
        CALL thermcell_dq_tl(klon,llm,pdtphys,fm_therm,entr_therm,&
                          zmasse,tr_seri(1:klon,1:llm,it),&
                          tr_seri_tl(1:klon,1:llm,it),d_tr_th,&
                          d_tr_th_tl)
        tr_seri_tl(:,:,it)=tr_seri_tl(:,:,it)+d_tr_th_tl(:,:)*pdtphys
        tr_seri(:,:,it)=tr_seri(:,:,it)+d_tr_th(:,:)*pdtphys
      endif
    ENDDO
  ENDIF
!   CALL writefield_phy('tr_seri_apres',tr_seri(:,:,1),llm)
    
  IF (pblayer) THEN
    DO it=1, iqmax
      source_tl(:) = eflux_tl(:,it)
      source(:) = eflux(:,it)
      
      ! deposition
      DO iq=1, idepmax
        IF (dep_species(iq)%name == species(it)%name) THEN
          source_tl(:)=source_tl(:) - depvel_loc(:,iq) * tr_seri_tl(:,1,it) &
              * 1.e-6 * (dry_mass*1.66e-23) * pmid(:,1) / (temp(:,1) * boltzman)
          source(:)=source(:) - depvel_loc(:,iq) * tr_seri(:,1,it) &
              * 1.e-6 * (dry_mass*1.66e-23) * pmid(:,1) / (temp(:,1) * boltzman)
        END IF
      END DO
      
      CALL cltrac_tl(pdtphys, coefh,t_seri, &
           tr_seri(1,1,it), source, &
           paprs, pplay, delp, &
           d_tr_cl, source_tl, tr_seri_tl(1,1,it),d_tr_cl_tl )
      tr_seri_tl(:,:,it) = tr_seri_tl(:,:,it) + d_tr_cl_tl(:,:)
      tr_seri(:,:,it) = tr_seri(:,:,it) + d_tr_cl(:,:)
    ENDDO
  ENDIF ! couche limite
    
  IF (do_chemistry) then
    CALL comp_chemistry_tl(tr_seri, refprod, refprescr, refjrates, temp, pmid, tr_seri_tl, &
                        refprod_tl, refprescr_tl, pdtphys,&
                        d_chem, d_chem_tl)
    tr_seri_tl(:,:,:) = tr_seri_tl(:,:,:) + d_chem_tl(:,:,:)
    tr_seri(:,:,:) = tr_seri(:,:,:) + d_chem(:,:,:)

    !!!!!!!!! WARNING, REMOVING NEGATIVE OR NULL VALUES !!!!
!    do i=1,klon
!      do l=1,llm
!        do iq=1,iqmax !!! iq=2 pour le CO
!          !if(tr_seri(i,l,iq) .le. 0.)tr_seri(i,l,iq)=1.e-10
!          !if(tr_seri(i,l,iq) .le. 0.)tr_seri_tl(i,l,iq)=0.0
!          IF(tr_seri(i,l,iq) .LE. 0.)tr_seri(i,l,iq)=1e-20
!        enddo
!      enddo
!    enddo
  ENDIF ! chemistry
  
  RETURN
END SUBROUTINE phytrac_tl
!---------------------------------------------------------------------
!---------------------------------------------------------------------
SUBROUTINE phytrac_ad ( pdtphys, t_seri,paprs,pplay, &
                        pmfu, pmfd, pde_u, pen_d, &
                        coefh, tr_seri, eflux, eflux_ad, tr_seri_ad, &
                        da,phi,mp,upd,dnd,wght,entr_therm,fm_therm, &
                        refprod,refprescr,temp,pmid,refjrates,depvel_loc, &
                        refprod_ad,refprescr_ad,nblossmax,nbprodmax)
  USE dimphy
  USE SPECIES_NAME
  USE CONSTANTS
  IMPLICIT NONE
  
  
  !======================================================================
  ! Auteur(s) FH
  ! Objet: Moniteur general des tendances traceurs
  !======================================================================
  INCLUDE "YOMCST.h"
  INCLUDE "dimensions.h"
  INCLUDE "indicesol.h"
  INCLUDE "temps.h"
  INCLUDE "control.h"
  INCLUDE "paramet.h"
 include "clesph0.h"
  !======================================================================
  
  ! Arguments:
  !
  !   EN ENTREE:
  !   ==========
  !
  !   divers:
  !   -------
  !
  ! INTEGER :: iqmax ! nombre de traceurs auxquels on applique la physique
  REAL :: pdtphys  ! pas d'integration pour la physique (seconde)
  REAL :: t_seri(klon,llm) ! temperature
  REAL :: tr_seri(klon,llm,iqmax) ! traceur  
  REAL :: tr_seri_ad(klon,llm,iqmax) ! traceur  
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tr_seri1, tr_seri2,tr_seri3
  REAL :: paprs(klon,llm+1)  ! pression pour chaque inter-couche (en Pa)
  REAL :: pplay(klon,llm)  ! pression pour le mileu de chaque couche (en Pa)
  !
  !   convection:
  !   -----------
  !
  REAL :: pmfu(klon,llm)  ! flux de masse dans le panache montant
  REAL :: pmfd(klon,llm)  ! flux de masse dans le panache descendant
  REAL :: pde_u(klon,llm) ! flux detraine dans le panache montant
  REAL :: pen_d(klon,llm) ! flux entraine dans le panache descendant

   ! Kerry Emanuel
  REAL :: da(klon,llm)  ! taux de detrainement de lascendance adiabatique
  REAL :: phi(klon,llm,llm)  ! flux de masse melange du a lenvironnement
  REAL :: mp(klon,llm) ! flux de masse dans la descente insaturee
  REAL :: upd(klon,llm) ! saturated updraft mass flux
  REAL :: dnd(klon,llm) ! saturated downdraft mass flux
  REAL :: wght(klon,llm) !poids des couches pour la couche dalimentation
 ! Thermiques
  REAL :: entr_therm(klon,llm)
  REAL :: fm_therm(klon,llm+1)
  
  !   Couche limite:
  !   --------------
  !
  REAL :: coefh(klon,llm) ! coeff melange CL
  !
  ! Sources et puits des traceurs:
  ! ------------------------------
  !
  REAL :: source(klon)       
  REAL :: source_ad(klon)   
  REAL :: eflux(klon,iqmax)
  REAL :: eflux_ad(klon,iqmax)
  
  !Chemistry
  REAL, DIMENSION(klon,llm,iprescrmax) :: refprescr, refprescr_ad
  REAL, DIMENSION(klon,llm,iprodmax) :: refprod, refprod_ad
  REAL, DIMENSION(klon,llm,ijratesmax) :: refjrates
  REAL, DIMENSION(klon,llm) :: temp,pmid
  REAL, DIMENSION(klon,idepmax) :: depvel_loc
  REAL :: boltzman
  REAL :: d_chem(klon,klev,iqmax) ! tendance de traceurs  chimie 
  REAL :: d_chem_ad(klon,klev,iqmax) ! tendance de traceurs  chimie
  INTEGER :: nblossmax,nbprodmax
  REAL :: sump1(klev-1,iqmax,nbprodmax,klon),suml1(klev-1,iqmax,nblossmax,klon)

!
  
  !======================================================================
  !
  ! Declaration des procedures appelees
  !
  INTEGER :: k, it,i,l,iq
  REAL :: delp(klon,llm)
  REAL :: zmasse(klon,llm)
  !
  ! Variables locales pour effectuer les appels en serie
  !----------------------------------------------------
  !
  REAL, ALLOCATABLE, DIMENSION(:,:) :: d_tr_cl, d_tr_cv, d_tr_cl_ad, d_tr_cv_ad
  REAL, ALLOCATABLE, DIMENSION(:,:) :: d_tr_th, d_tr_th_ad
  !
  !   Controles
  !-------------
  LOGICAL,SAVE :: pblayer,convection
  DATA pblayer,convection /.TRUE.,.TRUE./


  !
  !======================================================================
  
  ALLOCATE( tr_seri1(klon,llm,iqmax) )
  ALLOCATE( tr_seri2(klon,llm,iqmax) )
  ALLOCATE( tr_seri3(klon,llm,iqmax) )
  ALLOCATE( d_tr_cl(klon,llm) )
  ALLOCATE( d_tr_cv(klon,llm) )
  ALLOCATE( d_tr_cl_ad(klon,llm) )
  ALLOCATE( d_tr_cv_ad(klon,llm) )
  ALLOCATE( d_tr_th(klon,llm) )
  ALLOCATE( d_tr_th_ad(klon,llm) )

  boltzman=1.3806e-23 ! JK-1
  DO k = 1, llm
    delp(:,k) = paprs(:,k)-paprs(:,k+1)
  ENDDO
  zmasse(:,:) = delp(:,:) / rg
  tr_seri1 = tr_seri
  IF (convection) THEN
    DO it=1, iqmax
        IF (iflag_con .eq. 2) THEN
           CALL nflxtr(pdtphys, pmfu, pmfd, pde_u, pen_d, &
           paprs, tr_seri(1,1,it), d_tr_cv)
        ELSE IF ((iflag_con .eq. 30).or.(iflag_con .eq. 3)) then
          CALL cvltr(pdtphys, da, phi, mp, paprs, &
           tr_seri(:,:,it), upd, dnd, wght, d_tr_cv)
        ENDIF
        tr_seri(:,:,it) = tr_seri(:,:,it) + d_tr_cv(:,:)
    ENDDO
  ENDIF
  tr_seri2 = tr_seri
  DO it = 1,iqmax
    DO k = 1,llm
      DO i = 1,klon
        tr_seri(i,k,it)=MAX(tr_seri(i,k,it),0.)
        tr_seri(i,k,it)=MIN(tr_seri(i,k,it),1.e10)
      END DO
    END DO
  END DO

  tr_seri3 = tr_seri
  IF (iflag_con .eq. 30) then
   DO it=1,iqmax
    CALL thermcell_dq(klon,llm,pdtphys,fm_therm,entr_therm, &
       zmasse,tr_seri(:,:,it),d_tr_th)
    tr_seri(:,:,it)=tr_seri(:,:,it)+pdtphys*d_tr_th(:,:)
   ENDDO       !iqmax
  ENDIF
  IF (pblayer) THEN
    DO it=1, iqmax
      source(:) = eflux(:,it)
      
      ! deposition
      DO iq=1, idepmax
        IF (dep_species(iq)%name == species(it)%name) THEN
          source(:)=source(:)-depvel_loc(:,iq)*tr_seri(:,1,it) &
              * 1.e-6 * (dry_mass*1.66e-23) * pmid(:,1) / ( temp(:,1) * boltzman )
        END IF
      END DO
      
      CALL cltrac(pdtphys, coefh,t_seri, tr_seri(1,1,it), source, paprs, pplay, delp, d_tr_cl )
      tr_seri(:,:,it) = tr_seri(:,:,it) + d_tr_cl(:,:)
    ENDDO
  ENDIF

  !Start of AD
  IF (do_chemistry) THEN
    d_chem_ad(:,:,:) = tr_seri_ad(:,:,:)

    CALL comp_chemistry_ad(tr_seri, refprescr, refprod, refjrates, temp, pmid, pdtphys, &
                        d_chem_ad, tr_seri_ad, refprod_ad, refprescr_ad)

    CALL comp_chemistry(tr_seri, refprod, refprescr, refjrates, temp, pmid, pdtphys, &
                        d_chem, sump1, suml1, nbprodmax, nblossmax)
    tr_seri(:,:,:) = tr_seri(:,:,:) + d_chem(:,:,:)
    
    !!!!!!!!! WARNING: PREVENTING NEGATIVE CONCENTRATIONS!!!!
!    do i=1,klon
!      do l=1,llm
!        do iq=1,iqmax !!! iq=2 pour le CO
!        !if(tr_seri(i,l,iq) .le. 0.)tr_seri(i,l,iq)=1.e-10
!          IF(tr_seri(i,l,iq) <= 0.)tr_seri(i,l,iq)=1.e-20
!        enddo
!      enddo
!    enddo
  ENDIF !chemistry

  DO it=1, iqmax
    IF (pblayer) THEN
      
      d_tr_cl_ad(:,:) = tr_seri_ad(:,:,it)
      CALL cltrac_ad(pdtphys, coefh,t_seri, &
           tr_seri(1,1,it), source, &
           paprs, pplay, delp, &
           d_tr_cl,source_ad, tr_seri_ad(1,1,it),d_tr_cl_ad)
      
      ! deposition
      DO iq=1, idepmax
        IF (dep_species(iq)%name == species(it)%name) THEN
          tr_seri_ad(:,1,it)=tr_seri_ad(:,1,it)-depvel_loc(:,iq)*source_ad(:) &
            * 1.e-6 * (dry_mass*1.66e-23) * pmid(:,1) / (temp(:,1) * boltzman)
        END IF
      END DO
      
      eflux_ad(:,it) = eflux_ad(:,it) + source_ad(:)
      source_ad(:) = 0.
      tr_seri(:,:,it) = tr_seri(:,:,it) + d_tr_cl(:,:)
    
    ENDIF                   ! pblayer
    
    IF (convection) THEN
        IF (iflag_con .eq. 30) then
          d_tr_th_ad(:,:)=tr_seri_ad(:,:,it)*pdtphys
          CALL thermcell_dq_ad(klon,llm,pdtphys,fm_therm,entr_therm,&
          zmasse,tr_seri3(:,:,it),tr_seri_ad(:,:,it),d_tr_th,d_tr_th_ad)
          DO k = 1,llm
           DO i = 1,klon
            if (tr_seri2(i,k,it) .lt. 0.) tr_seri_ad(i,k,it)=0.
            if (tr_seri2(i,k,it) .gt. 1.e10) tr_seri_ad(i,k,it)=0.
           END DO
          END DO
        END IF ! iflagcon
        d_tr_cv_ad(:,:) = tr_seri_ad(:,:,it)
        IF (iflag_con .eq. 2) then
          CALL nflxtr_ad(pdtphys, pmfu, pmfd, pde_u, pen_d, &
          paprs, tr_seri1(1,1,it), d_tr_cv, &
          tr_seri_ad(1,1,it), d_tr_cv_ad)
        ELSE IF ((iflag_con .eq. 3).OR.(iflag_con .eq. 30)) THEN
          CALL cvltr_ad(pdtphys,da,phi,mp,paprs,tr_seri1(:,:,it),&
          tr_seri_ad(:,:,it),upd,dnd,wght,d_tr_cv,d_tr_cv_ad)
        ENDIF
    ENDIF ! convection
  ENDDO !it
  DEALLOCATE( tr_seri1 )
  DEALLOCATE( tr_seri2 )
  DEALLOCATE( tr_seri3 )
  DEALLOCATE( d_tr_cl )
  DEALLOCATE( d_tr_cv )
  DEALLOCATE( d_tr_cl_ad )
  DEALLOCATE( d_tr_cv_ad )
  DEALLOCATE( d_tr_th_ad )

  RETURN
END SUBROUTINE phytrac_ad
