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

SUBROUTINE timeloop(nobs_glo, time_0, iday_step, fwd, nsplit, nsplit_dyn, &
        nsplit_phy, day0, nrec, areafi, pphis, phis, q, q_ad, ps, &
        q_tl, phi, w, coefkz, tabobs_glo, nsec, &
        wfunc_int, indexq, itrajq, &
        periodflux, nbounds, bounddays, sflux_tl, sflux, sflux_ad, masse, &
        pbaru,pbarv,teta,t,zmfu, zmfd, zde_u,zen_d, yu1, yv1, ftsol, pctsrf, &
        vcov, ucov,tangent, &
        dake,mpke,phike,updke,dndke,wghtke, entr_therm, fm_therm, &
        refprescr_glo,refprod_glo,temp_glo,pmid_glo,refjrates_glo, &
        scale, scale_tl, scale_ad, pscale, pscale_tl, pscale_ad, &
        depvel_glo,month,year,idiagd,idiagl,idiagp,idiage,idiagb,outdiag,outwfunc,convoh)
  
  
  USE parallel
  USE mod_hallo
  USE vampir
  USE dimphy
  USE mod_interface_dyn_phys
  USE mod_phys_lmdz_para, ONLY : klon_mpi_para_nb,klon_mpi_para_end,scatter_phy=>scatter
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
  INCLUDE "mpif.h"
  
  INTEGER :: ig0,ii,jj,iii
  INTEGER :: ig,iq,j,i,l
  INTEGER*4 :: iday            ! Julian day
  REAL :: time_0               ! Time in the day for the beginning of the simulation (% day)
  INTEGER :: itau,irec,nrec
  INTEGER :: iday_step,day0
  LOGICAL :: fwd
  REAL :: dtdyn
  INTEGER :: isplit,nsplit_dyn,nsplit,nsplit_phy
  REAL :: masse(iip1,jjp1,llm)
  REAL :: pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm) ! Mass fluxes
  REAL, DIMENSION(iip1,jjp1,llm)     :: phi, w
  REAL :: teta(iip1,jjp1,llm)
  REAL :: areafi(klon)
  REAL :: pphis (klon)
  REAL :: t(klon,llm)
  REAL :: zmfd(klon,llm),zen_d(klon,llm)
  REAL :: zmfu(klon,llm),zde_u(klon,llm)
  REAL, DIMENSION(klon,llm)   ::  coefkz
  REAL :: dake(klon,llm),mpke(klon,llm)
  REAL :: phike(klon,llm,llm)
  REAL :: updke(klon,llm),dndke(klon,llm),wghtke(klon,llm)
  REAL :: fm_therm(klon,llm),entr_therm(klon,llm)
  REAL :: fm_thermp1(klon,llm+1)   !ajout d une dim pr dqthermcell

  REAL :: yu1(klon), yv1(klon)
  REAL :: ftsol(klon,nbsrf),pctsrf(klon,nbsrf)
  
  REAL :: phis(iip1,jjp1)      ! Geopotential at ground level
  REAL :: smass(iip1,jjp1)
  REAL :: ps(iip1,jjp1),pkf(ip1jmp1,llm)
  REAL :: p(iip1,jjp1,llmp1),pk(iip1,jjp1,llm)
  REAL :: alpha(ip1jmp1,llm),beta(ip1jmp1,llm)
  REAL :: pks(iip1,jjp1)       ! exner (f for filter)
  REAL, DIMENSION(llm,nsec,iip1,jjp1)  :: wfunc_int
  REAL, DIMENSION(llm)             :: p_mid, dp,  dpp
  INTEGER                          :: nsec, nobs_glo
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: tdyn
  INTEGER :: isec
  INTEGER :: days_done, indexq
  INTEGER :: itrajq(iqmax)
  INTEGER :: jjb,jje
  INTEGER :: nobs, nobs_para(0:mpi_size-1)
  REAL, DIMENSION(iip1,jjp1,llm,iqmax)   :: q, q_tl, q_ad
  INTEGER :: itausplit
  REAL, DIMENSION(nobs_glo,10):: tabobs_glo
  INTEGER                          :: i1, i2, periodflux, nbounds
  REAL                             :: t1, t2
  INTEGER, DIMENSION(nday+1)     :: bounddays
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: source, source_tl
  REAL, DIMENSION(iip1,jjp1,nbounds,iqmax)   :: sflux, sflux_tl, sflux_ad
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: eflux, eflux_tl
  REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  efluxstotraj
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: qstotraj
  LOGICAL :: tangent
  REAL, PARAMETER                  :: pente_max = 2.
  ! Van Leer I : 2.
  ! Godunov    : 0.
  REAL, ALLOCATABLE, DIMENSION(:,:) :: paprs, pplay
  REAL :: unskap,pksurcp
  REAL, DIMENSION(iip1,jjm,llm)     :: vcov
  REAL, DIMENSION(iip1,jjp1,llm)    :: ucov
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: qfistotraj
  REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: qfi, qfi_tl, qfi_ad
  REAL                             :: mr_tl, ppmv_tl
  REAL                             :: mr, ppmv
  REAL :: massebx(iip1,jjp1,llm),masseby(iip1,jjm,llm)
  TYPE(Request) :: mpi_request
  REAL, ALLOCATABLE :: tabobs(:,:)
  ! gch
  REAL                                    :: mr_log, mr_log_tl
  REAL, DIMENSION(ngridmx,llm,iprescrmax) :: refprescr_glo
  REAL, DIMENSION(ngridmx,llm,iprodmax)   :: refprod_glo
  REAL, DIMENSION(ngridmx,llm,ijratesmax)    :: refjrates_glo
  REAL, DIMENSION(klon,llm,iprescrmax)    :: refprescr, refprescr_tl
  REAL, DIMENSION(klon,llm,iprodmax)      :: refprod, refprod_tl
  REAL, DIMENSION(klon,llm,ijratesmax)    :: refjrates
  REAL, DIMENSION(klon,llm)               :: var_scatter
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: prescrstotraj, prodstotraj
  
  REAL, DIMENSION(ngridmx,llm) :: temp_glo, pmid_glo
  REAL, DIMENSION(ngridmx,idepmax) :: depvel_glo
  REAL, DIMENSION(klon,llm) :: temp,pmid
  REAL, DIMENSION(klon,idepmax) :: depvel_loc
  REAL, DIMENSION(iip1,jjp1,nbounds,iprescrmax) :: scale, scale_tl, scale_ad
  REAL, DIMENSION(iip1,jjp1,nbounds,iprodmax) :: pscale, pscale_tl, pscale_ad
  REAL, DIMENSION(iip1,jjp1,iprescrmax) :: intscale, intscale_tl !SACS
  REAL, DIMENSION(iip1,jjp1,iprodmax) :: intpscale, intpscale_tl !SACS
  REAL, DIMENSION(klon,iprescrmax) :: physcale, physcale_tl
  REAL, DIMENSION(klon,iprodmax) :: phypscale, phypscale_tl
  REAL, DIMENSION(ngridmx,iqmax):: eflux_glo
  REAL, DIMENSION(ngridmx,iprescrmax) :: physcale_glo
  REAL, DIMENSION(ngridmx,iprodmax) :: phypscale_glo
  INTEGER, DIMENSION(0:mpi_size-1) :: recv_count
  INTEGER, DIMENSION(0:mpi_size-1) :: displ
  INTEGER :: ierr
  INTEGER :: month,month1, month2,dec,year
  REAL :: frac1, frac2
  REAL, DIMENSION(klon) :: tmpeflux
  REAL, DIMENSION(ngridmx) :: eflux_glo_tmp
  INTEGER idiagd,idiagl,idiagp,idiage,idiagb
  INTEGER :: outdiag
  LOGICAL ::  outwfunc
  LOGICAL :: convoh
  REAL masseair(ngridmx,llm) ! masse d'air ds chaque maille en vecteur
  REAL massecourt(iip1-1,llm,jj_nb) ! masse d'air pr transfert en vect
  REAL massecourt_glo(iip1-1,llm,jjp1) ! masse d'air pr transfert en vect
  INTEGER, PARAMETER      :: nbprodmax=50 !nb max voies differentes de prod
  INTEGER, PARAMETER      :: nblossmax=50 ! et loss chimique
  INTEGER  :: nbprod(iqmax),nbloss(iqmax)
  REAL :: sommedep(2),sommesrc(iqmax),sommeburd(llm,iqmax)
  REAL :: sommeloss(llm,iqmax),sommeprod(llm,iqmax),dtps
  !REAL :: sommesrcmap(ngridmx),sommelossmap(llm-1,ngridmx)
  INTEGER :: nr,nl
  REAL :: sumd0(2,klon),sume0(iqmax,klon)
  REAL :: suml0(llm-1,iqmax,nblossmax,klon),sump0(llm-1,iqmax,nbprodmax,klon)
  REAL :: sumd0_glo(2,ngridmx),sume0_glo(iqmax,ngridmx)
  REAL :: suml0_glo(llm-1,iqmax,nblossmax,ngridmx)
  REAL :: sump0_glo(llm-1,iqmax,nbprodmax,ngridmx)
  REAL :: qfi_tmp(llm,iqmax,klon),qfi_glo(llm,iqmax,ngridmx)
  INTEGER, DIMENSION(100) :: nlevavgmax,nbformul
  character*7 :: str
  integer :: nbsat,idsat,nlevavg
  real :: nlevavglu
  REAL :: coefffc
  CHARACTER*132 :: ficinfo,ficdefsat
  REAL, ALLOCATABLE, DIMENSION(:)  :: pavg_mid,pavg, qa0, ak, qa0lu,dpavg
  INTEGER, DIMENSION(100) :: nbinfofile
  integer :: infos
  INTEGER :: iq1,nblinfo,nlev,nbs,chosenlev
  REAL, ALLOCATABLE, DIMENSION(:) :: qintavg, qintavg_tl, tmp_tl
  REAL, DIMENSION(llm)     :: veci, veci_tl
  REAL, ALLOCATABLE, DIMENSION(:) :: vecf, vecf_tl
  REAL, ALLOCATABLE, DIMENSION(:) :: amat,deltap
  REAL, ALLOCATABLE, DIMENSION(:) :: qprim
  CHARACTER*132 :: ficqprim
  REAL :: convpres,mair,conv
  CHARACTER*2 :: monthstr,idstr,nbsstr
  CHARACTER*4 :: yearstr
  !
  
  EXTERNAL :: exner_hyb
  !
  unskap   = 1./ kappa
  
  ALLOCATE( tdyn(iip1,jjp1,llm) )
  ALLOCATE( qfi(klon,llm,iqmax) )
  ALLOCATE( qfi_tl(klon,llm,iqmax) )
  ALLOCATE( qfi_ad(klon,llm,iqmax) )
  ALLOCATE( eflux(klon,iqmax) )
  ALLOCATE( eflux_tl(klon,iqmax) )
  ALLOCATE( source(iip1,jjp1,iqmax) )
  ALLOCATE( source_tl(iip1,jjp1,iqmax) )
  ALLOCATE( paprs(klon,llm+1) )
  ALLOCATE( pplay(klon,llm) )
  IF (.NOT. fwd) THEN
    ALLOCATE(qstotraj(iip1,jjp1,llm,iqmax,nsplit*nsplit_dyn))
    ALLOCATE(qfistotraj(klon,llm,iqmax,nsplit*nsplit_phy))
    ALLOCATE(efluxstotraj(klon,iqmax,nsplit*nsplit_phy))
    ALLOCATE(prescrstotraj(klon,llm,nsplit*nsplit_phy,iprescrmax))
    ALLOCATE(prodstotraj(klon,llm,nsplit*nsplit_phy,iprodmax))
  ENDIF
  
  IF (fwd) THEN
    eflux_tl(:,:) = 0.
    eflux(:,:) = 0.
  ENDIF
  
  ! open information on satellites
  ficdefsat='infousedsat.txt'
  open(70,file=ficdefsat,status='old',form='formatted',action='read')
  read(70,*)nbsat
  print*,'Number of satellites used:',nbsat
  WRITE(monthstr,'(i2.2)')month
  WRITE(yearstr,'(i4.4)')year
  IF (nbsat > 0) THEN
    do nbs=1,nbsat
      read(70,*)str,idsat
      read(70,*)nlevavgmax(idsat),nbformul(idsat)
      write(idstr,'(i2.2)')idsat
      ficinfo='infos_files/infos'//idstr//yearstr//monthstr//'.bin'
      infos = 70+nbs !OK pr le moment...
      PRINT*,'INFOS file',ficinfo
      OPEN(infos,file=ficinfo,form='unformatted',access='direct',recl=8*(1+nlevavgmax(idsat)*3+1 ))
      nbinfofile(idsat)=infos
    enddo
    1050 close(70)
  ENDIF
  !
  
  
  DO j=0,mpi_size-1
    nobs_para(j)=0
    DO i=1,nobs_glo
      IF (tabobs_glo(i,3)>=jj_begin_para(j) .AND. tabobs_glo(i,3)<=jj_END_para(j)) &
              nobs_para(j)=nobs_para(j)+1
    END DO
  END DO
  nobs=nobs_para(mpi_rank)
  
  ALLOCATE(tabobs(nobs,10))
  nobs=0
  DO i=1,nobs_glo
    IF (tabobs_glo(i,3)>=jj_begin .AND. tabobs_glo(i,3)<=jj_end) THEN
      nobs=nobs+1
      tabobs(nobs,:)=tabobs_glo(i,:)
    ENDIF
  ENDDO
  
  ! SAVE concentrations for adjoint and satellites
  WRITE(nbsstr,'(i2.2)')mpi_rank
  ficqprim='qprim'//nbsstr//yearstr//monthstr//'.bin'
  
  if (nbsat > 0) then
    nbs = 60+mpi_rank
    IF (fwd) THEN
      OPEN(nbs, FILE=ficqprim,FORM='unformatted', access='direct',recl=8*llm)  ! 64 bits
    ELSE
      OPEN(nbs, FILE=ficqprim,form='unformatted',access='direct',action='read',recl=8*llm)
    ENDIF
  end if

  !-----------------------------------------------------------------------
  !     Starting the time loop:
  !     ----------------------------
  !     MG initialisation of iday and time
  iday=day_ini+1
  !gch
  IF(.not.fwd)iday=day_ini+nday
  IF(time_0 >  1.) iday = iday+1
  !
  
  CALL InitTime
  CALL Init_timer
  CALL start_timer(timer_vanleer)
  CALL suspend_timer(timer_vanleer)
  
  DO itau=1,nday*iday_step
    
    IF(mpi_rank==0) PRINT*,'ITAU =',itau, nday, iday_step

    ! CHEMISTRY
    IF(do_chemistry)THEN
      IF(mod(itau-1,iday_step) == 0) THEN
        CALL readchem(iday-day_ini, &
                refprescr_glo,refprod_glo, refjrates_glo, &
                temp_glo, pmid_glo, convoh)
        
        CALL scatter_phy(refprod_glo, refprod)
        CALL scatter_phy(refprescr_glo, refprescr)
        CALL scatter_phy(refjrates_glo, refjrates)
        CALL scatter_phy(temp_glo,temp)
        CALL scatter_phy(pmid_glo,pmid)
        CALL scatter_phy(depvel_glo,depvel_loc)

      ENDIF
    END IF

    ! Split handling
    dtdyn=dtvr/float(nsplit*nsplit_dyn)
    IF (fwd) THEN
      irec=(day0-1)*iday_step+itau
    ELSE
      ! for the adjoint
      irec=day0*iday_step-itau+1
    ENDIF
    
    !-----------------------------------------------------------------------
    !     Reading mass fluxes and physics
    !     -----------------------------------------------------
    CALL VTb(VTlecfluxnc)
    CALL readfluxnc_p(irec,masse,pbaru,pbarv,w,teta, phi, &
            nrec,areafi,pphis, &
            t,zmfu, zmfd, zde_u,zen_d, coefkz, &
            dake, mpke, phike, updke, dndke, wghtke,&
            entr_therm,fm_therm,&
            yu1, yv1, ftsol, pctsrf, phis)
    CALL VTe(VTlecfluxnc)
    
    DO j=jj_begin,jj_end
      DO i=1,iip1
        smass(i,j)=0.
      ENDDO
    ENDDO
    DO l=1,llm
      DO j=jj_begin,jj_end
        DO i=1,iip1
          smass(i,j)=smass(i,j)+masse(i,j,l)
        ENDDO
      ENDDO
      DO j=1,jj_nb
        DO i=1,iip1-1
          massecourt(i,l,j)=masse(i,jj_begin+j-1,l)
        ENDDO
      ENDDO
    ENDDO
    
    displ(:)=0
    DO i=0,mpi_size-1
      recv_count(i)=jj_nb_para(i)*(iip1-1)*llm
      IF (i /= 0)displ(i)=displ(i-1)+recv_count(i-1)
    ENDDO
    
    CALL mpi_gatherv(massecourt,(iip1-1)*llm*jj_nb,mpi_real8,massecourt_glo,&
            recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
    IF(mpi_rank==0) then
      ! transfert massecourt en vecteur
      DO l=1,llm
        masseair(1,l) = massecourt_glo(1,l,1)
        masseair(ngridmx,l) = massecourt_glo(1,l,jjp1)
        DO j = 2, jjp1-1
          ig = 2+(j-2)*(iip1-1)
          DO i = 1, iip1-1
            masseair(ig-1+i,l) = massecourt_glo(i,l,j)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    
    DO j=jj_begin,jj_end
      DO i=1,iip1
        ps (i,j)=smass(i,j)/area(i,j)*g
      END DO
    END DO
    
    !     
    !print*,'     calcul de la pression et de pk (fonction d Exner)'
    !
    CALL pression( ip1jmp1, ap, bp, ps, p       )
    CALL exner_hyb( ip1jmp1, ps, p,alpha,beta, pks, pk, pkf )
    !
    !     Calcul de paprs et pplay 
    !       paprs  definis aux (llm +1) interfaces des couches
    !       pplay  definis aux (  llm )    milieux des couches
    !     
    DO l = 1, llmp1
      DO ig0=1,klon
        i=index_i(ig0)
        j=index_j(ig0)
        paprs(ig0,l)=p(i,j,l)
      ENDDO
    ENDDO
    !     
    DO l=1,llm
      DO ig0=1,klon
        i=index_i(ig0)
        j=index_j(ig0)
        pksurcp        = pk(i,j,l) / cpp
        pplay(ig0,l)   = preff * pksurcp ** unskap
      ENDDO
    ENDDO
    !
    DO jj=jj_begin,jj_end
      DO ii=1,iip1
        phis(ii,jj)= phi(ii,jj,1)-teta(ii,jj,1)*(pks(ii,jj)-pk(ii,jj,1))
      END DO
    END DO
    
    ! print*,'     Initialization of weighting function'
    IF (fwd .AND. itau == 1 .AND. outwfunc) THEN
      wfunc_int(:,:,:,:) = 0.
      
      DO jj=jj_begin,jj_end
        DO ii=1,iip1
          
          DO i = 1, llm
            p_mid(llm+1-i) = (p(ii,jj,i)+p(ii,jj,i+1))/2.
            dp(i) = p(ii,jj,i)-p(ii,jj,i+1)
          ENDDO
          
          isec = 1
          DO i = 1, llm
            wfunc_int(i,isec,ii,jj) = dp(i)
          ENDDO
          wfunc_int(:,isec,ii,jj) = wfunc_int(:,isec,ii,jj) / SUM( wfunc_int(:,isec,ii,jj) )
        ENDDO
      ENDDO
      
      CALL write_func_int( wfunc_int,nsec)
    ENDIF
    
    !=======================================================================
    ! STARTING A SUB TIMESTEP (inside an archiving step)
    !=======================================================================
    
    IF (fwd) THEN
      indexq = indexq + 1
      CALL gather_field(q,ip1jmp1,llm*iqmax,0)
      IF (mpi_rank==0) THEN
        DO i = 1, iqmax
          WRITE(itrajq(i),rec=indexq) q(:,:,:,i)
        END DO
      ENDIF
    ELSE
      DO i = 1, iqmax
        READ(itrajq(i),rec=indexq) q(:,:,:,i)
      END DO
      indexq = indexq - 1
    ENDIF
    
    !For tests on a single day
    !   IF (fwd .AND. itau > 2 * iday_step) GOTO 130
    !   IF ((.NOT. fwd) .AND. itau < (iday_step*(nday-1)+1)) GOTO 130
    !
    
    DO isplit=1,nsplit

!      WRITE(*,*) 'Beginning of isplit step'
      
      IF (fwd) THEN
        itausplit=(itau-1)*nsplit+isplit ! "smaller" time step
        
        ! Extract inputs for j
        ! Assumes daily inputs at the moment
        days_done = (itau-1) / iday_step
        j = MIN( days_done/periodflux+1, nbounds )
        i1 = j
        i2 = j+1
        t1 = bounddays(j)*nsplit*iday_step+1.
        t2 = bounddays(j+1)*nsplit*iday_step+1.
      ELSE
        ! backward for itau, forward for isplit...
        itausplit=(itau-1)*nsplit+nsplit-isplit+1
        
        ! Extract j
        days_done = nday - itau / iday_step
        j = MIN( days_done/periodflux+1, nbounds )
      ENDIF
      
      ! fluxes from control vector
      days_done = (itau-1) / iday_step
      j = MIN( days_done/periodflux+1, nbounds )
      i1 = j
      i2 = j+1
      t1 = bounddays(j)*nsplit*iday_step+1.
      t2 = bounddays(j+1)*nsplit*iday_step+1.
      
!        DO i = 1, iip1
      DO iq = 1, iqmax
        source_tl(:,jj_begin:jj_end,iq) = sflux_tl(:,jj_begin:jj_end,i1,iq) !&
!                    + (sflux_tl(i,jj_begin:jj_end,i2,iq)-sflux_tl(i,jj_begin:jj_end,i1,iq)) * (itausplit - t1) &
!                            /(t2 - t1)
        source(:,jj_begin:jj_end,iq) = sflux(:,jj_begin:jj_end,i1,iq) !&
!                    + (sflux(i,jj_begin:jj_end,i2,iq)-sflux(i,jj_begin:jj_end,i1,iq)) * (itausplit - t1) &
!                            /(t2 - t1)
      ENDDO
      !gch
      ! Scale factors
      IF(iprodmax > 0)THEN
        DO iq = 1, iprodmax
          intpscale_tl(:,jj_begin:jj_end, iq) = pscale_tl(:,jj_begin:jj_end,i1, iq) !&
!                      + (pscale_tl(i,jj_begin:jj_end,i2, iq)-pscale_tl(i,jj_begin:jj_end,i1, iq)) &
!                              * (itausplit - t1)/(t2 - t1)
          intpscale(:,jj_begin:jj_end, iq) = pscale(:,jj_begin:jj_end,i1, iq) !&
!                      + (pscale(i,jj_begin:jj_end,i2, iq)-pscale(i,jj_begin:jj_end,i1, iq)) &
!                              * (itausplit - t1) /(t2 - t1)
        END DO
      END IF
      
      IF(iprescrmax > 0)THEN
        DO iq = 1, iprescrmax
          intscale_tl(:,jj_begin:jj_end, iq) = scale_tl(:,jj_begin:jj_end,i1, iq) !&
!                      + (scale_tl(i,jj_begin:jj_end,i2, iq)-scale_tl(i,jj_begin:jj_end,i1, iq)) &
!                              * (itausplit - t1)/(t2 - t1)
          intscale(:,jj_begin:jj_end, iq) = scale(:,jj_begin:jj_end,i1, iq) !&
!                      + (scale(i,jj_begin:jj_end,i2, iq)-scale(i,jj_begin:jj_end,i1, iq)) &
!                              * (itausplit - t1) /(t2 - t1)
        END DO
      END IF
!        ENDDO
      
      CALL gr_dyn_fi_p(iqmax,iip1,jjp1,klon,source_tl,eflux_tl)
      CALL gr_dyn_fi_p(iqmax,iip1,jjp1,klon,source,eflux)
      
      ! SACS
      if (iprodmax > 0) then
        CALL gr_dyn_fi_p(iprodmax,iip1,jjp1,klon,intpscale_tl,phypscale_tl)
        CALL gr_dyn_fi_p(iprodmax,iip1,jjp1,klon,intpscale,phypscale)
      endif
      if (iprescrmax > 0) then
        CALL gr_dyn_fi_p(iprescrmax,iip1,jjp1,klon,intscale_tl,physcale_tl)
        CALL gr_dyn_fi_p(iprescrmax,iip1,jjp1,klon,intscale,physcale)
      endif

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
      dec=0 ! Shift global sflx_ch2o
      IF (mpi_rank>0)dec=klon_mpi_para_end(mpi_rank-1)
      DO i = 1, klon
        if (iprodmax > 0) then
          DO iq = 1,iprodmax
            refprod_tl(i,:,iq) = phypscale_tl(i,iq) * refprod(i,:,iq)
            refprod(i,:,iq)    = phypscale(i,iq)    * refprod(i,:,iq)
          END DO
        end if
        
        if (iprescrmax > 0) then
          DO iq = 1,iprescrmax
            refprescr_tl(i,:,iq) = physcale_tl(i,iq) * refprescr(i,:,iq)
            refprescr(i,:,iq)    = physcale(i,iq)    * refprescr(i,:,iq)
          END DO
        end if
      ENDDO
      
!
!        !
!        CALL gather_field(source,ip1jmp1,iqmax,0)
!
!        ! SACS
!        displ(:)=0
!        DO i=0,mpi_size-1
!          recv_count(i)=klon_mpi_para_nb(i)*iqmax
!          if(i /= 0)displ(i)=displ(i-1)+recv_count(i-1)
!        ENDDO
!
!        ! When gathering, it is the last dimension that must be gathered!!
!        ! Fluxes
!        DO iq=1,iqmax
!          tmpeflux(:)=eflux(:,iq)
!          eflux_glo_tmp(:)=0.
!          displ(:)=0
!          DO i=0,mpi_size-1
!            recv_count(i)=klon_mpi_para_nb(i)
!            IF (i /= 0)displ(i)=displ(i-1)+recv_count(i-1)
!          ENDDO
!          call mpi_gatherv(tmpeflux,klon,mpi_real8,eflux_glo_tmp, &
!                  recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
!          eflux_glo(:,iq)=eflux_glo_tmp(:)
!        ENDDO !iqmax
!
!        ! Production
!        if (iprodmax > 0) then
!          DO iq=1,iprodmax
!            tmpeflux(:)=phypscale(:,iq)
!            eflux_glo_tmp(:)=0.
!            displ(:)=0
!            DO i=0,mpi_size-1
!              recv_count(i)=klon_mpi_para_nb(i)
!              IF (i /= 0)displ(i)=displ(i-1)+recv_count(i-1)
!            ENDDO
!            call mpi_gatherv(tmpeflux,klon,mpi_real8,eflux_glo_tmp, &
!                    recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
!            phypscale_glo(:,iq)=eflux_glo_tmp(:)
!          END DO
!        end if
!
!        ! Prescribed fields
!        if (iprescrmax > 0) then
!          DO iq=1,iprescrmax
!            tmpeflux(:)=physcale(:,iq)
!            eflux_glo_tmp(:)=0.
!            displ(:)=0
!            DO i=0,mpi_size-1
!              recv_count(i)=klon_mpi_para_nb(i)
!              IF (i /= 0)displ(i)=displ(i-1)+recv_count(i-1)
!            ENDDO
!            call mpi_gatherv(tmpeflux,klon,mpi_real8,eflux_glo_tmp, &
!                    recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
!            physcale_glo(:,iq)=eflux_glo_tmp(:)
!          END DO
!        end if
!
!        ! Saving temporal values for later computation by the adjoint
!        IF (mpi_rank==0) THEN
!          DO iq = 1, iqmax
!            WRITE(itrajflx(iq),rec=itausplit)eflux_glo(:,iq)
!          END DO
!
!          IF(iprescrmax > 0) THEN
!            DO iq = 1, iprescrmax
!              WRITE(itrajscale(iq),rec=itausplit)physcale_glo(:,iq)
!            END DO
!          ENDIF
!
!          IF(iprodmax > 0) THEN
!            DO iq = 1, iprodmax
!              WRITE(itrajpscale(iq),rec=itausplit)phypscale_glo(:,iq)
!            END DO
!          ENDIF
!        ENDIF
!
!      ELSE
!        CALL Vtb(VTread_source)
!
!        ! Reading temporal inputs as computed by forward mode
!        IF (mpi_rank==0) THEN
!          DO iq = 1, iqmax
!            READ(itrajflx(iq),rec=(nday * nsplit * iday_step + 1 - itausplit)) &
!                    eflux_glo(:,iq)
!          END DO
!
!          IF(iprescrmax > 0) THEN
!            DO iq = 1, iprescrmax
!              READ(itrajscale(iq),rec=(nday * nsplit * iday_step + 1 - itausplit)) &
!                      physcale_glo(:,iq)
!            END DO
!          ENDIF
!
!          IF(iprodmax > 0) THEN
!            DO iq = 1, iprodmax
!              READ(itrajpscale(iq),rec=(nday * nsplit * iday_step + 1 - itausplit)) &
!                      phypscale_glo(:,iq)
!            END DO
!          ENDIF
!        ENDIF
!
!        CALL scatter_phy(eflux_glo,eflux)
!        IF(iprodmax > 0)CALL scatter_phy(phypscale_glo,phypscale)
!        IF(iprescrmax > 0)CALL scatter_phy(physcale_glo,physcale)
!
!        DO i = 1,klon
!          IF(iprodmax > 0) THEN
!            DO iq = 1,iprodmax
!              refprod(i,:,iq) = phypscale(i,iq) * refprod(i,:,iq)
!            END DO
!          END IF
!          IF(iprescrmax > 0) THEN
!            DO iq = 1,iprescrmax
!              refprescr(i,:,iq) = physcale(i,iq) * refprescr(i,:,iq)
!            END DO
!          END IF
!        ENDDO
!      ENDIF
      !
      !     
      ! print*,'     FWD statement of advection'
      !     
      
      DO iii=1,nsplit_dyn
        IF (.NOT. fwd) &
                qstotraj(:,jj_begin:jj_end,:,:,(isplit-1)*nsplit_dyn+iii) = q(:,jj_begin:jj_end,:,:)
        DO iq=1,iqmax
          IF (tangent) THEN
            CALL VTb(VTadvection)
            CALL vlsplt_tl_p(q(1,1,1,iq),pente_max,masse,w,pbaru,pbarv,dtdyn,q_tl(1,1,1,iq))
            CALL VTe(VTadvection)
          ELSE
            CALL VTb(VTadvection)
            CALL vlsplt_p(q(1,1,1,iq),pente_max,masse,w,pbaru,pbarv,dtdyn)
            CALL VTe(VTadvection)
          ENDIF
        ENDDO
      ENDDO


      !-----------------------------------------------------------------------
      ! print*,'     PHYSIQUE'
      !-----------------------------------------------------------------------
      IF(physic) THEN
      
        !     Translating tracer to physical map
        IF (tangent) CALL gr_dyn_fi_p(llm*iqmax,iip1,jjp1,klon,q_tl,qfi_tl)
        CALL gr_dyn_fi_p(llm*iqmax,iip1,jjp1,klon,q,qfi)
        
        !-----------------------------------------------------------------------
        !     Physical Loop
        !-----------------------------------------------------------------------
        DO iii=1,nsplit_phy
          !
          !     FWD statement of physics
          ! print*,'     (used in AD)'
          !

!          WRITE(*,*) 'Beginning of physical step', iii, nsplit_phy
          IF (.NOT. fwd) THEN
            efluxstotraj(1:klon,:,(isplit-1)*nsplit_phy+iii) = eflux(1:klon,:)
            qfistotraj(1:klon,:,:,(isplit-1)*nsplit_phy+iii) = qfi(1:klon,:,:)
            !gch
            IF(iprescrmax > 0) prescrstotraj(:,:,(isplit-1)*nsplit_phy+iii,:) = refprescr(:,:,:)
            IF(iprodmax > 0) prodstotraj(:,:,(isplit-1)*nsplit_phy+iii,:) = refprod(:,:,:)
            !
          ENDIF
          
          CALL VTb(VTphysiq)
          IF (tangent) THEN
            CALL phytrac_tl(dtphys, t, paprs, pplay, &
                    zmfu, zmfd, zde_u, zen_d, coefkz, qfi, eflux, eflux_tl, qfi_tl, &
                    dake,phike,mpke,updke,&
                    dndke,wghtke,entr_therm,fm_thermp1, &
                    refprod,refprescr,temp, pmid,refjrates,depvel_loc, &
                    refprod_tl,refprescr_tl)
          ELSE
            
!            IF (isplit==1) THEN !! Works only with nsplit_phy=1
!              do iq=1,iqmax
!                sommesrc(iq)=0.
!                do l=1,llm
!                  sommeprod(l,iq)=0.
!                  sommeloss(l,iq)=0.
!                  sommeburd(l,iq)=0.
!                enddo
!              enddo
!
!              dtps=0.
!            ENDIF
            
            CALL phytrac(dtphys, t, paprs, pplay, areafi, zmfu, zmfd, zde_u, zen_d, coefkz, qfi, eflux, &
                    entr_therm,fm_thermp1,&
                    dake,phike,mpke,updke,dndke,wghtke, &
                    refprod,refprescr,temp, pmid,refjrates,depvel_loc, &
                    sumd0,sump0,suml0,sume0,nbprodmax,nblossmax)

            dtps=dtps+dtphys
            displ(:)=0
            DO i=0,mpi_size-1
              recv_count(i)=klon_mpi_para_nb(i)*2
              IF (i /= 0) displ(i)=displ(i-1)+recv_count(i-1)
            ENDDO
            call mpi_gatherv(sumd0,klon*2,mpi_real8,sumd0_glo, &
                    recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
            displ(:)=0
            DO i=0,mpi_size-1
              recv_count(i)=klon_mpi_para_nb(i)*iqmax
              IF (i /= 0) displ(i)=displ(i-1)+recv_count(i-1)
            ENDDO
            call mpi_gatherv(sume0,klon*iqmax,mpi_real8,sume0_glo, &
                    recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
            displ(:)=0
            DO i=0,mpi_size-1
              recv_count(i)=klon_mpi_para_nb(i)*iqmax*nbprodmax*(llm-1)
              IF (i /= 0) displ(i)=displ(i-1)+recv_count(i-1)
            ENDDO
            call mpi_gatherv(sump0,klon*iqmax*nbprodmax*(llm-1),mpi_real8,sump0_glo, &
                    recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
            displ(:)=0
            DO i=0,mpi_size-1
              recv_count(i)=klon_mpi_para_nb(i)*iqmax*nblossmax*(llm-1)
              IF (i /= 0) displ(i)=displ(i-1)+recv_count(i-1)
            ENDDO
            call mpi_gatherv(suml0,klon*iqmax*nblossmax*(llm-1),mpi_real8,suml0_glo, &
                    recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
            displ(:)=0
            DO i=0,mpi_size-1
              recv_count(i)=klon_mpi_para_nb(i)*iqmax*llm
              IF (i /= 0) displ(i)=displ(i-1)+recv_count(i-1)
            ENDDO
            DO l=1,llm
              DO iq=1,iqmax
                DO nl=1,klon
                  qfi_tmp(l,iq,nl)=qfi(nl,l,iq)
                ENDDO
              ENDDO
            ENDDO
            call mpi_gatherv(qfi_tmp,klon*iqmax*llm,mpi_real8,qfi_glo, &
                    recv_count,displ,mpi_real8,0,comm_lmdz,ierr)
            IF (mpi_rank==0) then
              IF (outdiag==1) THEN
                do iq=1,iqmax
                  DO nl=1,ngridmx
                    sommesrc(iq)=sommesrc(iq)+sume0_glo(iq,nl)*1.e-9
                    do l=1,llm-1
                      sommeburd(l,iq)=sommeburd(l,iq)+qfi_glo(l,iq,nl)*masseair(nl,l)*1.e-9
                      do nr=1,nbprod(iq)
                        sommeprod(l,iq)=sommeprod(l,iq)+sump0_glo(l,iq,nr,nl)*masseair(nl,l)*1.e-9
                      enddo
                      do nr=1,nbloss(iq)
                        sommeloss(l,iq)=sommeloss(l,iq)+suml0_glo(l,iq,nr,nl)*masseair(nl,l)*1.e-9
                      enddo
                    
                    enddo !l
                  
                  enddo !nl
                enddo !iq

                ! TODO : deal with it
!                if(isplit==nsplit) then
!                  !print*, 'DTPS',dtps
!                  write(idiagl,*)itau,(sommeloss(l,id_CH4),l=1,llm)
!                  write(idiagp,*)itau,(sommeprod(l,id_CH4),l=1,llm)
!                  write(idiage,*)itau,sommesrc(id_CH4)
!                  write(idiagb,*)itau,(sommeburd(l,id_CH4),l=1,llm)
!                endif !isplit
              ENDIF ! outdiag
            ENDIF !mpi_rank
          
          ENDIF !
          CALL VTe(VTphysiq)
!          WRITE(*,*) 'END of physical step', iii, nsplit_phy

          !
        ENDDO             ! nsplit_phy
        !
        ! print*,'     On passe le traceur sur la grille dynamique.'
        IF (tangent) CALL gr_fi_dyn_p(llm*iqmax,klon,iip1,jjp1,qfi_tl,q_tl)
        CALL gr_fi_dyn_p(llm*iqmax,klon,iip1,jjp1,qfi,q)
      
      ENDIF               ! if (physic)


      !PRINT*, 'av DEBUT', tabobs(i,5)
      
      !     FWD statement of comparison with obs
      IF (fwd) THEN
      IF (ANY( tabobs(:,1) + tabobs(:,2) > itausplit ) &
               .AND. ANY( tabobs(:,1) <= itausplit )) THEN
        DO i = 1, nobs
          jj=tabobs(i,3)
          ii=tabobs(i,4)
          IF ( tabobs(i,1) <= itausplit .AND. tabobs(i,1) + tabobs(i,2) > itausplit) THEN
            iq = ABS (INT( tabobs(i,6) ))
            IF (iq <= iqmax .and. iq > 0) THEN
              IF (tabobs(i,6) < 0 ) THEN !  satellite
                ! iq = numero espece = la centaine de tabobs(i,5) pr les sat
                ! iq1 = numero du sat = le nombre d'unite de tabobs(i,5)
                iq = ABS (INT( tabobs(i,6)/100. ))
                iq1 = INT ( ABS(tabobs(i,6))-iq*100)
                ! lecture des pavg, ak et qa0 voulus
                ! decodage du numero de ligne voulu ds fic du mois
                nblinfo=int((int(tabobs(i,6))-tabobs(i,6))*10000000.) ! XXXX ATTENTION codage en dur,
                ! doit etre coherent avec codage dans python XXX
                IF( nbformul(iq1)==1 .OR. nbformul(iq1)==2 ) THEN !  general formula (Warning diff 1-2 in mair and equation used)
                  READ(nbinfofile(iq1),rec=nblinfo) nlevavglu
                  nlevavg=int(nlevavglu)
                  ALLOCATE(pavg(nlevavg+1)) ! pressure levels incl. surface pressure
                  ALLOCATE(pavg_mid(nlevavg)) ! pressure at the middle of the layer
                  ALLOCATE(qa0lu(nlevavg)) ! prior profile
                  ALLOCATE(ak(nlevavg))  ! averaging kernels
                  read(nbinfofile(iq1),rec=nblinfo)nlevavglu,ak(:),qa0lu(:),pavg(:)
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
                  ALLOCATE(qprim(llm))
                  DO l = 1, llm
                    veci_tl(llm+1-l) = q_tl(tabobs(i,4),tabobs(i,3),l,iq)
                    veci(llm+1-l) = q(tabobs(i,4),tabobs(i,3),l,iq)
                    qprim(l)= q(tabobs(i,4),tabobs(i,3),l,iq)
                  ENDDO
                  nbs=mpi_rank+60
                  WRITE(nbs,rec=i)qprim
                  ALLOCATE(vecf(nlevavg))
                  ALLOCATE(vecf_tl(nlevavg))
                  CALL intex_tl (llm,nlevavg,p_mid,pavg_mid,veci,vecf,veci_tl,vecf_tl)
                  ! correction by FC
                  coefffc = dot_product(veci,dpp) / dot_product(vecf,dpavg) &
                          / sum(dpp) * sum(dpavg)
                  vecf_tl(:) = vecf_tl(:) * coefffc
                  vecf(:) = vecf(:) * coefffc
                  ! PRINT*, 'coeff', coefffc,' obs numero',i ! veci, vecf
                  ! surface index 1 as in LMDz
                  ALLOCATE(qintavg(nlevavg))
                  ALLOCATE(qintavg_tl(nlevavg))
                  DO l = 1, nlevavg
                    qintavg_tl(nlevavg+1-l) = vecf_tl(l)
                    qintavg(nlevavg+1-l) = vecf(l)
                  ENDDO
                  !model-equivalent of observation, in kg/kg
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
                  mr=0.
                  mr_tl=0.
                  DO nlev=1,nlevavg
                    IF (nbformul(iq1)==1) THEN
                      mair=mair+deltap(nlev)
                      mr_tl = mr_tl + qintavg_tl(nlev)*amat(nlev) * deltap(nlev)
                      mr = mr+ (qa0(nlev) + (qintavg(nlev)-qa0(nlev))*amat(nlev))*deltap(nlev)
                    ENDIF
                    IF (nbformul(iq1)==2) THEN
                      mair=mair+amat(nlev)*deltap(nlev)
                      mr_tl = mr_tl + qintavg_tl(nlev)*amat(nlev) * deltap(nlev)
                      mr = mr+ qintavg(nlev)*amat(nlev)*deltap(nlev)
                    ENDIF
                  ENDDO
                  mr_tl=mr_tl/mair
                  mr=mr/mair

                  DEALLOCATE(vecf)
                  DEALLOCATE(vecf_tl)
                  DEALLOCATE(qintavg)
                  DEALLOCATE(qintavg_tl)
                  DEALLOCATE(amat)
                  DEALLOCATE(deltap)
                  DEALLOCATE(qa0lu)
                  !               DEALLOCATE(tmp_tl)
                  DEALLOCATE(qprim)
                  DEALLOCATE(dpavg)
                  deallocate(pavg)
                  deallocate(pavg_mid)
                  deallocate(qa0)
                  deallocate(ak)
                ENDIF ! general formula
              ELSE      ! surface
                l = nint( (tabobs(i,6) - iq)*100 ) ! extract level information
                mr_tl = q_tl(tabobs(i,4),tabobs(i,3),l,iq)
                mr = q(tabobs(i,4),tabobs(i,3),l,iq)
              ENDIF ! tabobs < 0
              ppmv_tl = mr_tl / adv_mass(iq) * dry_mass *1.e6
              ppmv = mr / adv_mass(iq) * dry_mass *1.e6
              tabobs(i,8) = tabobs(i,8) + ppmv_tl
              tabobs(i,7) = tabobs(i,7) + ppmv
            
              ! Save pressure
              l = nint( (tabobs(i,6) - iq)*100 )
              tabobs(i,9) = tabobs(i,9) + (p(ii,jj,l)+p(ii,jj,l+1))/2.
              tabobs(i,10) = tabobs(i,10) + p(ii,jj,l)-p(ii,jj,l+1)
            END IF
          ENDIF ! tabobs(i,1) == itausplit
        ENDDO ! i = 1, nobs
      ENDIF ! ( fwd .AND. ANY( tabobs(:,1) == itausplit ) )
      END IF
      !

    ENDDO                   ! isplit=1,nsplit


    IF (.NOT. fwd) THEN
      CALL timeloop_ad(iday_step, itau, nsplit, nsplit_dyn, &
              nsplit_phy, dtdyn, p, paprs, pplay, q, q_ad, qfi, &
              eflux, qstotraj, qfistotraj, efluxstotraj,&
              w, coefkz, nsec, wfunc_int, &
              periodflux, nbounds, bounddays, sflux_ad, masse, &
              pbaru, pbarv, t, zmfu, zmfd, zde_u, zen_d, pente_max, nobs, tabobs, &
              dake,mpke,phike,updke,dndke,wghtke,&
              entr_therm,fm_thermp1, &
              nbsat, nlevavgmax, nbinfofile,nbformul, ps, &
              prescrstotraj, prodstotraj, &
              refprescr, refprod, temp, pmid, refjrates, depvel_loc, &
              nblossmax, nbprodmax, iday, month, &
              pscale_ad, scale_ad)
    
    ENDIF
    
    !-----------------------------------------------------------------------
    !     Computes ucov and vcov for outputs
    !-----------------------------------------------------------------------
    
    CALL Register_Hallo(masse,ip1jmp1,llm,0,1,1,0,mpi_request)
    CALL SendRequest(mpi_request)
    CALL WaitRequest(mpi_request)
    
    CALL massbar(masse, massebx, masseby )
    
    DO l=1,llm
      DO j=jj_begin,jj_end
        DO i=1,iip1
          ucov(i,j,l)=pbaru(i,j,l)/massebx(i,j,l)*cu(i,j)*cu(i,j) &
                  /istdyn
        ENDDO
      ENDDO
    ENDDO
    
    jjb=jj_begin
    jje=jj_end
    IF (south_pole) jje=jj_end-1
    
    DO l=1,llm
      DO j=jjb,jje
        DO i=1,iip1
          vcov(i,j,l)=pbarv(i,j,l)/masseby(i,j,l)*cv(i,j)*cv(i,j) &
                  /istdyn
        ENDDO
      ENDDO
    ENDDO
    
    !gch
    130 CONTINUE
    IF(mod(itau,iday_step) == 0) THEN
      IF (fwd) THEN
        iday=iday+1
      ELSE
        iday=iday-1
      ENDIF
    ENDIF
    !
    !
    !=======================================================================
  ENDDO                     ! itau
  
  PRINT *, 'Total time spent on parallelization :',DiffTime()
  PRINT *, 'CPU time spent on parallelization :',DiffCpuTime()
  PRINT *, 'CPU time spent on vlsplt ',timer_running(timer_vanleer)
  
  DEALLOCATE( tdyn )
  DEALLOCATE( qfi )
  DEALLOCATE( qfi_tl )
  DEALLOCATE( qfi_ad )
  DEALLOCATE( eflux )
  DEALLOCATE( eflux_tl )
  DEALLOCATE( source )
  DEALLOCATE( source_tl )
  DEALLOCATE( paprs )
  DEALLOCATE( pplay )
  IF (.NOT. fwd) THEN
    DEALLOCATE(qstotraj)
    DEALLOCATE(qfistotraj)
    DEALLOCATE(efluxstotraj)
    DEALLOCATE(prescrstotraj)
    DEALLOCATE(prodstotraj)
  ENDIF
  
  IF (fwd) THEN
    !cla Special treatment for tabobs
    CALL gather_tabobs(tabobs,nobs,nobs_para,nobs_glo,tabobs_glo)
  
  ELSE
    !gch
    IF(iqmax > 0) CALL gather_field(q_ad,ip1jmp1,llm*iqmax,0)
    IF(iqmax > 0) CALL gather_field(sflux_ad,ip1jmp1,iqmax*nbounds,0)
    IF(iprescrmax > 0) CALL gather_field(scale_ad,ip1jmp1,iprescrmax * nbounds,0)
    IF(iprodmax > 0) CALL gather_field(pscale_ad,ip1jmp1,iprodmax * nbounds,0)
    !
  ENDIF
  
  RETURN
END SUBROUTINE timeloop
