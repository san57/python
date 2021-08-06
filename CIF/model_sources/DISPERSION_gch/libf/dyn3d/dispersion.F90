PROGRAM dispersion
  
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
  
  !
  !   LMDZT atmospheric tracer transport model, forward, tangent-linear and adjoint 
  !
  !   Original version           : F Hourdin, A Idelkadi (LMD)
  !   Tangent-linear and adjoint : F Chevallier (LSCE)
  !   Parallelization            : Y Meurdesoif (LSCE)
  !   Fortran90                  : C Carouge    (LSCE)
  !   Clean-up                   : F Chevallier (LSCE)
  !
  
  USE mod_const_mpi, ONLY: init_const_mpi
  USE parallel
  USE dimphy
  USE mod_interface_dyn_phys
  USE comgeomphy
  USE mod_hallo
  USE Bands
  ! ##MS## USE WRITE_field
  ! chez IP WRITE_field_phy
  USE WRITE_field_p
  USE SPECIES_NAME
  
  IMPLICIT NONE
  
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "comconst.h"
  INCLUDE "comdissnew.h"
  INCLUDE "comvert.h"
  INCLUDE "comgeom2.h"
  INCLUDE "logic.h"
  INCLUDE "temps.h"
  INCLUDE "control.h"
  INCLUDE "ener.h"
  INCLUDE "description.h"
  INCLUDE "serre.h"
  INCLUDE "tracstoke.h"
  INCLUDE "ajout.h"
  INCLUDE "clesph0.h"
  
  INTEGER :: iday ! Julian day
  INTEGER :: iday_step,day0
  REAL    :: time_0 ! Time of the day (% of day) at the beginning of the simulations
  REAL,dimension(20) :: clesphy0
  
  !   Dynamical variables
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: vcov, ucov
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)   :: q, q_tl, q_ad
  REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: phi, w
  REAL                      :: ps(iip1,jjp1)
  REAL,DIMENSION(iip1,jjp1) :: phis            ! Geopotential at ground level
  !
  !   Physical variables
  REAL, ALLOCATABLE, DIMENSION(:,:) :: zmfd, zen_d, zmfu, zde_u
  REAL, ALLOCATABLE, DIMENSION(:,:) :: t, coefkz, ftsol, pctsrf
  REAL, ALLOCATABLE, DIMENSION(:)   :: pphis, areafi, yu1, yv1
  REAL, ALLOCATABLE, DIMENSION(:,:) :: dake, mpke, updke,dndke, wghtke
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: phike
  REAL, ALLOCATABLE, DIMENSION(:,:) :: entr_therm,fm_therm

  ! ##MS## zde_d(:,:), zen_u(:,:), fracimpa, frac_nucl, ?? readfluxnc?
  !
  ! Intermediate dynamical variables for transport
  REAL :: pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm) ! mass fluxes
  REAL :: masse(iip1,jjp1,llm)
  REAL :: teta(iip1,jjp1,llm)
  INTEGER :: irec,nrec
  
  INTEGER :: nsplit_dyn,nsplit_phy,nsplit
  INTEGER :: iq,j,i
  LOGICAL :: readstart, fwd, tangent
  EXTERNAL :: iniconst, inifilr
  EXTERNAL :: inigeom
  
  ! Link with PYVAR
  REAL, ALLOCATABLE, DIMENSION(:,:):: tabobs
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: qscale, qscale_tl, qscale_ad
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: sflux, sflux_tl, sflux_ad, sflx
  INTEGER :: ndayloc, periodflux, indexq, nbounds
  INTEGER, ALLOCATABLE, DIMENSION(:) :: itrajq
  INTEGER :: outdiag
  LOGICAL :: outwfunc, footprint
  LOGICAL :: convoh
  LOGICAL :: file_exists
  
  ! For satellites
  INTEGER                                :: nsec, nobs
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:)  :: wfunc_int
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: bounddays
  
  ! For budget calculation
  INTEGER :: idiagp,idiagl,idiagd,idiage,idiagb
  CHARACTER*2 monthstr
  CHARACTER*4 yearstr
  !gch
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: scale, scale_tl, scale_ad
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: pscale, pscale_tl, pscale_ad
  REAL, DIMENSION(ngridmx,llm) :: depvel
  REAL, ALLOCATABLE, DIMENSION(:,:) :: depvel_glo
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: refprescr_glo, refprod_glo, refjrates_glo
  INTEGER :: month,year
  REAL, DIMENSION(ngridmx,llm) :: temp_glo, pmid_glo
  REAL :: sflux_max
  ! ##MS##
  ! itrajfp ?
  ! wfunc_intc, wfunc, wfunc_tmp, dimwf
  ! INTEGER :: ii, jj
  !
  CALL init_const_mpi
  CALL init_parallel
  CALL Read_Distrib
  CALL Init_Phys_lmdz(iim,jjp1,llm,mpi_size,distrib_phys)
  CALL set_bands
  CALL Init_interface_dyn_phys
  CALL SetDistrib(jj_Nb_Caldyn)
  CALL Init_Mod_hallo
  CALL InitComgeomphy
  
  !-----------------------------------------------------------------------
  !   Initialisations:
  !   ----------------
  descript = 'Run GCM LMDZ'
  
  ! Various parameters in totinput
  OPEN(1,file='totinput',form='formatted',action='read')
  PRINT*,'Effective number of tracers'
  READ(1,*) iqmax
  WRITE(*,*) iqmax
  PRINT*,'Activating SACS chemistry (T) or not (F)'
  READ(1,'(l1)') do_chemistry
  WRITE(*,'(l1)') do_chemistry
  PRINT*,'Splitting for dynamical time step'
  READ(1,*) nsplit_dyn
  WRITE(*,*) nsplit_dyn
  PRINT*,'Splitting for physical time step'
  READ(1,*) nsplit_phy
  WRITE(*,*) nsplit_phy
  PRINT*,'Reading start file (T)?'
  READ(1,'(l1)') readstart
  WRITE(*,'(l1)') readstart
  PRINT*,'Forward=T or backward=F?'
  READ(1,'(l1)') fwd
  WRITE(*,'(l1)') fwd
  PRINT*, 'Write output diagnostics?'
  READ(1,*) outdiag
  WRITE(*,*) outdiag
  PRINT*, 'Write output wfunc?'
  READ(1,*) outwfunc
  WRITE(*,*) outwfunc
  WRITE(*,*) 'Running footprint?'
  READ(1,*) footprint
  WRITE(*,*) footprint
  WRITE(*,*) 'which NDAYLOC, i.e. number of days in the months?'
  READ(1,*) ndayloc
  WRITE(*,*) ndayloc
  WRITE(*,*) 'Need OH conversion?'
  READ(1,*) convoh
  WRITE(*,*) convoh
  WRITE(*,*) 'Period length (in days) for fluxes'
  READ(1,*) periodflux
  WRITE(*,*) periodflux
  WRITE(*,*) 'Number of days in the simulations'
  READ(1,*) nday
  WRITE(*,*) nday
  WRITE(*,*) 'Month of the year'
  READ(1,*) month
  WRITE(*,*) month
  WRITE(*,*) 'Year'
  READ(1,*) year
  WRITE(*,*) year
  WRITE(*,*) 'Flag Convection'
  READ(1,*) iflag_con
  WRITE(*,*) iflag_con
  WRITE(*,*) 'Run the physics'
  READ(1,*) physic
  WRITE(*,*) physic
  CLOSE(1)
  
  ! Initializes transported species
  CALL init_chem()
  
  !
  ! Allocations
  !
  ALLOCATE( q(iip1,jjp1,llm,iqmax) )
  ALLOCATE( q_tl(iip1,jjp1,llm,iqmax) )
  ALLOCATE( q_ad(iip1,jjp1,llm,iqmax) )
  ALLOCATE( vcov(iip1,jjm,llm) )
  ALLOCATE( ucov(iip1,jjp1,llm) )
  ALLOCATE( phi(iip1,jjp1,llm) )
  ALLOCATE( w(iip1,jjp1,llm) )
  ALLOCATE( coefkz(klon,llm) )
  ALLOCATE( qscale(iip1,jjp1,iqmax) )
  ALLOCATE( qscale_ad(iip1,jjp1,iqmax) )
  ALLOCATE( qscale_tl(iip1,jjp1,iqmax) )
  ALLOCATE( zmfd(klon,llm) )
  ALLOCATE( zen_d(klon,llm) )
  ALLOCATE( zmfu(klon,llm) )
  ALLOCATE( zde_u(klon,llm) )
  ALLOCATE( dake(klon,llm) )
  ALLOCATE( mpke(klon,llm) )
  ALLOCATE( phike(klon,llm,llm) )
  ALLOCATE( updke(klon,llm) )
  ALLOCATE( dndke(klon,llm) )
  ALLOCATE( wghtke(klon,llm) )
  ALLOCATE( fm_therm(klon,llm) )
  ALLOCATE( entr_therm(klon,llm) )
  ALLOCATE( t(klon,llm) )
  ALLOCATE( pphis (klon) )
  ALLOCATE( areafi(klon) )
  ALLOCATE( ftsol(klon,nbsrf))
  ALLOCATE( pctsrf(klon,nbsrf))
  ALLOCATE( yu1(klon) )
  ALLOCATE( yv1(klon) )

  ! ##MS## nsec =1 toujours ?
  nsec  = 1
  ALLOCATE( wfunc_int(llm,nsec,iip1,jjp1) )

  !
  ! Read observations (FWD or AD statement)
  !
  INQUIRE(FILE='obs.bin', EXIST=file_exists)
  IF(file_exists)THEN
    OPEN(1,file='obs.txt',form='formatted',action='read')
    READ(1,*) nobs
    CLOSE(1)
    
    OPEN(UNIT=1, FILE="obs.bin", form='unformatted',&
            ACCESS="direct", action='read', recl=8 * 6 * nobs)
    
    ALLOCATE ( tabobs(nobs, 10) )
    tabobs(:,:) = 0.
    READ(1,rec=1) tabobs(:, 1:6)
  ELSE
    nobs=0
    ALLOCATE ( tabobs(nobs,10) )
    tabobs(:,:) = 0.
  END IF
  
  ! Number of periods in the state vector constraining the simulations
  ALLOCATE ( bounddays(nday+1) )
  nbounds = 0
  !ndayloc = 28
  !IF (footprint) THEN
  !  ndayloc = 30
  !ENDIF
  DO iday = 1, ndayloc, periodflux
    nbounds = nbounds + 1
    bounddays(nbounds) = iday-1
  ENDDO
  nbounds = nbounds + 1
  bounddays(nbounds) = nday

  !
  ! Archive of trajectory (for write in FWD mode, for read in AD mode)
  !
  ALLOCATE(itrajq(iqmax))
  IF (fwd) THEN
    IF (mpi_rank==0) THEN
      DO iq = 1, iqmax
        itrajq(iq) = 900 + iq
        OPEN(itrajq(iq),file='trajq_'//trim(species(iq)%name)//'.bin',form='unformatted',access='direct', &
                recl=8*iip1*jjp1*llm)  ! 64 bits
      END DO
    ENDIF
  ELSE
    DO iq = 1, iqmax
      itrajq(iq) = 900 + iq
      OPEN(itrajq(iq),file='trajq_'//trim(species(iq)%name)//'.bin',form='unformatted',access='direct', &
              recl=8*iip1*jjp1*llm,action='read')  ! 64 bits
    END DO
  ENDIF

  ! DIAGNOSTICS on the fly
  idiagp  = 93
  idiagl  = 94
  idiagd  = 95
  idiage  = 96
  idiagb  = 85
  WRITE(monthstr,'(i2.2)')month
  WRITE(yearstr,'(i4.4)')year
  IF(mpi_rank==0) then
    IF (outdiag==1) THEN
      OPEN(idiagp,file='diagprod.'//monthstr//'.'//yearstr//'.txt',form='formatted')
      OPEN(idiagl,file='diagloss.'//monthstr//'.'//yearstr//'.txt',form='formatted')
      OPEN(idiagd,file='diagdepo.'//monthstr//'.'//yearstr//'.txt',form='formatted')
      OPEN(idiage,file='diagemis.'//monthstr//'.'//yearstr//'.txt',form='formatted')
      OPEN(idiagb,file='diagburd.'//monthstr//'.'//yearstr//'.txt',form='formatted')
    ENDIF
  ENDIF

  ! ALLOCATE INPUTS
  ALLOCATE(sflx(2,iip1,jjp1,nbounds))
  ALLOCATE(sflux(iip1,jjp1,nbounds,iqmax))
  ALLOCATE(sflux_tl(iip1,jjp1,nbounds,iqmax))
  ALLOCATE(sflux_ad(iip1,jjp1,nbounds,iqmax))

  ! SACS variables
  ALLOCATE(scale(iip1,jjp1,nbounds, iprescrmax))
  ALLOCATE(scale_tl(iip1,jjp1,nbounds, iprescrmax))
  ALLOCATE(scale_ad(iip1,jjp1,nbounds, iprescrmax))
  ALLOCATE(pscale(iip1,jjp1,nbounds, iprodmax))
  ALLOCATE(pscale_tl(iip1,jjp1,nbounds, iprodmax))
  ALLOCATE(pscale_ad(iip1,jjp1,nbounds, iprodmax))

  ALLOCATE(refprescr_glo(ngridmx,llm,iprescrmax))
  ALLOCATE(refprod_glo(ngridmx,llm,iprodmax))
  ALLOCATE(refjrates_glo(ngridmx,llm,ijratesmax))
  
  ALLOCATE(depvel_glo(ngridmx, idepmax))

  ! Read scale factors, surface fluxes and satellite weighting functions
  ! as defined in PYVAR
  sflux(:,:,:,:) = 0.
  sflux_tl(:,:,:,:) = 0.
  sflux_ad(:,:,:,:) = 0.
  scale(:,:,:,:) = 1.
  scale_tl(:,:,:,:) = 0.
  scale_ad=0.
  pscale(:,:,:,:) = 1.
  pscale_ad=0.
  pscale_tl(:,:,:,:) = 0.
  qscale(:,:,:) = 1.
  qscale_tl(:,:,:) = 0.
  qscale_ad=0.

  ! LOOP ON ACTIVE SPECIES = flux + initial conditions
  ! IF mod.bin species, initializes fluxes from it
  DO iq = 1, iqmax
    INQUIRE(FILE='mod_'//trim(species(iq)%name)//'.bin', EXIST=file_exists)
    IF (file_exists) THEN
      sflx(:,:,:,:) = 0.
      OPEN(1,file='mod_'//trim(species(iq)%name)//'.bin',form='unformatted',access='direct', &
                recl=8*iip1*jjp1*nbounds*2)
      READ(1,rec=1) sflx(:,:,:,:)
      sflux(:,:,:,iq) = sflx(1,:,:,:)
      sflux_tl(:,:,:,iq) = sflx(2,:,:,:)
      CLOSE(1)
      
    END IF
    
!    ! 2D map of scaling factors on initial state
!    ! Note: consider integrating 3D scaling factors if necessary
!    INQUIRE(FILE='init_'//trim(species(iq)%name)//'.bin', EXIST=file_exists)
!    IF (file_exists) THEN
!      OPEN(1,file='init_'//trim(species(iq)%name)//'.bin',form='unformatted',action='read')
!      DO i = 1, iip1*jjp1
!        READ(1) qscale(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,iq), &
!                qscale_tl(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,iq)
!      ENDDO
!      CLOSE(1)
!    END IF
  ENDDO

  ! LOOP ON PRESCRIBED SPECIES = scaling factors
  IF(do_chemistry.and.iprescrmax > 0)THEN
    DO iq = 1, iprescrmax
      INQUIRE(FILE='mod_scale_'//trim(prescr_species(iq)%name)//'.bin', EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(1,file='mod_scale_'//trim(prescr_species(iq)%name)//'.bin',form='unformatted',action='read')
        j = 0
        DO iday = 1, ndayloc + periodflux, periodflux
          j = j + 1
          DO i = 1, iip1*jjp1
            READ(1) scale(   i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq), &
                    scale_tl(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq)
            ! remove negative values (!! not in AD nor in TL!!)
            scale(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq) = &
                    MAX( 0., scale(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq) )
          ENDDO
        END DO
        CLOSE(1)
      END IF
    END DO
  END IF


  ! LOOP ON PROD/LOSS SPECIES = scaling factors
  IF(do_chemistry.and.iprodmax > 0)THEN
    DO iq = 1, iprodmax
      INQUIRE(FILE='mod_prodscale_'//trim(prodloss_species(iq)%name)//'.bin', EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(1,file='mod_prodscale_'//trim(prodloss_species(iq)%name)//'.bin',form='unformatted',action='read')
        j = 0
        DO iday = 1, ndayloc + periodflux, periodflux
          j = j + 1
          DO i = 1, iip1*jjp1
            READ(1) pscale(   i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq), &
                    pscale_tl(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq)
            ! remove negative values (!! not in AD nor in TL!!)
            pscale(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq) = &
                    MAX( 0., pscale(i-((i-1)/iip1)*iip1,(i-1)/iip1+1,j, iq) )
          ENDDO
        END DO
        CLOSE(1)
      END IF
    END DO
  END IF

  ! Reading initial conditions if needed
  time_0=0.
  q(:,:,:,:) = 0.
  q_ad(:,:,:,:) = 0.
  q_tl(:,:,:,:) = 0.

  IF (readstart) THEN
    IF (fwd) THEN
      CALL dynstate0("start.nc",q,time_0)
      
      INQUIRE(FILE='start_tl.nc', EXIST=file_exists)
      IF (file_exists) THEN
        CALL dynstate0("start_tl.nc",q_tl,time_0)
      ELSE
        INQUIRE(FILE='start_tl.bin', EXIST=file_exists)
        IF (file_exists) THEN
          OPEN(1,file='start_tl.bin',form='unformatted',status='unknown',action='read')
          READ(1,END=1) q_tl
1         CLOSE(1)
        END IF
      END IF
      
    ELSE
      CALL dynstate0("start.nc",q_ad,time_0)
    ENDIF
    PRINT*,'READING start'
    PRINT*,'YEAR_INI= ',year_ini
    PRINT*,'DAY_INI= ',day_ini
  ELSE
    day_ini=0
    year_ini=0
  ENDIF

  IF (fwd) THEN

    tangent = .TRUE.
    IF (.NOT. ( ANY(sflux_tl /= 0.) &
            .OR. ANY(qscale_tl /= 0.) &
            .OR. ANY(q_tl /= 0.) &
            .OR. ANY(scale_tl /= 0.) &
            .OR. ANY(pscale_tl /= 0.)))THEN
      tangent = .FALSE.
    END IF

  ELSE
    tangent = .FALSE.
    sflux_ad(:,:,:,:) = 0.
    qscale_ad(:,:,:) = 0.
    scale_ad(:,:,:,:) = 0.
    IF (outwfunc) THEN
      OPEN(1,file='wfunc.bin',form='unformatted',action='read')
      READ(1) wfunc_int
      CLOSE(1)
    ENDIF
  ENDIF
  
  !
  !-- READ DEPOSITION VELOCITY
  ! ngridmx,2
  IF(do_chemistry .AND. idepmax > 0)THEN

    depvel_glo(:,:)=0.
    DO iq = 1, idepmax
      INQUIRE(FILE='dep_'//trim(dep_species(iq)%name)//'.nc', EXIST=file_exists)
      IF (file_exists) THEN
        CALL readsdepvel(trim(dep_species(iq)%name), &
                'dep_'//trim(dep_species(iq)%name)//'.nc',depvel)

        ! Scaling CH2O deposition
        if (dep_species(iq)%name == 'CH2O') THEN
          depvel = 1.4 * depvel
        end if

        ! Deposition files are assumed to be monthly files
        depvel_glo(:,iq)=depvel(:,month)
      END IF
    END DO
  END IF
  
  !   Re-defining the splitting step
  IF (MOD(nsplit_dyn,nsplit_phy)==0) THEN
    nsplit=nsplit_phy
  ELSEIF (MOD(nsplit_phy,nsplit_dyn)==0) THEN
    nsplit=nsplit_dyn
  ELSE
    STOP'One of the two time splits must divide the other'
  ENDIF
  nsplit_dyn=nsplit_dyn/nsplit
  nsplit_phy=nsplit_phy/nsplit
  PRINT*,'nsplit,nsplit_dyn,nsplit_phy', nsplit,nsplit_dyn,nsplit_phy

  ! Writing first step in trajq
  IF (fwd) THEN
    IF (mpi_rank==0) THEN
      DO i=1, iqmax
        WRITE(itrajq(i),rec=1) q(:,:,:,i)
      END DO
    END IF
!    DO iq=1,iqmax
!      DO i = 1, llm
!        q_tl(:,:,i,iq) = q_tl(:,:,i,iq)*qscale(:,:,iq) + q(:,:,i,iq)*qscale_tl(:,:,iq)
!        q(:,:,i,iq) = q(:,:,i,iq) * qscale(:,:,iq)
!      ENDDO
!    ENDDO
  ENDIF

  ! Initialize constants
  g = 9.8
  rad = 6400000
  omeg = 7.272205e-05
  kappa = 0.285716
  daysec = 86400
  cpp = 1004.70885
  preff = 101325.
  pa= 50000.
  ucov = 0.
  vcov = 0.
  teta = 0.
  masse = 0.
  ps = 0.
  phis = 0.

  CALL dynstate1("vcoord.nc")

  CALL iniconst
  CALL inigeom

  day_end=day_ini+nday
  IF (fwd) THEN
    day0=1
  ELSE
    day0=nday
  ENDIF

  PRINT*, 'Day0 ', day0
  
  !  Reading headers from flux files
  PRINT*,'Reading headers from flux files (readfluxnc_p)'
  irec = 0
  CALL readfluxnc_p(irec,masse,pbaru,pbarv,w,teta,phi, &
       nrec,areafi,pphis,t,zmfu, zmfd, zde_u,zen_d, coefkz, &
       dake, mpke, phike, updke, dndke, wghtke, &
       entr_therm,fm_therm,&
       yu1, yv1, ftsol, pctsrf, phis)

  ! Initializes chemistry
  PRINT*, 'Initializing chemistry (readchem)'
  IF(do_chemistry)THEN
    CALL readchem(irec, &
            refprescr_glo,refprod_glo, refjrates_glo, &
            temp_glo, pmid_glo, convoh)
  END IF

  ! on recalcule le nombre de pas de temps dans une journee
  !   istdyn est la frequence de stokage en "pas de temps dynamique" dans
  !   le modele on-line. dtvr est le pas de temps du meme modele en s.
  ! iday-step est le nombre de lectures par jour (grand pas de temps).

  iday_step=daysec/(istdyn*dtvr)
  PRINT*,'iday_step=',iday_step
  !
  CALL iniconst

  ! Computing time step
  dtvr    = daysec/FLOAT(iday_step)
  dtphys  = dtvr/(nsplit*nsplit_phy) ! Should be after iniconst!
  PRINT*,'Physics time step ', dtphys
  PRINT *,'Dynamics time steps', dtvr
  !
  CALL gr_dyn_fi_p(1,iip1,jjp1,klon,area,areafi)
  !
  CALL suphec
  CALL inifilr
  PRINT*,'Physics time step ',dtphys

  !   Initializing output file
  IF (mpi_rank==0) CALL dynrestart0("restart.nc",day_end,year_ini,phis)

  !
  IF (.NOT. fwd) THEN
    indexq = nday*iday_step + 1
  ELSE
    indexq = 1
  ENDIF

  CALL timeloop(nobs, time_0, iday_step, fwd, nsplit, nsplit_dyn, &
          nsplit_phy, day0, nrec, areafi, pphis, phis, q, q_ad, ps, &
          q_tl, phi, w, coefkz, tabobs, nsec, &
          wfunc_int, indexq, itrajq,  &
          periodflux, nbounds, bounddays, sflux_tl, sflux, sflux_ad, masse, &
          pbaru,pbarv,teta,t,zmfu, zmfd, zde_u, zen_d, yu1, yv1, ftsol, pctsrf, &
          vcov, ucov, tangent, &
          dake,mpke,phike,updke,dndke,wghtke,entr_therm,fm_therm, &
          refprescr_glo,refprod_glo,temp_glo,pmid_glo,refjrates_glo, &
          scale, scale_tl, scale_ad, pscale, pscale_tl, pscale_ad, &
          depvel_glo,month,year,idiagd,idiagl,idiagp,idiage,idiagb,outdiag,outwfunc,convoh)

  IF (mpi_rank==0) WRITE(*,*)'Timeloop is over. Gathering fields'

  CALL gather_field(ucov,ip1jmp1,llm,0)
  CALL gather_field(vcov,ip1jm,llm,0)
  CALL gather_field(teta,ip1jmp1,llm,0)
  CALL gather_field(ps,ip1jmp1,1,0)
  CALL gather_field(masse,ip1jmp1,llm,0)

  IF (fwd) THEN
    CALL gather_field(q,ip1jmp1,llm*iqmax,0)
    CALL gather_field(q_tl,ip1jmp1,llm*iqmax,0)
    IF (mpi_rank==0) CALL dynrestart1("restart.nc",0.0,vcov,ucov,teta,q,masse,ps)
  ELSE
    IF(iqmax > 0) CALL gather_field(q_ad,ip1jmp1,llm*iqmax,0)
    IF(iprodmax > 0) CALL gather_field(pscale_ad,ip1jmp1,iprodmax*nbounds,0)
    IF(iprescrmax > 0) CALL gather_field(scale_ad,ip1jmp1,iprescrmax*nbounds,0)
    IF(iqmax > 0) CALL gather_field(sflux_ad,ip1jmp1,iqmax*nbounds,0)
    IF (mpi_rank==0)  CALL dynrestart1("restart.nc",0.0,vcov,ucov,teta,q_ad,masse,ps)
  ENDIF

  !----------------------------------------------------------------------
  !
  ! Write output for Python
  !
  IF (mpi_rank==0) THEN
    WRITE(*,*)'WRITE OUTPUT FOR PYTHON'

    IF (fwd) THEN

      IF (nobs > 0) THEN
        OPEN(1,file='obs_out.bin',form='unformatted',&
                access='direct', recl=8*4*nobs)
        WRITE(1, rec=1) tabobs(:,7:10)
        CLOSE(1)
      END IF
      IF (tangent) THEN
        OPEN(1,file='restart_tl.bin',form='unformatted')
        WRITE(1) q_tl
        CLOSE(1)
      ENDIF

    ELSE

      DO i = 1, iqmax
        READ(itrajq(i),rec=1) q(:,:,:,i)
      END DO

      ! WRITING ADJOINT OUTPUTS
      IF (mpi_rank==0) THEN
        ! fluxes
        DO iq = 1, iqmax
          OPEN(1666,file='mod_fluxes_'//trim(species(iq)%name)//'_out.bin',&
                  form='unformatted',access='direct', recl=8*iip1*jjp1)
          j = 0
          DO iday = 1, ndayloc, periodflux
            j = j + 1
            WRITE(1666, rec=j) sflux_ad(:,:,j,iq)
            WRITE(*,*) 'sflux', trim(species(iq)%name), j, MAXVAL(sflux_ad(:,:,j,iq))
          END DO
          CLOSE(1666)
        END DO

        ! initial condition sensitivity
        DO iq = 1, iqmax
          OPEN(1667,file='mod_init_'//trim(species(iq)%name)//'_out.bin',&
                  form='unformatted',access='direct', recl=8*iip1*jjp1*llm)
          WRITE(1667, rec=1) q_ad(:,:, :, iq)
          WRITE(*,*) 'init', trim(species(iq)%name), MAXVAL(q_ad(:,:,:, iq))
          CLOSE(1667)
        END DO

        ! prescribed fields
        DO iq = 1, iprescrmax
          OPEN(1668,file='mod_scale_'//trim(prescr_species(iq)%name)//'_out.bin',&
                  form='unformatted',access='direct', recl=8*iip1*jjp1)
          j = 0
          DO iday = 1, ndayloc, periodflux
            j = j + 1
            WRITE(1668, rec=j) scale_ad(:,:,j,iq)
            WRITE(*,*) 'scale', trim(prescr_species(iq)%name), j, MAXVAL(scale_ad(:,:,j,iq))
          END DO
          CLOSE(1668)
        END DO

        ! prod/loss fields
        DO iq = 1, iprodmax
          OPEN(1669,file='mod_prodscale_'//trim(prodloss_species(iq)%name)//'_out.bin',&
                  form='unformatted',access='direct', recl=8*iip1*jjp1)
          j = 0
          DO iday = 1, ndayloc, periodflux
            j = j + 1
            WRITE(1669, rec=j) pscale_ad(:,:,j,iq)
            WRITE(*,*) 'pscale', trim(prodloss_species(iq)%name), j, MAXVAL(pscale_ad(:,:,j,iq))
          END DO
          CLOSE(1669)
        END DO

      ENDIF
    ENDIF

    ! Closing trajectory files
    DO iq=1, iqmax
      CLOSE(itrajq(iq))
    END DO
    
    ! Closing diagnostic files
    IF(mpi_rank==0) then
      IF (outdiag==1) then
        CLOSE(idiagp)
        CLOSE(idiagl)
        CLOSE(idiagd)
        CLOSE(idiage)
        CLOSE(idiagb)
      ENDIF
    ENDIF
  
  ENDIF

  if (mpi_rank==0)THEN
    WRITE(*,*)'DEALLOCATE VARIABLES'
    DEALLOCATE( sflux )
    DEALLOCATE( sflux_tl )
    DEALLOCATE( sflux_ad )
    DEALLOCATE( scale )
    DEALLOCATE( scale_tl )
    DEALLOCATE( scale_ad )
    DEALLOCATE( q )
    DEALLOCATE( q_tl )
    DEALLOCATE( q_ad )
    DEALLOCATE( vcov )
    DEALLOCATE( ucov )
    DEALLOCATE( phi )
    DEALLOCATE( w )
    DEALLOCATE( coefkz )
    DEALLOCATE( qscale )
    DEALLOCATE( qscale_ad )
    DEALLOCATE( qscale_tl )
    DEALLOCATE( pscale )
    DEALLOCATE( pscale_ad )
    DEALLOCATE( pscale_tl )
    DEALLOCATE( tabobs )
    DEALLOCATE( bounddays )
    DEALLOCATE( zmfd )
    DEALLOCATE( zen_d )
    DEALLOCATE( zmfu )
    DEALLOCATE( zde_u )
    DEALLOCATE( t )
    DEALLOCATE( dake )
    DEALLOCATE( mpke )
    DEALLOCATE( phike )
    DEALLOCATE( updke )
    DEALLOCATE( dndke )
    DEALLOCATE( wghtke )
    DEALLOCATE( fm_therm )
    DEALLOCATE( entr_therm )
    DEALLOCATE( pphis )
    DEALLOCATE( areafi )
    DEALLOCATE( yu1 )
    DEALLOCATE( yv1 )
    
    ! End properly
    OPEN(1, file='all_good')
    CLOSE(1)
    
  end if
  
  WRITE(*,*)'FINALIZING PARALLEL'
  CALL finalize_parallel
  
  STOP
END PROGRAM dispersion
