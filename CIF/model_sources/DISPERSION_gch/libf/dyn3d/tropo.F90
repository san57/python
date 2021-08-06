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

SUBROUTINE tropo(temp, nlon, nlat, nlev, pres, plimu, pliml, plimlex, dofill, tp, tperr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! determination of tropopause height from gridded temperature data
  !
  ! reference: Reichler, T., M. Dameris, and R. Sausen (2003)
  !
  ! input:    temp(nlon,nlat,nlev)    3D-temperature field
  !           nlon                    grid points in x
  !           nlat                    grid points in y
  !           pres(nlev)              pressure levels in hPa
  !           plimu                   upper limit for tropopause pressure in Pa, usually 45000.
  !           pliml                   lower limit for tropopause pressure in Pa, usually 7500.
  !           plimlex                 lower limit in extratropics, usually same as pliml, i.e., 7500.
  !           dofill                  fill undefined values with neighboring points if .true.
  !
  ! output:   tp(nlon, nlat)          tropopause pressure in Pa
  !           tperr                   # of undetermined values
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  
  INTEGER,INTENT(in)                        :: nlon, nlat, nlev
  REAL,INTENT(in),DIMENSION(nlon,nlat,nlev) :: temp
  REAL,INTENT(in),DIMENSION(nlev)           :: pres
  REAL, INTENT(in)                          :: plimu, pliml, plimlex
  LOGICAL, INTENT(in)                       :: dofill
  REAL,INTENT(out),DIMENSION(nlon,nlat)     :: tp
  INTEGER,INTENT(out)                       :: tperr
  
  INTEGER                                   :: i, invert, ifil
  INTEGER                                   :: lon, lat
  REAL,DIMENSION(nlev)                      :: t
  REAL,DIMENSION(nlev)                      :: p
  REAL                                      :: trp
  
  REAL, PARAMETER                           :: gamma=-0.002 ! K/m
  
  ! check vertical orientation of data
  IF (pres(1) > pres(2)) THEN
      invert=1
      DO i=1,nlev
        p(i)=pres(nlev+1-i)*100.  ! hPa > Pa
      ENDDO
  ELSE
      invert=0
      DO i=1,nlev
        p(i)=pres(i)*100.         ! hPa > Pa
      ENDDO
  ENDIF
  
  tperr = 0
  DO lon=1,nlon
    DO lat=1,nlat
      IF (invert==1) THEN
          DO i=1,nlev
            t(i)=temp(lon,lat,nlev+1-i)
          ENDDO
      ELSE
          DO i=1,nlev
            t(i)=temp(lon,lat,i)
          ENDDO
      ENDIF
      CALL twmo(nlev, t, p, plimu, pliml, gamma, trp)
      IF (lat.LT..15*nlat.AND.trp.LT.plimlex) trp=-99.
      IF (lat.GT..85*nlat.AND.trp.LT.plimlex) trp=-99.
      tp(lon,lat)=trp
      IF(trp.LT..0) THEN
          tperr = tperr+1    
      ENDIF
    END DO
  END DO
  
  ! fill holes
  IF (dofill) THEN
      CALL fill(tp, nlon, nlat, ifil)
      IF (ifil/=tperr) THEN
          PRINT*, 'Inconsistent'
          STOP
      ENDIF
  ENDIF
  
  RETURN
END SUBROUTINE tropo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! twmo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE twmo(level, t, p, plimu, pliml, gamma, trp)
  
  IMPLICIT NONE
  INTEGER,INTENT(in)                  :: level
  REAL,INTENT(in),DIMENSION(level)    :: t, p
  REAL,INTENT(in)                     :: plimu, pliml, gamma
  REAL,INTENT(out)                    :: trp
  
  REAL,PARAMETER                      :: kap=0.286
  REAL,PARAMETER                      :: faktor = -9.81/287.0
  REAL,PARAMETER                      :: deltaz = 2000.0
  REAL,PARAMETER                      :: ka1=kap-1.
  
  REAL                                :: pmk, pm, a, b, tm, dtdp, dtdz
  REAL                                :: ag, bg, ptph
  REAL                                :: pm0, pmk0, dtdz0
  REAL                                :: p2km, asum, aquer
  REAL                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
  INTEGER                             :: icount, jj
  INTEGER                             :: j
      
  trp=-99.0                           ! negative means not valid
  DO j=level,2,-1
    
    ! dt/dz
    pmk= .5 * (p(j-1)**kap+p(j)**kap)
    pm = pmk**(1/kap)              
    a = (t(j-1)-t(j))/(p(j-1)**kap-p(j)**kap)
    b = t(j)-(a*p(j)**kap)
    tm = a * pmk + b              
    dtdp = a * kap * (pm**ka1)
    dtdz = faktor*dtdp*pm/tm
    
    ! dt/dz valid?
    IF (j==level)    go to 999     ! no, start level, initialize first
    IF (dtdz<=gamma) go to 999     ! no, dt/dz < -2 K/km
    IF (pm>plimu)   go to 999     ! no, too low
    
    ! dtdz is valid, calculate tropopause pressure
    IF (dtdz0<gamma) THEN
        ag = (dtdz-dtdz0) / (pmk-pmk0)    
        bg = dtdz0 - (ag * pmk0)         
        ptph = EXP(LOG((gamma-bg)/ag)/kap)
    ELSE
        ptph = pm
    ENDIF
    
    IF (ptph<pliml) go to 999    
    IF (ptph>plimu) go to 999          
    
    ! 2nd test: dtdz above 2 km must not exceed gamma
    p2km = ptph + deltaz*(pm/tm)*faktor          ! p at ptph + 2km
    asum = 0.0                                   ! dtdz above
    icount = 0                                   ! number of levels above
    
    ! test until apm < p2km
    DO jj=j,2,-1
      
      pmk2 = .5 * (p(jj-1)**kap+p(jj)**kap)    ! p mean ^kappa
      pm2 = pmk2**(1/kap)                      ! p mean
      IF(pm2>ptph) go to 110                ! doesn't happen
      IF(pm2<p2km) go to 888                ! ptropo is valid
      
      a2 = (t(jj-1)-t(jj))                     ! a
      a2 = a2/(p(jj-1)**kap-p(jj)**kap)
      b2 = t(jj)-(a2*p(jj)**kap)               ! b
      tm2 = a2 * pmk2 + b2                     ! T mean
      dtdp2 = a2 * kap * (pm2**(kap-1))        ! dt/dp
      dtdz2 = faktor*dtdp2*pm2/tm2
      asum = asum+dtdz2
      icount = icount+1
      aquer = asum/float(icount)               ! dt/dz mean
      
      ! discard ptropo ?
      IF (aquer<=gamma) go to 999           ! dt/dz above < gamma
      
110   CONTINUE
    ENDDO                           ! test next level
    
888 CONTINUE                        ! ptph is valid
    trp = ptph
    RETURN
    
999 CONTINUE                        ! continue search at next higher level
    pm0 = pm
    pmk0 = pmk
    dtdz0  = dtdz
    
  ENDDO
  
  ! no tropopouse found
  RETURN
END SUBROUTINE twmo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE fill(dat, ix, iy, ir)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE parallel
  
  INTEGER, INTENT(in)                 :: ix, iy
  INTEGER, INTENT(out)                :: ir
  REAL, DIMENSION(ix,iy)              :: dat
  REAL, DIMENSION(4)                  :: help
  INTEGER                             :: jx, jy, icount, ipk, ic
  REAL                                :: drop, sum
  
  icount = 0
  DO jy=jj_begin,jj_end
    DO jx=1,ix
      IF (loch(dat(jx,jy))) icount = icount+1
    ENDDO
  ENDDO
  
  IF (icount>(ix*jj_nb)/2) STOP 'ERROR: Too many holes (>50%)'
  ir = icount
  IF (icount==0) RETURN
  
  ipk = 0
10 CONTINUE
  DO jy=jj_begin,jj_end
    DO jx=1,ix
      
      IF(loch(dat(jx,jy))) THEN
          drop = dat(jx,jy)
          
          ! left edge
          IF (jx==1) THEN
              IF (jy==1) THEN
                  help(1) = dat(jx,jy+1)
                  help(2) = dat(jx+1,jy)
                  help(3) = drop
                  help(4) = drop
                  go to 200
              ENDIF
              IF (jy==iy) THEN
                  help(1) = drop
                  help(2) = dat(jx+1,jy)
                  help(3) = dat(jx,jy-1)
                  help(4) = drop
                  go to 200
              ENDIF
              help(1) = dat(jx,jy+1)
              help(2) = dat(jx+1,jy)
              help(3) = dat(jx,jy-1)
              help(4) = drop
              go to 200
          ENDIF
          
          ! right edge
          IF (jx==ix) THEN
              IF (jy==1) THEN
                  help(1) = dat(jx,jy+1)
                  help(2) = drop
                  help(3) = drop
                  help(4) = dat(jx-1,jy)
                  go to 200
              ENDIF
              IF (jy==iy) THEN
                  help(1) = drop
                  help(2) = drop
                  help(3) = dat(jx,jy-1)
                  help(4) = dat(jx-1,jy)
                  go to 200
              ENDIF
              help(1) = dat(jx,jy+1)
              help(2) = drop
              help(3) = dat(jx,jy-1)
              help(4) = dat(jx-1,jy)
              go to 200
          ENDIF
          
          ! bottom edge
          IF (jy==1) THEN
              help(1) = dat(jx,jy+1)
              help(2) = dat(jx+1,jy)
              help(3) = drop
              help(4) = dat(jx-1,jy)
              go to 200
          ENDIF
          
          ! upper edge
          IF(jy==iy) THEN
              help(1) = drop
              help(2) = dat(jx+1,jy)
              help(3) = dat(jx,jy-1)
              help(4) = dat(jx-1,jy)
              go to 200
          ENDIF
          
          ! no edge
          help(1) = dat(jx,jy+1)
          help(2) = dat(jx+1,jy)
          help(3) = dat(jx,jy-1)
          help(4) = dat(jx-1,jy)
          
200       CONTINUE
          
          ic = 0
          sum = 0.0
          DO jj=1,4
            IF(.NOT.loch(help(jj))) THEN
                sum = sum+help(jj)
                ic = ic+1
            ENDIF
          ENDDO
          
          IF (ic>0) THEN
              dat(jx,jy) = sum/float(ic) ! fill with mean of valid
              ipk = ipk+1       ! neighbourpoints
          ENDIF
          
      ENDIF
      IF (ipk >= icount) RETURN ! until all filled
    ENDDO
  ENDDO
  go to 10
  
CONTAINS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL FUNCTION loch(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    REAL, INTENT(in)        :: x
    
    edge = -98.0
    IF (x<edge) THEN
        loch = .TRUE.
    ELSE
        loch = .FALSE.
    ENDIF
    RETURN
  END FUNCTION loch
  
END SUBROUTINE fill

