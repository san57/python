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

SUBROUTINE read_pstoke0(irec, zrec, zkon, zkev, areafi, phisfi, &
                        t, mfu, mfd, de_u, en_d, coefh, dake, mpke, &
                        phike, updke, dndke, wghtke, entr_therm, fm_therm, &
                        pyu1,pyv1,ftsol,psrf)
  
  USE dimphy
  USE mod_phys_lmdz_para
  IMPLICIT NONE
  
  INCLUDE "netcdf.inc"
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "comconst.h"
  INCLUDE "comgeom.h"
  INCLUDE "temps.h"
  INCLUDE "ener.h"
  INCLUDE "logic.h"
  INCLUDE "description.h"
  INCLUDE "serre.h"
  INCLUDE "indicesol.h"
  INCLUDE "control.h"
  INCLUDE "clesph0.h"
  
  INTEGER :: zkon,zkev
  INTEGER :: klev_limconv
  REAL :: phisfi(klon)
  REAL :: phisfi2(iim,jjm+1),areafi2(iim,jjm+1)
  
  REAL :: mfu(klon,klev), mfd(klon,klev)
  REAL :: de_u(klon,klev)
  REAL :: en_d(klon,klev)
  REAL :: coefh(klon,klev)
  REAL :: t(klon,klev)
  REAL :: dake(klon,klev)   !Variables liees au schema de Kerry Emanuel finissent par ke
  REAL :: mpke(klon,klev)
  REAL :: updke(klon,klev)
  REAL :: dndke(klon,klev)
  REAL :: wghtke(klon,klev)
  REAL :: fm_therm(klon,klev)
  REAL :: entr_therm(klon,klev)
  REAL :: phike(klon,klev,klev)

  
  REAL :: mfu2(iim,jj_nb,klev), mfd2(iim,jj_nb,klev)
  REAL :: de_u2(iim,jj_nb,klev)
  REAL :: en_d2(iim,jj_nb,klev)
  REAL :: coefh2(iim,jj_nb,klev)
  REAL :: t2(iim,jj_nb,klev)
  REAL :: dake2(iim,jj_nb,klev)
  REAL :: mpke2(iim,jj_nb,klev)
  REAL :: updke2(iim,jj_nb,klev)
  REAL :: dndke2(iim,jj_nb,klev)
  REAL :: wghtke2(iim,jj_nb,klev)
  REAL :: fm_therm2(iim,jj_nb,klev)
  REAL :: entr_therm2(iim,jj_nb,klev)
  REAL :: phike2(iim,jj_nb,klev,klev)

  REAL :: zerovector(iim,jj_nb,klev)

  
  REAL :: pl(klev)
  INTEGER :: irec
  INTEGER :: xid,yid,zid,tid
  INTEGER :: zrec,zim,zjm
  INTEGER :: ncrec,nckev,ncim,ncjm
  SAVE :: zim,zjm
  
  REAL :: areafi(klon)
  CHARACTER(len=NF_MAX_NAME) ::  namedim

  REAL :: pyu1(klon), pyv1(klon)
  REAL :: pyu12(iim,jj_nb), pyv12(iim,jj_nb)
  REAL :: ftsol(klon,nbsrf)
  REAL :: psrf(klon,nbsrf)
  REAL :: ftsol1(klon),ftsol2(klon),ftsol3(klon),ftsol4(klon)
  REAL :: psrf1(klon),psrf2(klon),psrf3(klon),psrf4(klon)
  REAL :: ftsol12(iim,jj_nb),ftsol22(iim,jj_nb),ftsol32(iim,jj_nb),ftsol42(iim,jj_nb)
  REAL :: psrf12(iim,jj_nb),psrf22(iim,jj_nb),psrf32(iim,jj_nb), psrf42(iim,jj_nb)

  INTEGER,SAVE :: ncidp
  INTEGER,SAVE :: varidmfu, varidmfd, varidps, variddeu
  INTEGER,SAVE :: varidda, varidmp, varidupd, variddnd, varidwght
  INTEGER,SAVE :: varidfmt,varident
  INTEGER,ALLOCATABLE,SAVE :: varidphi(:)
  INTEGER,SAVE :: varidt
  INTEGER,SAVE :: varidend,varidch
  INTEGER,SAVE :: varidyu1,varidyv1,varidpl,varidai
  INTEGER,SAVE :: varidfts1,varidfts2,varidfts3,varidfts4
  INTEGER,SAVE :: varidpsr1,varidpsr2,varidpsr3,varidpsr4
  
  INTEGER :: i, start(4),COUNT(4),status
  INTEGER :: k
  INTEGER :: i_iim,i_jjnb,iklev
  REAL :: rcode
  CHARACTER(len=12) :: nvar

  !Definition du vecteur nul
  zerovector(:,:,:)=0.
  !Defintion du niveau limite pour la convection de KE
  klev_limconv=25



  IF (.not. allocated(varidphi)) allocate(varidphi(klev))


  ! ---------------------------------------------
  !   Initialisation de la lecture des fichiers
  ! ---------------------------------------------
  
  IF (irec == 0) THEN
      
      ncidp=NCOPN('phystoke.nc',NCNOWRIT,rcode)
      
      varidps=NCVID(ncidp,'phis',rcode)
      
      varidpl=NCVID(ncidp,'sig_s',rcode)
      
      varidai=NCVID(ncidp,'aire',rcode)
      
      varidmfu=NCVID(ncidp,'mfu',rcode)
      
      varidt=NCVID(ncidp,'t',rcode)
  
      IF (iflag_con.EQ.2) THEN       !Tiedkte scheme

        varidmfu=NCVID(ncidp,'mfu',rcode)
    
        varidmfd=NCVID(ncidp,'mfd',rcode)
   
        variddeu=NCVID(ncidp,'de_u',rcode)
        
        varidend=NCVID(ncidp,'en_d',rcode)
    
      ENDIF
      
      varidch=NCVID(ncidp,'coefh',rcode)            !BL coeff

      IF (iflag_con.EQ.3 .OR. iflag_con.EQ.30) THEN         !KE scheme
     
        varidda=NCVID(ncidp,'da',rcode)

        varidmp=NCVID(ncidp,'mp',rcode)

        varidupd=NCVID(ncidp,'upwd',rcode)
      
        variddnd=NCVID(ncidp,'dnwd',rcode)

        varidwght=NCVID(ncidp,'wght',rcode)
  
!traitement special pour la variable phi (convection KE)
        DO k=1,klev_limconv
         IF (k<10) THEN
           WRITE(nvar,'(i1)') k
         ELSE IF (k<100) THEN
           WRITE(nvar,'(i2)') k
         ELSE
           WRITE(nvar,'(i3)') k
         END IF
         nvar='phi_lev'//trim(nvar)
         varidphi(k)=NCVID(ncidp,nvar,rcode)
        END DO
 
      ENDIF    !end loop on KE scheme
   
      IF (iflag_con .EQ. 30) THEN

       varidfmt=NCVID(ncidp,'fm_th',rcode)

       varident=NCVID(ncidp,'en_th',rcode)
      ENDIF
      
      varidch=NCVID(ncidp,'coefh',rcode)
      
      varidyu1=NCVID(ncidp,'pyu1',rcode)
      
      varidyv1=NCVID(ncidp,'pyv1',rcode)

      varidfts1=NCVID(ncidp,'ftsol1',rcode)

      varidfts2=NCVID(ncidp,'ftsol2',rcode)

      varidfts3=NCVID(ncidp,'ftsol3',rcode)

      varidfts4=NCVID(ncidp,'ftsol4',rcode)

      varidpsr1=NCVID(ncidp,'psrf1',rcode)

      varidpsr2=NCVID(ncidp,'psrf2',rcode)

      varidpsr3=NCVID(ncidp,'psrf3',rcode)

      varidpsr4=NCVID(ncidp,'psrf4',rcode)

      !       ID pour les dimensions
      
      status = nf_inq_dimid(ncidp,'lat',yid)
      status = nf_inq_dimid(ncidp,'lon',xid)
      status = nf_inq_dimid(ncidp,'sig_s',zid)
      status = nf_inq_dimid(ncidp,'time_counter',tid)
      
      !       lecture des dimensions
      
      status = nf_inq_dim(ncidp,yid,namedim,ncjm)
      status = nf_inq_dim(ncidp,xid,namedim,ncim)
      status = nf_inq_dim(ncidp,zid,namedim,nckev)
      status = nf_inq_dim(ncidp,tid,namedim,ncrec)
      
      zrec=ncrec
      zkev=nckev
      zim=ncim
      zjm=ncjm
      
      zkon=zim*(zjm-2)+2
      
      WRITE(*,*) 'read_pstoke : zrec = ', zrec
      WRITE(*,*) 'read_pstoke : kev = ', zkev
      WRITE(*,*) 'read_pstoke : zim = ', zim 
      WRITE(*,*) 'read_pstoke : zjm = ', zjm
      WRITE(*,*) 'read_pstoke : kon = ', zkon
      
      !       niveaux de pression
      
      status=NF_GET_VARA_REAL(ncidp,varidpl,1,klev,pl)
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF
      
      !       lecture de aire et phis
      
      start(1)=1
      start(2)=jj_begin
      start(3)=1
      start(4)=0
      
      COUNT(1)=iim
      COUNT(2)=jj_nb
      COUNT(3)=1
      COUNT(4)=0
      
      !       
      !       phis
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidps,start,count,phisfi2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidps,start,count,phisfi2)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF
      
      CALL grid2Dto1D_mpi(phisfi2,phisfi)
      !       aire
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidai,start,count,areafi2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidai,start,count,areafi2)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF
      
      CALL grid2Dto1D_mpi(areafi2,areafi)
  ELSE
      
      
      !       ---------------------
      !       lecture des champs
      !       ---------------------
      
      start(1)=1
      start(2)=jj_begin
      start(3)=1
      start(4)=irec
      
      COUNT(1)=iim
      COUNT(2)=jj_nb
      COUNT(3)=klev
      COUNT(4)=1
      
      !       abder t
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidt,start,count,t2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidt,start,count,t2)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF
      
      CALL grid2Dto1D_mpi(t2,t)
      
      IF(iflag_con.EQ.2) THEN       !Tiedtke scheme
        !        mfu
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidmfu,start,count,mfu2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidmfu,start,count,mfu2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF
        CALL grid2Dto1D_mpi(mfu2,mfu)
      
        !       mfd
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidmfd,start,count,mfd2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidmfd,start,count,mfd2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF
        CALL grid2Dto1D_mpi(mfd2,mfd)
      
        !       de_u
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,variddeu,start,count,de_u2)
#else
        status=NF_GET_VARA_REAL(ncidp,variddeu,start,count,de_u2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF
      
        CALL grid2Dto1D_mpi(de_u2,de_u)
      
        !       en_d
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidend,start,count,en_d2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidend,start,count,en_d2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF
      
        CALL grid2Dto1D_mpi(en_d2,en_d)
      ENDIF
      
      !       coefh 
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidch,start,count,coefh2)
#else
      status=NF_GET_VARA_REAL(ncidp,varidch,start,count,coefh2)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF
      
      CALL grid2Dto1D_mpi(coefh2,coefh)

      IF (iflag_con.EQ.3 .OR. iflag_con.EQ.30) THEN    !KE convection scheme
        !       da (convection Kerry Emanuel)
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidda,start,count,dake2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidda,start,count,dake2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF

        CALL grid2Dto1D_mpi(dake2,dake)

        !       mp (convection Kerry Emanuel)
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidmp,start,count,mpke2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidmp,start,count,mpke2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF

        CALL grid2Dto1D_mpi(mpke2,mpke)

        !       upd (convection Kerry Emanuel)
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidupd,start,count,updke2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidupd,start,count,updke2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF

        CALL grid2Dto1D_mpi(updke2,updke)

        !       dnd (convection Kerry Emanuel)
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,variddnd,start,count,dndke2)
#else
        status=NF_GET_VARA_REAL(ncidp,variddnd,start,count,dndke2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF
      
        CALL grid2Dto1D_mpi(dndke2,dndke)

        !       wght (convection Kerry Emanuel)
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidwght,start,count,wghtke2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidwght,start,count,wghtke2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF

        CALL grid2Dto1D_mpi(wghtke2,wghtke)

       !       phi (convection Kerry Emanuel)
       DO k=1,klev_limconv              !niveaux où il y a de la convection
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidphi(k),start,count,phike2(:,:,:,k))
#else
        status=NF_GET_VARA_REAL(ncidp,varidphi(k),start,count,phike2(:,:,:,k))
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF

        CALL grid2Dto1D_mpi(phike2(:,:,:,k),phike(:,:,k))
       END DO
       DO k=klev_limconv+1,klev
        CALL grid2Dto1D_mpi(zerovector(:,:,:),phike(:,:,k))
        phike(:,:,k)=phike(:,:,k)
       END DO
      ENDIF

      IF (iflag_con .EQ. 3) THEN
 
       !       entr_therm
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varident,start,count,entr_therm2)
#else
        status=NF_GET_VARA_REAL(ncidp,varident,start,count,entr_therm2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF

        CALL grid2Dto1D_mpi(entr_therm2,entr_therm)

       !       fm_therm
#ifdef NC_DOUBLE
        status = NF_GET_VARA_DOUBLE(ncidp,varidfmt,start,count,fm_therm2)
#else
        status=NF_GET_VARA_REAL(ncidp,varidfmt,start,count,fm_therm2)
#endif
        IF (status /= 0) THEN
            WRITE(0,*) 'Netcdf error', nf_strerror(status)
            STOP
        ENDIF

        CALL grid2Dto1D_mpi(fm_therm2,fm_therm)
      ENDIF

      
      start(3)=irec
      start(4)=0
      COUNT(3)=1
      COUNT(4)=0
      
      !       pyu1
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidyu1,start,count,pyu12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidyu1,start,count,pyu12)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF
      
      CALL grid2Dto1D_mpi(pyu12,pyu1)
      
      !       pyv1
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidyv1,start,count,pyv12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidyv1,start,count,pyv12)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF
      
      CALL grid2Dto1D_mpi(pyv12,pyv1)
	   
      !       ftsol1
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidfts1,start,count,ftsol12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts1,start,count,ftsol12)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      CALL grid2Dto1D_mpi(ftsol12,ftsol1)

      !       ftsol2
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidfts2,start,count,ftsol22)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts2,start,count,ftsol22)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      CALL grid2Dto1D_mpi(ftsol22,ftsol2)

      !       ftsol3
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidfts3,start,count,ftsol32)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts3,start,count,ftsol32)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      CALL grid2Dto1D_mpi(ftsol32,ftsol3)

      !       ftsol4
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidfts4,start,count,ftsol42)
#else
      status=NF_GET_VARA_REAL(ncidp,varidfts4,start,count,ftsol42)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      CALL grid2Dto1D_mpi(ftsol42,ftsol4)

      !       psrf1 
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidpsr1,start,count,psrf12)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr1,start,count,psrf12)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      !       CALL dump2d(iip1-1,jjm+1,psrf12,'PSRF1NC')
      CALL grid2Dto1D_mpi(psrf12,psrf1)
      !       psrf2 
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidpsr2,start,count,psrf22)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr2,start,count,psrf22)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      !       CALL dump2d(iip1-1,jjm+1,psrf22,'PSRF2NC')
      CALL grid2Dto1D_mpi(psrf22,psrf2)

      !       psrf3 
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidpsr3,start,count,psrf32)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr3,start,count,psrf32)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      CALL grid2Dto1D_mpi(psrf32,psrf3)

      !       psrf4 
#ifdef NC_DOUBLE
      status = NF_GET_VARA_DOUBLE(ncidp,varidpsr4,start,count,psrf42)
#else
      status=NF_GET_VARA_REAL(ncidp,varidpsr4,start,count,psrf42)
#endif
      IF (status /= 0) THEN
          WRITE(0,*) 'Netcdf error', nf_strerror(status)
          STOP
      ENDIF

      CALL grid2Dto1D_mpi(psrf42,psrf4)

      DO i = 1,klon
        psrf(i,1) = psrf1(i)
        psrf(i,2) = psrf2(i)
        psrf(i,3) = psrf3(i)
        psrf(i,4) = psrf4(i)

        ftsol(i,1) = ftsol1(i)
        ftsol(i,2) = ftsol2(i)
        ftsol(i,3) = ftsol3(i)
        ftsol(i,4) = ftsol4(i)

      ENDDO
      
  ENDIF
  
  RETURN
  
END SUBROUTINE read_pstoke0

