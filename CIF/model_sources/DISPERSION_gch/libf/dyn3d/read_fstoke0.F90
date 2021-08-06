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

SUBROUTINE read_fstoke0(irec, zrec,zim,zjm,zlm, aready,phis, masse,pbaru,pbarv,w,teta,phi)
  USE parallel
  IMPLICIT NONE
  
  INCLUDE "netcdf.inc"
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "comgeom.h"
  INCLUDE "comvert.h"

  INTEGER :: irec
  INTEGER,PARAMETER ::nlevnc=39
  
  REAL :: pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm)
  REAL :: teta(iip1,jjp1,llm),phis(iip1,jjp1),phi(iip1,jjp1,llm)
  REAL :: masse(iip1,jjp1,llm),w(iip1,jjp1,llm)
  REAL :: aready(iip1,jjp1)
  INTEGER :: ncrec,ncim,ncjm,nclm
  INTEGER :: zrec,zim,zjm,zlm
  INTEGER :: xid,yid,zid,tid

  INTEGER,SAVE :: ncidf,ncidfv
  INTEGER,SAVE :: varidpu,varidpv,varidt,varidw,varidps,varidph,varidai
  INTEGER,SAVE :: varidpl,varidm

  REAL :: pl(nlevnc)

  INTEGER :: start(4),COUNT(4),status

  REAL :: rcode

  CHARACTER(len=NF_MAX_NAME) ::  namedim
  INTEGER :: jje

  !       ---------------------------------------------
  !       Initialisation de la lecture des fichiers
  !       ---------------------------------------------
  
  IF (irec == 0) THEN
   
      ncidf=NCOPN('fluxstoke.nc',NCNOWRIT,rcode)
      
      varidps=NCVID(ncidf,'phis',rcode)
      
      varidpl=NCVID(ncidf,'sig_s',rcode)
      
      varidai=NCVID(ncidf,'aire',rcode)
      
      varidm=NCVID(ncidf,'masse',rcode)
      
      varidpu=NCVID(ncidf,'pbaru',rcode)
      
      varidw=NCVID(ncidf,'w',rcode)
      
      varidt=NCVID(ncidf,'teta',rcode)
      
      varidph=NCVID(ncidf,'phi',rcode)
      
      ncidfv=NCOPN('fluxstokev.nc',NCNOWRIT,rcode)
      
      varidpv=NCVID(ncidfv,'pbarv',rcode)
      
      !       ID pour les dimensions
      
      status = nf_inq_dimid(ncidf,'lat',yid)
      status = nf_inq_dimid(ncidf,'lon',xid)
      status = nf_inq_dimid(ncidf,'sig_s',zid) 
      status = nf_inq_dimid(ncidf,'time_counter',tid)
	   
      !       lecture des dimensions
      
      status = nf_inq_dim(ncidf,yid,namedim,ncjm)
      status = nf_inq_dim(ncidf,xid,namedim,ncim)
      status = nf_inq_dim(ncidf,zid,namedim,nclm)
      status = nf_inq_dim(ncidf,tid,namedim,ncrec)
	   
      zjm=ncjm
      zim=ncim
      zlm=nclm
      zrec=ncrec
      
      WRITE(*,*) 'read_fstoke : zrec = ', zrec
      WRITE(*,*) 'read_fstoke : zlm = ', zlm
      WRITE(*,*) 'read_fstoke : zim = ', zim
      WRITE(*,*) 'read_fstoke : zjm = ', zjm
      
      !       niveaux de pression
	   
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidpl,1,zlm,pl)
#else
      status=NF_GET_VARA_REAL(ncidf,varidpl,1,zlm,pl)
#endif
      
      !       Lecture de phis et aire
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=0
      
      COUNT(1)=zim
      COUNT(2)=zjm
      COUNT(3)=1
      COUNT(4)=0
      
      !       phis
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidps,start,count,phis)
#else
      status=NF_GET_VARA_REAL(ncidf,varidps,start,count,phis)
#endif
      
      !       aire
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidai,start,count,aready)
#else
      status=NF_GET_VARA_REAL(ncidf,varidai,start,count,aready)
#endif
      
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      
      COUNT(1)=zim
      COUNT(2)=zjm
      COUNT(3)=1
      COUNT(4)=1
      
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidm,start,count,masse)
#else
      status=NF_GET_VARA_REAL(ncidf,varidm,start,count,masse)
#endif
      
  ELSE 
      
      !       ---------------------
      !       lecture des champs 
      !       ---------------------
	   
      
      start(1)=1
      start(2)=jj_begin
      start(3)=1
      start(4)=irec
      
      COUNT(1)=iip1
      COUNT(2)=jj_nb
      COUNT(3)=llm
      COUNT(4)=1
      
      !       masse 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidm,start,count,masse(:,jj_begin:jj_end,:))
#else
      status=NF_GET_VARA_REAL(ncidf,varidm,start,count,masse(:,jj_begin:jj_end,:))
#endif
      
      !       pbaru
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidpu,start,count,pbaru(:,jj_begin:jj_end,:))
#else
      status=NF_GET_VARA_REAL(ncidf,varidpu,start,count,pbaru(:,jj_begin:jj_end,:))
#endif
      
      !       w 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidw,start,count,w(:,jj_begin:jj_end,:))
#else
      status=NF_GET_VARA_REAL(ncidf,varidw,start,count,w(:,jj_begin:jj_end,:))
#endif
      
      !       teta 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidt,start,count,teta(:,jj_begin:jj_end,:))
#else
      status=NF_GET_VARA_REAL(ncidf,varidt,start,count,teta(:,jj_begin:jj_end,:))
#endif
	   
      !       phi 
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidf,varidph,start,count,phi(:,jj_begin:jj_end,:))
#else
      status=NF_GET_VARA_REAL(ncidf,varidph,start,count,phi(:,jj_begin:jj_end,:))
#endif
      
      
      IF (south_pole) THEN
          COUNT(2) = jj_nb-1 
          jje=jj_end-1
      ELSE
          jje=jj_end
      ENDIF
      
      !       pbarv
#ifdef NC_DOUBLE
      status=NF_GET_VARA_DOUBLE(ncidfv,varidpv,start,count,pbarv(:,jj_begin:jje,:))
#else
      status=NF_GET_VARA_REAL(ncidfv,varidpv,start,count,pbarv(:,jj_begin:jje,:))
#endif
      
      start(3)=irec
      start(4)=0
      COUNT(2)=jjp1
      COUNT(3)=1
      COUNT(4)=0
      
  ENDIF
  
  RETURN
END SUBROUTINE read_fstoke0
