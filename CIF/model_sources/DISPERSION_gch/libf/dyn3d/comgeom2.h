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

!CDK comgeom2
      COMMON/comgeom/ &
      cu(iip1,jjp1),cv(iip1,jjm),unscu2(iip1,jjp1),unscv2(iip1,jjm)  ,&
      area(iip1,jjp1),areasurg(iip1,jjp1),areau(iip1,jjp1)           ,&
      areav(iip1,jjm),unsarea(iip1,jjp1),apoln,apols                 , &
      unsareaz(iip1,jjm),airuscv2(iip1,jjm),airvscu2(iip1,jjm)       ,&
      areaij1(iip1,jjp1),areaij2(iip1,jjp1),areaij3(iip1,jjp1)       ,&
      areaij4(iip1,jjp1),alpha1(iip1,jjp1),alpha2(iip1,jjp1)         ,&
      alpha3(iip1,jjp1),alpha4(iip1,jjp1),alpha1p2(iip1,jjp1)        ,&
      alpha1p4(iip1,jjp1),alpha2p3(iip1,jjp1),alpha3p4(iip1,jjp1)    ,&
      fext(iip1,jjm),constang(iip1,jjp1), rlatu(jjp1),rlatv(jjm),&
      rlonu(iip1),rlonv(iip1),cuvsurcv(iip1,jjm),cvsurcuv(iip1,jjm)  ,&
      cvusurcu(iip1,jjp1),cusurcvu(iip1,jjp1)                        ,&
      cuvscvgam1(iip1,jjm),cuvscvgam2(iip1,jjm),cvuscugam1(iip1,jjp1),&
      cvuscugam2(iip1,jjp1),cvscuvgam(iip1,jjm),cuscvugam(iip1,jjp1) ,&
      unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2                ,&
      unsair_gam1(iip1,jjp1),unsair_gam2(iip1,jjp1)                  ,&
      unsairz_gam(iip1,jjm),aivscu2gam(iip1,jjm),aiuscv2gam(iip1,jjm) &
      , xprimu(iip1),xprimv(iip1)

!
      REAL :: &
      cu,cv,unscu2,unscv2,area,areasurg,areau,areav,apoln,apols,unsarea &
      ,unsareaz,airuscv2,airvscu2,areaij1,areaij2,areaij3,areaij4     ,&
      alpha1,alpha2,alpha3,alpha4,alpha1p2,alpha1p4,alpha2p3,alpha3p4 ,&
      fext,constang,rlatu,rlatv,rlonu,rlonv,cuvscvgam1,cuvscvgam2     ,&
      cvuscugam1,cvuscugam2,cvscuvgam,cuscvugam,unsapolnga1           , &
      unsapolnga2,unsapolsga1,unsapolsga2,unsair_gam1,unsair_gam2     ,&
      unsairz_gam,aivscu2gam,aiuscv2gam,cuvsurcv,cvsurcuv,cvusurcu    ,&
      cusurcvu,xprimu,xprimv
!
