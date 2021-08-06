SUBROUTINE thermcell_dq(ngrid,nlay,ptimestep,fm,entr,masse,q,dq)


 USE dimphy
 IMPLICIT NONE
 INCLUDE "dimensions.h"

!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
! Modif 2013/01/04 (FH hourdin@lmd.jussieu.fr)
!  Introduction of an implicit computation of vertical advection in
!  the environment of thermal plumes in thermcell_dq
!  Schéma implicite
!
!=======================================================================

INTEGER :: ngrid,nlay

REAL :: ptimestep
REAL :: masse(ngrid,nlay),fm(ngrid,nlay+1)
REAL :: entr(ngrid,nlay)
REAL :: q(ngrid,nlay),qtemp(ngrid,nlay)
REAL :: dq(ngrid,nlay)
REAL :: qa(ngrid,nlay),detr(ngrid,nlay)
REAL :: zzm
INTEGER :: ig,k
REAL :: cfl
REAL :: qold(ngrid,nlay),fqa(ngrid,nlay+1)
INTEGER :: niter,iter


! Calcul du critere CFL pour l'advection dans la subsidence
! cfl = 0.
! do k=1,nlay
!  do ig=1,ngrid
!   zzm=masse(ig,k)/ptimestep
!   cfl=max(cfl,fm(ig,k)/zzm)
!   if (entr(ig,k).gt.zzm) then
!    print*,'entr dt > m ',entr(ig,k)*ptimestep,masse(ig,k)
!   endif
!  enddo
! enddo

!Stockage de q dans qold
 qold=q

! Calcul du detrainement
 do k=1,nlay
  do ig=1,ngrid
   detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
   if (detr(ig,k).lt.0.) then
    entr(ig,k)=entr(ig,k)-detr(ig,k)
    detr(ig,k)=0.
   endif
  enddo
 enddo

! Computation of tracer concentrations in the ascending plume
 do ig=1,ngrid
  qa(ig,1)=q(ig,1)
 enddo

 do k=2,nlay
  do ig=1,ngrid
   if ((fm(ig,k+1)+detr(ig,k))*ptimestep .gt. 1.e-5*masse(ig,k)) then
    qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))/(fm(ig,k+1)+detr(ig,k))
   else
    qa(ig,k)=q(ig,k)
   endif
  enddo
 enddo

! Plume vertical flux
 do k=2,nlay-1
  fqa(:,k)=fm(:,k)*qa(:,k-1)
 enddo
 fqa(:,1)=0. 
 fqa(:,nlay)=0.


! Trace species evolution
 do k=nlay-1,1,-1
   q(:,k)=(q(:,k)+ptimestep/masse(:,k)*(fqa(:,k)-fqa(:,k+1)+fm(:,k+1)*q(:,k+1)))/(1.+fm(:,k)*ptimestep/masse(:,k))
 enddo

! Tendencies
 do k=1,nlay
  do ig=1,ngrid
   dq(ig,k)=(q(ig,k)-qold(ig,k))/ptimestep
   q(ig,k)=qold(ig,k)
  enddo
 enddo

RETURN
END SUBROUTINE
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
SUBROUTINE thermcell_dq_tl(ngrid,nlay,ptimestep,fm,entr,masse,q,q_tl,dq,dq_tl)


 USE dimphy
 IMPLICIT NONE
 INCLUDE "dimensions.h"

!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
! Modif 2013/01/04 (FH hourdin@lmd.jussieu.fr)
!  Introduction of an implicit computation of vertical advection in
!  the environment of thermal plumes in thermcell_dq
!  Schéma implicite
!
!=======================================================================

INTEGER :: ngrid,nlay

REAL :: ptimestep
REAL :: masse(ngrid,nlay),fm(ngrid,nlay+1)
REAL :: entr(ngrid,nlay)
REAL :: q(ngrid,nlay),q_tl(ngrid,nlay)
REAL :: qtemp(ngrid,nlay),qtemp_tl(ngrid,nlay)
REAL :: dq(ngrid,nlay),dq_tl(ngrid,nlay)
REAL :: qa(ngrid,nlay),qa_tl(ngrid,nlay)
REAL :: detr(ngrid,nlay)
INTEGER :: ig,k
REAL :: qold(ngrid,nlay),qold_tl(ngrid,nlay)
REAL :: fqa(ngrid,nlay+1),fqa_tl(ngrid,nlay+1)
INTEGER :: niter,iter


!Stockage de q dans qold
 qold_tl=q_tl
 qold=q


! Calcul du detrainement
 do k=1,nlay
  do ig=1,ngrid
   detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
   if (detr(ig,k).lt.0.) then
    entr(ig,k)=entr(ig,k)-detr(ig,k)
    detr(ig,k)=0.
   endif
  enddo
 enddo

! Computation of tracer concentrations in the ascending plume
 do ig=1,ngrid
  qa_tl(ig,1)=q_tl(ig,1)
  qa(ig,1)=q(ig,1)
 enddo

 do k=2,nlay
  do ig=1,ngrid
   if ((fm(ig,k+1)+detr(ig,k))*ptimestep .gt. 1.e-5*masse(ig,k)) then
    qa_tl(ig,k)=(fm(ig,k)*qa_tl(ig,k-1)+entr(ig,k)*q_tl(ig,k))/(fm(ig,k+1)+detr(ig,k))
    qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))/(fm(ig,k+1)+detr(ig,k))
   else
    qa_tl(ig,k)=q_tl(ig,k)
    qa(ig,k)=q(ig,k)
   endif
  enddo
 enddo

! Plume vertical flux
 do k=2,nlay-1
  fqa_tl(:,k)=fm(:,k)*qa_tl(:,k-1)
  fqa(:,k)=fm(:,k)*qa(:,k-1)
 enddo
 fqa_tl(:,1)=0. 
 fqa_tl(:,nlay)=0. 
 fqa(:,1)=0. 
 fqa(:,nlay)=0.


! Trace species evolution
 do k=nlay-1,1,-1
  q(:,k)=(q(:,k)+ptimestep/masse(:,k)*(fqa(:,k)-fqa(:,k+1)+fm(:,k+1)*q(:,k+1)))/(1.+fm(:,k)*ptimestep/masse(:,k))
  q_tl(:,k)=(q_tl(:,k)+ptimestep/masse(:,k)*(fqa_tl(:,k)-fqa_tl(:,k+1)+fm(:,k+1)*q_tl(:,k+1)))/(1.+fm(:,k)*ptimestep/masse(:,k))
 enddo



! Tendencies
 do k=1,nlay
  do ig=1,ngrid
   dq_tl(ig,k)=(q_tl(ig,k)-qold_tl(ig,k))/ptimestep
   q_tl(ig,k)=qold_tl(ig,k)
   dq(ig,k)=(q(ig,k)-qold(ig,k))/ptimestep
   q(ig,k)=qold(ig,k)
  enddo
 enddo

RETURN
END SUBROUTINE
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
SUBROUTINE thermcell_dq_ad(ngrid,nlay,ptimestep,fm,entr,masse,q,q_ad,dq,dq_ad)


 USE dimphy
 IMPLICIT NONE
 INCLUDE "dimensions.h"

!=======================================================================
!
!   Calcul du transport verticale dans la couche limite en presence
!   de "thermiques" explicitement representes
!   calcul du dq/dt une fois qu'on connait les ascendances
!
! Modif 2013/01/04 (FH hourdin@lmd.jussieu.fr)
!  Introduction of an implicit computation of vertical advection in
!  the environment of thermal plumes in thermcell_dq
!  Schéma implicite
!
!=======================================================================

INTEGER :: ngrid,nlay

REAL :: ptimestep
REAL :: masse(ngrid,nlay),fm(ngrid,nlay+1)
REAL :: entr(ngrid,nlay)
REAL :: q(ngrid,nlay),q_ad(ngrid,nlay)
REAL :: qtemp(ngrid,nlay),qtemp_ad(ngrid,nlay)
REAL :: dq(ngrid,nlay),dq_ad(ngrid,nlay)
REAL :: qa(ngrid,nlay),qa_ad(ngrid,nlay)
REAL :: detr(ngrid,nlay)
INTEGER :: ig,k
REAL :: qold(ngrid,nlay),qold_ad(ngrid,nlay)
REAL :: fqa(ngrid,nlay+1),fqa_ad(ngrid,nlay+1)
INTEGER :: niter,iter

!Initialisation variables locales + q_ad
qold_ad(:,:)=0.
fqa_ad(:,:)=0.
qa_ad(:,:)=0.
qtemp_ad(:,:)=0.


!!! !Debut du direct
!!! qold=q

! Calcul du detrainement
 do k=1,nlay
  do ig=1,ngrid
   detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
   if (detr(ig,k).lt.0.) then
    entr(ig,k)=entr(ig,k)-detr(ig,k)
    detr(ig,k)=0.
   endif
  enddo
 enddo

! Tendencies
 do k=1,nlay
  do ig=1,ngrid

   !Start of Adjoint (ICI)
   qold_ad(ig,k) = q_ad(ig,k)
   q_ad(ig,k) = 0.
   q_ad(ig,k) = dq_ad(ig,k)/ptimestep
   qold_ad(ig,k) = qold_ad(ig,k) - dq_ad(ig,k)/ptimestep
   dq_ad(ig,k)=0.
  enddo
 enddo

! Trace species evolution
 do k=1,nlay-1
  fqa_ad(:,k)=fqa_ad(:,k)+q_ad(:,k)*ptimestep/masse(:,k)/(1.+fm(:,k)*ptimestep/masse(:,k))
  fqa_ad(:,k+1)=fqa_ad(:,k+1)-q_ad(:,k)*ptimestep/masse(:,k)/(1.+fm(:,k)*ptimestep/masse(:,k))
  q_ad(:,k+1)=q_ad(:,k+1)+fm(:,k+1)*q_ad(:,k)*ptimestep/masse(:,k)/(1.+fm(:,k)*ptimestep/masse(:,k))
  q_ad(:,k)=q_ad(:,k)/(1.+fm(:,k)*ptimestep/masse(:,k))
 enddo

! Plume vertical flux
 fqa_ad(:,1)=0. 
 fqa_ad(:,nlay)=0.
 do k=nlay-1,2,-1
  qa_ad(:,k-1)=qa_ad(:,k-1)+fm(:,k)*fqa_ad(:,k)
  fqa_ad(:,k)=0.
 enddo

! Computation of tracer concentrations in the ascending plume 
 do k=nlay,2,-1
  do ig=1,ngrid
   if ((fm(ig,k+1)+detr(ig,k))*ptimestep .gt. 1.e-5*masse(ig,k)) then
    q_ad(ig,k)=q_ad(ig,k)+entr(ig,k)*qa_ad(ig,k)/(fm(ig,k+1)+detr(ig,k))
    qa_ad(ig,k-1)=qa_ad(ig,k-1)+fm(ig,k)*qa_ad(ig,k)/(fm(ig,k+1)+detr(ig,k))
    qa_ad(ig,k)=0.
   else
    q_ad(ig,k)=q_ad(ig,k)+qa_ad(ig,k)
    qa_ad(ig,k)=0.
   endif
  enddo
 enddo 

 do ig=1,ngrid
  q_ad(ig,1)=qa_ad(ig,1)+q_ad(ig,1)
  qa_ad(ig,1)=0.
 enddo

 q_ad=qold_ad+q_ad
 qold_ad=0.


RETURN
END SUBROUTINE

