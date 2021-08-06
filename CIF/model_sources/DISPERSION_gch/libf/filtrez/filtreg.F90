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

SUBROUTINE filtreg ( champ, ibeg, iend, nlat, nbniv,  &
   ifiltre, iarea, griscal ,iter)
  USE Parallel, ONLY : OMP_CHUNK 
  IMPLICIT NONE

  !=======================================================================
  !
  !   Auteur: P. Le Van        07/10/97
  !   ------
  !
  !   Objet: filtre matriciel longitudinal ,avec les matrices precalculees
  !                     pour l'operateur  Filtre    .
  !   ------
  !
  !   Arguments:
  !   ----------
  !
  !      
  !      ibeg..iend            lattitude a filtrer
  !      nlat                  nombre de latitudes du champ
  !      nbniv                 nombre de niveaux verticaux a filtrer
  !      champ(iip1,nblat,nbniv)  en entree : champ a filtrer
  !                            en sortie : champ filtre
  !      ifiltre               +1  Transformee directe
  !                            -1  Transformee inverse
  !                            +2  Filtre directe
  !                            -2  Filtre inverse
  !
  !      iaire                 1   si champ intensif
  !                            2   si champ extensif (pondere par les aires)
  !
  !      iter                  1   filtre simple
  !
  !=======================================================================
  !
  !
  !                      Variable Intensive
  !                ifiltre = 1     filtre directe
  !                ifiltre =-1     filtre inverse
  !
  !                      Variable Extensive
  !                ifiltre = 2     filtre directe
  !                ifiltre =-2     filtre inverse
  !
  !
  INCLUDE "dimensions.h"
  INCLUDE "paramet.h"
  INCLUDE "parafilt.h"
  INCLUDE "coefils.h"
  !
  INTEGER :: ibeg,iend,nlat,nbniv,ifiltre,iter
  INTEGER :: i,j,l,k
  INTEGER :: iim2,immjm
  INTEGER :: jdfil1,jdfil2,jffil1,jffil2,jdfil,jffil
  
  REAL  :: champ( iip1,nlat,nbniv)
  REAL :: matriceun,matriceus,matricevn,matricevs,matrinvn,matrinvs
  COMMON/matrfil/matriceun(iim,iim,nfilun),matriceus(iim,iim,nfilus) &
     , matricevn(iim,iim,nfilvn),matricevs(iim,iim,nfilvs) &
     ,  matrinvn(iim,iim,nfilun),matrinvs (iim,iim,nfilus)
  REAL :: eignq(iim), sdd1(iim),sdd2(iim)
  LOGICAL ::   griscal
  INTEGER ::   hemisph, iarea
  !
  
  IF(ifiltre==1.OR.ifiltre==-1) & 
     STOP'Pas de transformee simple dans cette version'
  
  IF( iter== 2 )  THEN
      PRINT *,' Pas d iteration du filtre dans cette version !' &
         , ' Utiliser old_filtreg et repasser !'
      STOP
  ENDIF
  
  IF( ifiltre== -2 .AND..NOT.griscal )     THEN
      PRINT *,' Cette routine ne calcule le filtre inverse que ', &
         ' sur la grille des scalaires !'
      STOP
  ENDIF

  IF( ifiltre/=2 .AND.ifiltre/= - 2 )  THEN
      PRINT *,' Probleme dans filtreg car ifiltre NE 2 et NE -2' &
         ,' corriger et repasser !'
      STOP
  ENDIF
  !
  
  iim2   = iim * iim
  immjm  = iim * jjm
  !
  !
  IF( griscal )   THEN
      IF( nlat/= jjp1 )  THEN
          PRINT  1111
          STOP
      ELSE
          !
          IF( iarea==1 )  THEN
              CALL SCOPY(  iim,    sddv, 1,  sdd1, 1 ) 
              CALL SCOPY(  iim,  unsddv, 1,  sdd2, 1 )
          ELSE
              CALL SCOPY(  iim,  unsddv, 1,  sdd1, 1 )
              CALL SCOPY(  iim,    sddv, 1,  sdd2, 1 )
          END IF
          !
          jdfil1 = 2
          jffil1 = jfiltnu
          jdfil2 = jfiltsu
          jffil2 = jjm
      END IF
  ELSE
      IF( nlat/=jjm )  THEN
          PRINT  2222
          STOP
      ELSE
          !
          IF( iarea==1 )  THEN
              CALL SCOPY(  iim,    sddu, 1,  sdd1, 1 ) 
              CALL SCOPY(  iim,  unsddu, 1,  sdd2, 1 )
          ELSE
              CALL SCOPY(  iim,  unsddu, 1,  sdd1, 1 )
              CALL SCOPY(  iim,    sddu, 1,  sdd2, 1 )
          END IF
          !
          jdfil1 = 1
          jffil1 = jfiltnv
          jdfil2 = jfiltsv
          jffil2 = jjm
      END IF
  END IF
  !
  !
  DO hemisph = 1, 2
    !
    IF ( hemisph==1 )  THEN
        jdfil = MAX(jdfil1,ibeg)
        jffil = MIN(jffil1,iend)
    ELSE
        jdfil = MAX(jdfil2,ibeg)
        jffil = MIN(jffil2,iend)
    END IF
    
    !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
    DO l = 1, nbniv
      DO j = jdfil,jffil
        
        
        DO  i = 1, iim
          champ(i,j,l) = champ(i,j,l) * sdd1(i)
        ENDDO
        !
        
        IF( hemisph== 1 )      THEN
            
            IF( ifiltre== -2 )   THEN
#ifdef CRAY
                CALL MXVA( matrinvn(1,1,j), 1, iim, champ(1,j,l), 1, eignq  ,  &
                   1, iim, iim                         )
#else
#ifdef BLAS
                CALL SGEMV("N", iim,iim, 1.0, matrinvn(1,1,j),iim, &
                   champ(1,j,l), 1, 0.0, eignq, 1)
#else
                DO k = 1, iim
                  eignq(k) = 0.0
                ENDDO
                DO k = 1, iim
                  DO i = 1, iim
                    eignq(k) = eignq(k) + matrinvn(k,i,j)*champ(i,j,l)
                  ENDDO
                ENDDO
#endif
#endif
            ELSE IF ( griscal )     THEN
#ifdef CRAY
                CALL MXVA( matriceun(1,1,j), 1, iim, champ(1,j,l), 1, eignq , &
                   1, iim, iim                         )
#else
#ifdef BLAS
                CALL SGEMV("N", iim,iim, 1.0, matriceun(1,1,j),iim, &
                   champ(1,j,l), 1, 0.0, eignq, 1)
#else
                DO k = 1, iim
                  eignq(k) = 0.0
                ENDDO
                DO i = 1, iim
                  DO k = 1, iim
                    eignq(k) = eignq(k) + matriceun(k,i,j)*champ(i,j,l)
                  ENDDO
                ENDDO
#endif
#endif
            ELSE 
#ifdef CRAY
                CALL MXVA( matricevn(1,1,j), 1, iim, champ(1,j,l), 1, eignq ,  &
                   1, iim, iim                         )
#else
#ifdef BLAS
                CALL SGEMV("N", iim,iim, 1.0, matricevn(1,1,j),iim, &
                   champ(1,j,l), 1, 0.0, eignq, 1)
#else
                DO k = 1, iim
                  eignq(k) = 0.0
                ENDDO
                DO i = 1, iim
                  DO k = 1, iim
                    eignq(k) = eignq(k) + matricevn(k,i,j)*champ(i,j,l)
                  ENDDO
                ENDDO
#endif
#endif
            ENDIF

        ELSE

            IF( ifiltre== -2 )   THEN
#ifdef CRAY
                CALL MXVA( matrinvs(1,1,j-jfiltsu+1),  1, iim, champ(1,j,l),1 ,   &
                   eignq,  1, iim, iim                    )
#else
#ifdef BLAS
                CALL SGEMV("N", iim,iim, 1.0, matrinvs(1,1,j-jfiltsu+1),iim, &
                   champ(1,j,l), 1, 0.0, eignq, 1)
#else
                DO k = 1, iim
                  eignq(k) = 0.0
                ENDDO
                DO i = 1, iim
                  DO k = 1, iim
                    eignq(k) = eignq(k) + matrinvs(k,i,j-jfiltsu+1)*champ(i,j,l)
                  ENDDO
                ENDDO
#endif
#endif
            ELSE IF ( griscal )     THEN
#ifdef CRAY
                CALL MXVA( matriceus(1,1,j-jfiltsu+1), 1, iim, champ(1,j,l),1 ,  &
                   eignq,  1, iim, iim                    )
#else
#ifdef BLAS
                CALL SGEMV("N", iim,iim, 1.0, matriceus(1,1,j-jfiltsu+1),iim, &
                   champ(1,j,l), 1, 0.0, eignq, 1)
#else
                DO k = 1, iim
                  eignq(k) = 0.0
                ENDDO
                DO i = 1, iim
                  DO k = 1, iim
                    eignq(k) = eignq(k) + matriceus(k,i,j-jfiltsu+1)*champ(i,j,l)
                  ENDDO
                ENDDO
#endif
#endif
            ELSE 
#ifdef CRAY
                CALL MXVA( matricevs(1,1,j-jfiltsv+1), 1, iim, champ(1,j,l),1 ,  &
                   eignq,  1, iim, iim                    )
#else
#ifdef BLAS
                CALL SGEMV("N", iim,iim, 1.0, matricevs(1,1,j-jfiltsv+1),iim, &
                   champ(1,j,l), 1, 0.0, eignq, 1)
#else
                DO k = 1, iim
                  eignq(k) = 0.0
                ENDDO
                DO i = 1, iim
                  DO k = 1, iim
                    eignq(k) = eignq(k) + matricevs(k,i,j-jfiltsv+1)*champ(i,j,l)
                  ENDDO
                ENDDO
#endif
#endif
            ENDIF
            
        ENDIF
        !
        IF( ifiltre== 2 )  THEN
            DO i = 1, iim
              champ( i,j,l ) = ( champ(i,j,l) + eignq(i) ) * sdd2(i)
            ENDDO
        ELSE
            DO i=1,iim
              champ( i,j,l ) = ( champ(i,j,l) - eignq(i) ) * sdd2(i)
            ENDDO
        ENDIF
        !
        champ( iip1,j,l ) = champ( 1,j,l )
        !
      ENDDO
      !
    ENDDO
    !$OMP END DO NOWAIT
    !    
  ENDDO
  !
1111 FORMAT(//20x,'ERREUR dans le dimensionnement du tableau  CHAMP a &
        &filtrer, sur la grille des scalaires'/)
2222 FORMAT(//20x,'ERREUR dans le dimensionnement du tableau CHAMP a fi &
        &ltrer, sur la grille de V ou de Z'/)
  RETURN
END SUBROUTINE filtreg
