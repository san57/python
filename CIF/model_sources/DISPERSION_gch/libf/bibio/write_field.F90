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

MODULE write_field
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: MaxWriteField = 100
  INTEGER, DIMENSION(MaxWriteField),SAVE :: FieldId
  INTEGER, DIMENSION(MaxWriteField),SAVE :: FieldVarId
  INTEGER, DIMENSION(MaxWriteField),SAVE :: FieldIndex
  CHARACTER(len=255), DIMENSION(MaxWriteField) ::  FieldName 
  
  INTEGER,SAVE :: NbField = 0
  
  INTERFACE WriteField
    MODULE PROCEDURE WriteField3d,WriteField2d,WriteField1d
  END INTERFACE
CONTAINS
  
  FUNCTION GetFieldIndex(name)
    IMPLICIT NONE
    INTEGER          :: GetFieldindex
    CHARACTER(len=*) :: name
    
    CHARACTER(len=255) :: TrueName
    INTEGER            :: i
    
    
    TrueName=TRIM(ADJUSTL(name))
    
    GetFieldIndex=-1
    DO i=1,NbField
      IF (TrueName==FieldName(i)) THEN
          GetFieldIndex=i
          EXIT
      ENDIF
    ENDDO
  END FUNCTION GetFieldIndex
  
  SUBROUTINE WriteField3d(name,Field)
    IMPLICIT NONE
    CHARACTER(len=*) :: name
    REAL, DIMENSION(:,:,:) :: Field 
    INTEGER, DIMENSION(3) :: Dim
    
    Dim=SHAPE(Field)
    CALL WriteField_gen(name,Field,DIM(1),DIM(2),DIM(3))  
    
  END SUBROUTINE WriteField3d
  
  SUBROUTINE WriteField2d(name,Field)
    IMPLICIT NONE
    CHARACTER(len=*) :: name
    REAL, DIMENSION(:,:) :: Field 
    INTEGER, DIMENSION(2) :: Dim
    
    Dim=SHAPE(Field)
    CALL WriteField_gen(name,Field,DIM(1),DIM(2),1)  
    
  END SUBROUTINE WriteField2d
  
  SUBROUTINE WriteField1d(name,Field)
    IMPLICIT NONE
    CHARACTER(len=*) :: name
    REAL, DIMENSION(:) :: Field 
    INTEGER, DIMENSION(1) :: Dim
    
    Dim=SHAPE(Field)
    CALL WriteField_gen(name,Field,DIM(1),1,1)  
    
  END SUBROUTINE WriteField1d
  
  SUBROUTINE WriteField_gen(name,Field,dimx,dimy,dimz)
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    CHARACTER(len=*) :: name
    INTEGER :: dimx,dimy,dimz
    REAL,DIMENSION(dimx,dimy,dimz) :: Field
    INTEGER :: status
    INTEGER :: index
    INTEGER :: start(4)
    INTEGER :: COUNT(4)
    
    
    Index=GetFieldIndex(name)
    IF (Index==-1) THEN
        CALL CreateNewField(name,dimx,dimy,dimz)
	Index=GetFieldIndex(name)
    ELSE
        FieldIndex(Index)=FieldIndex(Index)+1.
    ENDIF
    
    start(1)=1
    start(2)=1
    start(3)=1
    start(4)=FieldIndex(Index)
    
    COUNT(1)=dimx
    COUNT(2)=dimy
    COUNT(3)=dimz
    COUNT(4)=1
    
    status = NF_PUT_VARA_DOUBLE(FieldId(Index),FieldVarId(Index),start,count,Field)
    status = NF_SYNC(FieldId(Index))
    
  END SUBROUTINE WriteField_gen
  
  SUBROUTINE CreateNewField(name,dimx,dimy,dimz)
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'  
    CHARACTER(len=*) :: name
    INTEGER :: dimx,dimy,dimz
    INTEGER :: TabDim(4)
    INTEGER :: status
    
    
    NbField=NbField+1
    FieldName(NbField)=TRIM(ADJUSTL(name))
    FieldIndex(NbField)=1
    
    
    status = NF_CREATE(TRIM(ADJUSTL(name))//'.nc', NF_CLOBBER, FieldId(NbField))
    status = NF_DEF_DIM(FieldId(NbField),'X',dimx,TabDim(1))
    status = NF_DEF_DIM(FieldId(NbField),'Y',dimy,TabDim(2))
    status = NF_DEF_DIM(FieldId(NbField),'Z',dimz,TabDim(3))
    status = NF_DEF_DIM(FieldId(NbField),'iter',NF_UNLIMITED,TabDim(4))
    status = NF_DEF_VAR(FieldId(NbField),FieldName(NbField),NF_DOUBLE,4,TabDim,FieldVarId(NbField))
    status = NF_ENDDEF(FieldId(NbField))
    
  END SUBROUTINE CreateNewField
  
  
  
  SUBROUTINE write_field1D(name,Field)
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: MaxDim=1
    CHARACTER(len=*)   :: name
    REAL, DIMENSION(:) :: Field
    INTEGER, DIMENSION(MaxDim) :: Dim
    INTEGER :: i,nb
    INTEGER, PARAMETER :: id=10
    INTEGER, PARAMETER :: NbCol=4
    INTEGER :: ColumnSize 
    INTEGER :: pos
    CHARACTER(len=255) :: form
    CHARACTER(len=255) :: MaxLen
    
    
    OPEN(unit=id,file=name//'.field',form='formatted',status='replace')
    WRITE (id,'("----- Field '//name//'",//)')
    Dim=SHAPE(Field)
    MaxLen=int2str(LEN(TRIM(int2str(DIM(1)))))
    ColumnSize=20+6+3+LEN(TRIM(int2str(DIM(1))))
    Nb=0
    Pos=2
    DO i=1,DIM(1)
      nb=nb+1
      
      IF (MOD(nb,NbCol)==0) THEN
          form='(t'//TRIM(int2str(pos))// ',i'//TRIM(MaxLen) //'," ---> ",g22.16,/)'
          Pos=2
      ELSE
          form='(t'//TRIM(int2str(pos))// ',i'//TRIM(MaxLen) //'," ---> ",g22.16," | ",)'
          Pos=Pos+ColumnSize
      ENDIF
      WRITE (id,form,advance='no') i,Field(i)
    ENDDO
    
    CLOSE(id)
    
  END SUBROUTINE write_field1D
  
  SUBROUTINE write_field2D(name,Field)
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: MaxDim=2
    CHARACTER(len=*)   :: name
    REAL, DIMENSION(:,:) :: Field
    INTEGER, DIMENSION(MaxDim) :: Dim
    INTEGER :: i,j,nb
    INTEGER, PARAMETER :: id=10
    INTEGER, PARAMETER :: NbCol=4
    INTEGER :: ColumnSize 
    INTEGER :: pos,offset
    CHARACTER(len=255) :: form
    CHARACTER(len=255) :: spacing
    
    OPEN(unit=id,file=name//'.field',form='formatted',status='replace')
    WRITE (id,'("----- Field '//name//'",//)')
    
    Dim=SHAPE(Field)
    offset=LEN(TRIM(int2str(DIM(1))))+LEN(TRIM(int2str(DIM(2))))+3
    ColumnSize=20+6+3+offset
    
    spacing='(t2,"'//REPEAT('-',ColumnSize*NbCol)//'")'
    
    DO i=1,DIM(2)
      nb=0
      Pos=2
      DO j=1,DIM(1)
        nb=nb+1
        
        IF (MOD(nb,NbCol)==0) THEN
            form='(t'//TRIM(int2str(pos))//            &
               ',"('//TRIM(int2str(j))//','          &
               //TRIM(int2str(i))//')",t'       & 
               //TRIM(int2str(pos+offset))     &    
               //'," ---> ",g22.16,/)'
            Pos=2
        ELSE
            form='(t'//TRIM(int2str(pos))//            &
               ',"('//TRIM(int2str(j))//','          &
               //TRIM(int2str(i))//')",t'       & 
               //TRIM(int2str(pos+offset))     &    
               //'," ---> ",g22.16," | ")'
            Pos=Pos+ColumnSize
        ENDIF
        WRITE (id,form,advance='no') Field(j,i)
      ENDDO
      IF (MOD(nb,NbCol)==0) THEN
          WRITE (id,spacing)
      ELSE
          WRITE (id,'')
          WRITE (id,spacing)
      ENDIF
    ENDDO
    
  END SUBROUTINE write_field2D
  
  SUBROUTINE write_field3D(name,Field)
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: MaxDim=3
    CHARACTER(len=*)   :: name
    REAL, DIMENSION(:,:,:) :: Field
    INTEGER, DIMENSION(MaxDim) :: Dim
    INTEGER :: i,j,k,nb
    INTEGER, PARAMETER :: id=10
    INTEGER, PARAMETER :: NbCol=4
    INTEGER :: ColumnSize 
    INTEGER :: pos,offset
    CHARACTER(len=255) :: form
    CHARACTER(len=255) :: spacing
    
    OPEN(unit=id,file=name//'.field',form='formatted',status='replace')
    WRITE (id,'("----- Field '//name//'"//)')
    
    Dim=SHAPE(Field)
    offset=LEN(TRIM(int2str(DIM(1))))+LEN(TRIM(int2str(DIM(2))))+LEN(TRIM(int2str(DIM(3))))+4
    ColumnSize=22+6+3+offset
    
    !    open(unit=id,file=name,form=formatted
    
    spacing='(t2,"'//REPEAT('-',ColumnSize*NbCol)//'")'
    
    DO i=1,DIM(3)
    
      DO j=1,DIM(2)
        nb=0
        Pos=2
        
        DO k=1,DIM(1)
          nb=nb+1
          
          IF (MOD(nb,NbCol)==0) THEN
              form='(t'//TRIM(int2str(pos))//            &
                 ',"('//TRIM(int2str(k))//','          &
                 //TRIM(int2str(j))//','          &
                 //TRIM(int2str(i))//')",t'       & 
                 //TRIM(int2str(pos+offset))      &    
                 //'," ---> ",g22.16,/)'
              Pos=2
          ELSE
              form='(t'//TRIM(int2str(pos))//            &
                 ',"('//TRIM(int2str(k))//','          &
                 //TRIM(int2str(j))//','          &
                 //TRIM(int2str(i))//')",t'       & 
                 //TRIM(int2str(pos+offset))      &    
                 //'," ---> ",g22.16," | ")'
              ! d�pent de l'impl�mention, sur compaq, c'est necessarea
              !            Pos=Pos+ColumnSize
          ENDIF
          WRITE (id,form,advance='no') Field(k,j,i)
        ENDDO
        IF (MOD(nb,NbCol)==0) THEN
            WRITE (id,spacing)
        ELSE
            WRITE (id,'("")')
            WRITE (id,spacing)
        ENDIF
      ENDDO
      WRITE (id,spacing)
    ENDDO
    
    CLOSE(id)
    
  END SUBROUTINE write_field3D
  
  FUNCTION int2str(int)
    IMPLICIT NONE
    INTEGER, PARAMETER :: MaxLen=10
    INTEGER,INTENT(in) :: int
    CHARACTER(len=MaxLen) :: int2str
    LOGICAL :: flag
    INTEGER :: i
    flag=.TRUE.
    
    i=int
    
    int2str=''
    DO WHILE (flag)
      int2str=CHAR(MOD(i,10)+48)//int2str
      i=i/10
      IF (i==0) flag=.FALSE.
    ENDDO
  END FUNCTION int2str
  
END MODULE write_field
  
