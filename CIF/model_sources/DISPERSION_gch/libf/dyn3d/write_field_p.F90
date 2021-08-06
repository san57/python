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

module write_field_p
implicit none
  
  interface WriteField_p
    module procedure Write_field3d_p,Write_Field2d_p,Write_Field1d_p
  end interface WriteField_p
  
  contains
  
  subroutine write_field1D_p(name,Field)
    USE parallel
    USE write_field
    implicit none
  
    integer, parameter :: MaxDim=1
    character(len=*)   :: name
    real, dimension(:) :: Field
    real, dimension(:),allocatable :: New_Field
    integer, dimension(MaxDim) :: Dim
    
    
    Dim=shape(Field)
    allocate(New_Field(Dim(1)))
    New_Field(:)=Field(:)
    call Gather_Field(New_Field,dim(1),1,0)
    
    if (MPI_Rank==0) call WriteField(name,New_Field)
    
    end subroutine write_field1D_p

  subroutine write_field2D_p(name,Field)
    USE parallel
    USE write_field
    implicit none
  
    integer, parameter :: MaxDim=2
    character(len=*)   :: name
    real, dimension(:,:) :: Field
    real, dimension(:,:),allocatable :: New_Field
    integer, dimension(MaxDim) :: Dim
    
    Dim=shape(Field)
    allocate(New_Field(Dim(1),Dim(2)))
    New_Field(:,:)=Field(:,:)
    call Gather_Field(New_Field(1,1),dim(1)*dim(2),1,0)
    
    if (MPI_Rank==0) call WriteField(name,New_Field)
    
     
  end subroutine write_field2D_p
  
  subroutine write_field3D_p(name,Field)
    USE parallel
    USE write_field
    implicit none
  
    integer, parameter :: MaxDim=3
    character(len=*)   :: name
    real, dimension(:,:,:) :: Field
    real, dimension(:,:,:),allocatable :: New_Field
    integer, dimension(MaxDim) :: Dim
    
    Dim=shape(Field)
    allocate(New_Field(Dim(1),Dim(2),Dim(3)))
    New_Field(:,:,:)=Field(:,:,:)
    call Gather_Field(New_Field(1,1,1),dim(1)*dim(2),dim(3),0)
    
   if (MPI_Rank==0) call WriteField(name,New_Field)
    
  end subroutine write_field3D_p  

end module write_field_p
  
