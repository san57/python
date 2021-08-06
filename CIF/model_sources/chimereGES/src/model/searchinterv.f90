subroutine searchinterv(val,vect,ndim,ind,sav)

  implicit none

  ! local variables
  integer :: ndim,ind,indic,i
  real(kind=8) :: val,sav
  real(kind=8),dimension(ndim) :: vect

  if(val.le.vect(1)) then
     ind=1
     sav=vect(1)
     goto 1001
  endif
  if(val.ge.vect(ndim)) then
     ind=ndim-1
     sav=vect(ndim)
     goto 1001
  endif
  !
  do i=ind,ndim-1
     if(val.gt.vect(i).and.val.le.vect(i+1)) then
        ind=i
        sav=val
        goto 1001
     endif
  enddo
  do i=ind-1,1,-1
     if(val.gt.vect(i).and.val.le.vect(i+1)) then
        ind=i
        sav=val
        goto 1001
     endif
  enddo
1001 continue
  !
  return


end subroutine searchinterv
