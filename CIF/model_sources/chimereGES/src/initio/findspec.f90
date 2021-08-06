subroutine findspec(charspec,nospec)

  !  Returns the index of the species, if found among active ones,
  !  and 0 if the species is not found.
  !  INPUT  : CHARSPEC         The species name
  !           SPECIES          Active species names array
  !  OUTPUT : NOSPEC           The species index

  use chimere_consts
  use chimere_common
  implicit none

  !*****************************************************************************************
  ! subroutine arguments
  character(len=*),intent(in) :: charspec
  integer,intent(out)             :: nospec
  ! local variables
  integer :: ns
  !*****************************************************************************************

  nospec = 0
  do ns=1,nspectot
     if(charspec.eq.species(ns)%name) then
        nospec = ns
        exit
     endif
  enddo

end subroutine findspec
