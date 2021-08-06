subroutine ainiaboun
 
  !  Initialization of adjointised immissions

  use chimere_consts
  use chimere_common
  use wholedomain_common

  implicit none
  
  
  allocate(aboundtop(nzonal_domain,nmerid_domain,nspec,0:nhourrun+1))
  allocate(aboundlat(nlatbound_domain,nspec,0:nhourrun+1))
  !! RQ: dimension temporelle programmee de facon trop brutale
  ! il faut jouer sur 1/2 et faire glisser normalement!!!

  aboundtop(:,:,:,:)=dzero
  aboundlat(:,:,:)=dzero


end subroutine ainiaboun
