subroutine ainiaemis

  !  Initialization of adjointised emissions

  use chimere_consts
  use chimere_common
  use wholedomain_common
  
  implicit none
  
  
  allocate(aemisb(0:nhourrun+1,nemisb,nzonal_domain,nmerid_domain))
  allocate(aemisaloc(0:nhourrun+1,nemisa,nzonal_domain,nmerid_domain,nlevemis))

  aemisb(:,:,:,:)=dzero
  aemisaloc(:,:,:,:,:)=dzero

end subroutine ainiaemis
