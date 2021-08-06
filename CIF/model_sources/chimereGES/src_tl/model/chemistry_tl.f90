subroutine chemistry_tl(ns,izo,ime,ivert,chpr,chlo) 

  !  Calculation of chemical production and loss terms given the species  
  !  (NS) for the current cell of coordinates izo,ime,ivert                                       
  !  INPUT : NS         Species number 
  !          izo,ime,ivert coordinates of the current point
  !          KREACL     Number of reactions consuming species NS          
  !          IREACL     The Adresses of reactions consuming species NS    
  !          RATE       The K or J values for all reactions               
  !          NREACTANTS Number of reactants                               
  !          CONC       Array of instantaneous concentrations             
  !          KREACP     Number of reactions producing species NS          
  !          IREACP     The Adresses of reactions producing species NS    
  !          STOI       Array of Stoichiometric coefficients              
  !          WGTSL      Lower weight for interpol. between tabulated tempe
  !          WGTSL      Upper weight for interpol. between tabulated tempe
  !  OUTPUT: CHPR       Chemical production of species NS in box NB       
  !          CHLO       Chemical loss       of species NS in box NB       

  use chimere_consts
  use worker_common

  implicit none

  include 'mpif.h'


  !*****************************************************************************************
  ! subroutine arguments
  integer      :: ns
  integer      :: izo,ime,ivert
  real(kind=8) :: chpr,chlo
  ! local variables
  integer      :: kr,it,ir,itt
  integer      :: nr,nrr
  real(kind=8) :: trat,prod

  !*****************************************************************************************
  !  Initializations                                                      

  chlo = dzero 
  chpr = dzero 

  !  Loss Terms                                                           

  do kr=1,kreacl(ns) 
    ir = ireacl(ns,kr)
    trat=0d0
    do it=1,nreactants(ir) 
      prod=1d0
      do itt=1,nreactants(ir)
        if(itt.ne.it) prod = prod*conc(irctt(ir,itt),izo,ime,ivert)
      enddo
      trat=trat+prod*conc_tl(irctt(ir,it),izo,ime,ivert)
    enddo
    trat = trat*rate(ir,izo,ime,ivert)
    chlo = chlo + trat 
  enddo

  !  Production Terms                                                     
  do kr=1,kreacp(ns) 
    ir = ireacp(ns,kr) 
    trat = 0d0
    do nr=1,nreactants(ir) 
       prod=1d0
       do nrr=1,nreactants(ir) 
         if(nrr.ne.nr)prod = prod*conc(irctt(ir,nrr),izo,ime,ivert)
       enddo
       trat = trat+prod*conc_tl(irctt(ir,nr),izo,ime,ivert) 
    enddo
    trat=trat*rate(ir,izo,ime,ivert)
    chpr = chpr + trat*(wgstl(izo,ime,ivert)*stoi(ns,ir,istoit(izo,ime,ivert))           &
          &                     + wgsth(izo,ime,ivert)*stoi(ns,ir,istoit(izo,ime,ivert)+1))    
  enddo
  
  return 
END subroutine chemistry_tl
