subroutine inichem

  !  Initialization of chemical arrays

  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

  !*****************************************************************************************
  real(kind=8),dimension(ntabuzenmax)   :: zenang
  real,dimension(ntemps)                :: fr
  real(kind=8)                          :: xmass
  real(kind=8)                          :: a1,a2,a3
  real              :: dum1,dum2,dum3
  character(len=15) :: dummy

  integer :: ns,nf,ne,nsp,nread,lt,ity,noreac,noprod
  integer :: i,k,l,ii
  integer :: np,nr,nre,it,nt,nl,nw
  integer :: nprods
  integer :: iread
  integer :: iloss,iprod
  integer :: nwavel
  integer :: ifnspec,ifnchem,ifnrates,ifnfamilies,ifnstoi,ifnphot
  integer :: ifnoutspec
!  integer :: ifnoutspec,ifnppmsp

  character(len=15),dimension(1000) :: iread1,iread2
!  character(len=15) :: snppm
  character(len=15) :: str
  !*****************************************************************************************

  tabtemp = (/260.,280.,300.,320./)
  !  The number of constants per rate type
  ltabrate(1:25) = (/1,2,3,7,1,4,4,8,8,4,2,8,1,7,1,6,3,4,2,5,4,5,5,5,1/)

!  print*,'  Reading active species names'

  call opfi(ifnspec,fnspec,'f','o',0)
  nread = 0
  inNH3 = 0
  inSO2 = 0
  do ns=1,1000000
    read(ifnspec,*,end=1001)nsp,species(ns)%name,species(ns)%transp,species(ns)%transpv,species(ns)%bounddry
    nread = nread + 1
    if(species(ns)%name.eq.'NH3') inNH3 = nread
    if(species(ns)%name.eq.'SO2') inSO2 = nread
    if(species(ns)%name.eq.'HNO3') inHNO3 = nread
  enddo
1001 continue
  if(nread.ne.nspec) then
    print *,'*** ERROR: Number of active species read in file'
    print *,'           ACTIVE_SPECIES different from NSPEC'
    print *,'           in file chimere_params.h'
    print *,'  ',nread, ' /= ', nspec
   call exit1('Exiting')
  endif
  close(ifnspec)

  !  Prescribed species names (add other here if desired but change
  !  NPRESC in file CHIMERE.H
  !!! Hmmm ... This seems to be inconsistent ! jlm. OK. lmbb

  species(nspec+1)%name = 'M'
  species(nspec+2)%name = 'O2'
  species(nspec+3)%name = 'N2'
  species(nspec+4)%name = 'H2O'
  if(nspresc.ne.4) then
    call exit1('*** ERROR : NSPRESC different from 4')
  endif

!  print*,'  Definition of families'

  call opfi(ifnfamilies,fnfamilies,'f','o',0)
  nread = 0
  do nf=1,nfam
    read(ifnfamilies,*,end=1007)nelem(nf),iread1(1)                &
          &  ,(iread2(ne),ne=1,nelem(nf))
    nread = nread + 1
    species(nspec+nspresc+nf)%name = iread1(1)
    do ne=1,nelem(nf)
        call findspec(iread2(ne),ifam(nf,ne))
        if(ifam(nf,ne).eq.0) then
           print *,'*** ERROR: Family element not given in'
           print *,'           file ACTIVE_SPECIES:',iread2(ne)
           call exit1('Exiting')
        endif
    enddo
  enddo
1007 continue
  close(ifnfamilies)

!  print*,'  Reading output species'

  call opfi(ifnoutspec,fnoutspec,'f','o',0)
  noutspec = 0
  do nsp=1,nspectot
    read(ifnoutspec,*,end=1005)iread1(1),xmass
    call findspec(iread1(1),ns)
    if(ns.gt.0) then
        noutspec = noutspec + 1
        output_species(noutspec)%name  = iread1(1)
        output_species(noutspec)%iaddr = ns
        species(ns)%fspec = xmass !*1.6603e-12
    else
        print *,'* OUTPUT_SPECIES: species ',iread1(1)               &
             &            ,' ignored'
    endif
  enddo
1005 continue
  close(ifnoutspec)

!  print*,'  Reading reaction addressing arrays of chemistry'

  do ns=1,nspectot
    kreacp(ns) = 0
    kreacl(ns) = 0
  enddo

  call opfi(ifnchem,fnchem,'f','o',0)
  
  nread = 0
  do nr=1,nreac

    read(ifnchem,*,end=1003)nreactants(nr)                         &
          &  ,(iread1(nre),nre=1,nreactants(nr))                             &
          &  ,nprods,(iread2(np),np=1,nprods)
    nread = nread + 1

 !   print*,' Addressing the reactant list'

    do it=1,nreactants(nr)
        call findspec(iread1(it),iloss)
        if(iloss.eq.0) then
           print *,'*** ERROR: Reactant species not found in'
           print *,'           ACTIVE_SPECIES file: ',iread1(it)
           call exit1('Exiting')
        endif
        kreacl(iloss) = kreacl(iloss) + 1
        ireacl(iloss,kreacl(iloss)) = nr
        irctt(nr,it) = iloss
    enddo

  !  print*,'  Addressing the product list'

    do np=1,nprods
        call findspec(iread2(np),iprod)
        if(iprod.gt.0) then
           kreacp(iprod) = kreacp(iprod) + 1
           ireacp(iprod,kreacp(iprod)) = nr
        endif
    enddo
  enddo
1003 continue
  close(ifnchem)

! lmbb flag photolysis and stoechiometry
  !  Initialization of stoichiometric coefficients
  do nr=1,nreac
    do ns=1,nspectot
        do nt=1,ntemps
           stoi(ns,nr,nt) = dun
        enddo
    enddo
  enddo
  !print*,'  Reading and addressing stoichiometric coefficients'

  call opfi(ifnstoi,fnstoi,'f','o',0)
  do nr=1,10000000
    read(ifnstoi,*,end=1004)iread1(1),fr,noreac
    call findspec(iread1(1),noprod)
    if(noprod.gt.0) then
        do nt=1,ntemps
           stoi(noprod,noreac,nt) = fr(nt)
        enddo
    endif
  enddo
1004 continue
  close(ifnstoi)

!  print*,'  Reading clear-sky photolysis parameters'

  call opfi(ifnphot,fnphot,'f','o',0)
  read(ifnphot,*)ntabuzen,nphot,nwavel,nlevphot                     &
       &              ,(altiphot(nl),nl=1,nlevphot)
  do nl=1,nlevphot
     altiphot(nl) = 100.*altiphot(nl)
  enddo
  do nt=1,ntabuzen
     !  First 3 lines reserved for future use
     do nw=1,3
        read(ifnphot,*)
     enddo
     !  Clear-sky photolysis
     do np=1,nphot
        read(ifnphot,*)str,zenang(nt),dum1,dum2,dum3                &
             &                    ,(photoj(nt,nl,np),nl=1,nlevphot)
     enddo
  enddo
  close(ifnphot)
  

  !  Calculation of Cosines of zenith angles

  do nt=1,ntabuzen
     zetaref(nt) = cos(zenang(nt)*pi/180d0)
  enddo

 ! print*,'  Initialization of rate values'

  do nr=1,nreac
     do nt=1,ntabmax
        tabrate(nt,nr) = dzero
     enddo
  enddo

  !print*,'  Reading the reaction rates constants'

  call opfi(ifnrates,fnrates,'f','o',0)
  do nr=1,nreac
     iphoto(nr) = 0
     read(ifnrates,*,end=1006)noreac,ity,(tabrate(lt,noreac)        &
          &  ,lt=1,ltabrate(ity))
     ityperate(noreac) = ity

     !  One also precalculates the slopes for comp. efficiency in case of
     !  photolytic reactions

     if(ity.eq.5.or.ity.eq.13) then
        iphoto(nr) = int(tabrate(1,noreac))
     endif
  enddo
1006 continue
  close(ifnrates)

  
end subroutine inichem
