subroutine awrite_adj
  
  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none
  
  !*****************************************************************************
  ! Local variables
  integer :: ncstat
  integer :: ne,ivert,ime,izo,nh,ns
  integer :: ide,idex,ilev,ilath,nl
  
  character(len=dlen) :: datebuf,datebuf2
  
  real(kind=8),dimension(nzonal_domain,nmerid_domain,nlevemis) :: buf3d_a !aemisa
  real(kind=8),dimension(nzonal_domain,nmerid_domain,1) :: buf2d_b !aemisb
  real(kind=8),dimension(nspecbounout,nzonal_domain,nmerid_domain) :: topbuf
  real(kind=8),dimension(nspecbounout,nhbound_domain,nverti) :: latbuf
  real(kind=8),dimension(nzonal_domain,nmerid_domain,nverti) :: buf3d_c ! ainiconc
  real(kind=8) ::  min,max
  
  ! Functions
  character(len=dlen) :: numeric2mm5date
  integer :: mm5date2numeric
  !*****************************************************************************
  
  min=1000.
  max=-9000.
  
  do nh=0,nhourrun
    
    ! TIMES
    datebuf=numeric2mm5date(idate(nh))
    ncstat=nf90_put_var(aoutea_ncid,out_times_varid_ad(1),datebuf,(/1,1+nh/),(/dlen,1/))
    if(optemisb.ne.0) then
      ncstat=nf90_put_var(aouteb_ncid,out_times_varid_ad(4),datebuf,(/1,1+nh/),(/dlen,1/))
    endif
    ncstat=nf90_put_var(aoutbc_ncid,out_times_varid_ad(2),datebuf,(/1,1+nh/),(/dlen,1/))
    
    ! emissions
    do ne=1,nemisa
      do ivert=1,nlevemis
        do ime=1,nmerid_domain
          do izo=1,nzonal_domain
            buf3d_a(izo,ime,ivert)=aemisaloc(nh,ne,izo,ime,ivert)
            if (buf3d_a(izo,ime,ivert).lt.min) then
              if(buf3d_a(izo,ime,ivert).ne.0.) min=buf3d_a(izo,ime,ivert)
            endif
            if (buf3d_a(izo,ime,ivert).gt.max) max=buf3d_a(izo,ime,ivert)
          
          end do
        end do
      end do
      ncstat=nf90_put_var(aoutea_ncid,emisa_species(ne)%varid,buf3d_a,&
              (/1,1,1,nh+1/),(/nzonal_domain,nmerid_domain,nlevemis,1/))
      !if(any(buf3d_a(:,:,1) .ne. 0.)) then
      !    print 101,emisa_species(ne)%name,nh 
      !    do ime=1,nmerid_domain
      !    do izo=1,51 !nzonal_domain
      !      print 102,izo,ime,buf3d_a(izo,ime,1)
      !    enddo
      !    enddo
      !endif
!101 format('AAAAAAAAAAAAA ',a23,' ',i2) !,' ',8585 x(e65.56,' '))  
!102 format(i3,' ',i3,' ',e65.56)
   end do
    
    if (optemisb.ne.0) then
      do ne=1,nemisb
        do ime=1,nmerid_domain
          do izo=1,nzonal_domain
            buf2d_b(izo,ime,1)=aemisb(nh,ne,izo,ime)
          enddo
        enddo
        ncstat=nf90_put_var(aouteb_ncid,emisb_species(ne)%varid,buf2d_b,&
                (/1,1,1,nh+1/),(/nzonal_domain,nmerid_domain,1,1/))
      enddo
    endif
    
    ! top conc
    min=1000.
    max=-9000.
    ne = 0
    do ns=1,nspecboun
      if(isboun(ns).gt.0) then
        ne = ne + 1
        do ime=1,nmerid_domain
          do izo=1,nzonal_domain
            topbuf(ne,izo,ime)=aboundtop(izo,ime,isboun(ns),nh)
          enddo
        enddo
      endif
    enddo
    ncstat=nf90_put_var(aoutbc_ncid,atopconc_conc_varid,topbuf,&
            (/1,1,1,nh+1/),(/nspecbounout,nzonal_domain,nmerid_domain,1/))
    
    ! lat conc
    min=1000.
    max=-9000.
    ne = 0
    do ns=1,nspecboun
      if(isboun(ns).gt.0) then
        ne = ne + 1
        nl=0
        do ivert=1,nverti
          ilath=0
          
          do izo=0,nzonal_domain+1,nzonal_domain+1
            do ime=1,nmerid_domain
              nl=nl+1
              ilath=ilath+1
              latbuf(ne,ilath,ivert)=aboundlat(nl,isboun(ns),nh)
            enddo
          enddo
          
          do ime=0,nmerid_domain+1,nmerid_domain+1
            do izo=1,nzonal_domain
              nl=nl+1
              ilath=ilath+1
              latbuf(ne,ilath,ivert)=aboundlat(nl,isboun(ns),nh)
            enddo
          enddo
        
        enddo ! ivert
      end if
    enddo !ns
    ncstat=nf90_put_var(aoutbc_ncid,alatconc_conc_varid,latbuf,&
            (/1,1,1,nh+1/),(/nspecbounout,nhbound_domain,nverti,1/))
    
    ! Synchronize disk to avoid data loss
    ncstat=nf90_sync(aoutea_ncid)
    if(optemisb.ne.0)ncstat=nf90_sync(aouteb_ncid)
    ncstat=nf90_sync(aoutbc_ncid)
  
  enddo ! nh
  
  datebuf=numeric2mm5date(idate(nhourrun))
  ncstat=nf90_put_var(aoutini_ncid,out_times_varid_ad(3),datebuf,(/1,1/),(/dlen,1/))
  ne = 0
  do ns=1,nspecboun
    if(isboun(ns).gt.0) then
      ne = ne + 1
      min=1000.
      max=-9000.
      do ivert=1,nverti
        do ime=1,nmerid_domain
          do izo=1,nzonal_domain
            buf3d_c(izo,ime,ivert)=aconcini(isboun(ns),izo,ime,ivert)
          end do
        end do
      end do !ivert
      
      ncstat=nf90_put_var(                                 &
              aoutini_ncid,                                       &
              aconc_species(ns)%varid,                           &
              buf3d_c,                                          &
              (/1,1,1/),(/nzonal_domain,nmerid_domain,nverti/) &
              )
    end if
  end do !ns
  
  ! Synchronize disk to avoid data loss
  ncstat=nf90_sync(aoutini_ncid)

end subroutine awrite_adj
