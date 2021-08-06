subroutine write_depo


  use netcdf
  use chimere_consts
  use chimere_common
  use wholedomain_common
  implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'write_depo.f90')

  !*****************************************************************************************
  ! Local variables
  integer :: ncstat
  integer :: ispec,ime,izo
  integer :: tsavedepos

  character(len=dlen) :: datebuf

  real,dimension(nzonal_domain,nmerid_domain) :: buf2d

  ! Functions
  character(len=dlen) :: numeric2mm5date
  !*****************************************************************************************


  ! TIMES
  datebuf=numeric2mm5date(idate(ihourrun))
  tsavedepos=ihourrun/nsavedepos
  ncstat=nf90_put_var(depo_ncid,depo_times_varid,datebuf,(/1,tsavedepos/),(/dlen,1/))
  NCERR(__LINE__)


  ! DEPOSITIONS
  do ispec=1,nspec
     do ime=1,nmerid_domain
        do izo=1,nzonal_domain
           buf2d(izo,ime)=drydep(ispec,izo,ime)
        end do
     end do
     ncstat=nf90_put_var(                                 &
          depo_ncid,                                      &
          output_depo(ispec)%varid,                       &
          buf2d,                                          &
          (/1,1,tsavedepos/),(/nzonal_domain,nmerid_domain,1/) &
          )
     NCERR(__LINE__)
  end do
  do ispec=nspec+1,2*nspec
     do ime=1,nmerid_domain
        do izo=1,nzonal_domain
           buf2d(izo,ime)=wetdep(ispec-nspec,izo,ime)
        end do
     end do
     ncstat=nf90_put_var(                                 &
          depo_ncid,                                      &
          output_depo(ispec)%varid,                       &
          buf2d,                                          &
          (/1,1,tsavedepos/),(/nzonal_domain,nmerid_domain,1/) &
          )
     NCERR(__LINE__)
  end do

  ! Synchronize disk to avoid data loss
  ncstat=nf90_sync(depo_ncid)
  NCERR(__LINE__)

end subroutine write_depo
