subroutine awrite_aconcs


    use netcdf
    use chimere_consts
    use chimere_common
    use wholedomain_common

    implicit none

#define NCERR(lnum) if(ncstat/=NF90_NOERR) call nc_err(ncstat,lnum,'write_concs.f90')

    !*****************************************************************************************
    ! Local variables
    integer :: ncstat
    integer :: ispec,ivert,ime,izo
    integer :: tsaveconcs

    character(len=dlen) :: datebuf

    real(kind=8),dimension(nzonal_domain+(6*nzdoms),nmerid_domain+(6*nmdoms),1:nverti+1) :: buf3d
    real(kind=8),dimension(nzonal_domain,nmerid_domain,nverti) :: buf3do


    ! Functions
    character(len=dlen) :: numeric2mm5date
    !*****************************************************************************************


    ! TIMES
    datebuf=numeric2mm5date(idate(ihourrun))
    ncstat=nf90_put_var(aend_ncid,end_times_varid,datebuf,(/1,1/),(/dlen,1/))
NCERR(__LINE__)


    !print*,' ACONCENTRATIONS'
    do ispec=1,nspectot
        do ivert=1,nverti+1
            do ime=1,nmerid_domain+(6*nmdoms)
                do izo=1,nzonal_domain+(6*nzdoms)
                    buf3d(izo,ime,ivert)=aconcsave(ispec,izo,ime,ivert)
                end do
            end do
        end do
        ncstat=nf90_put_var(                                 &
                aend_ncid,                                       &
                species(ispec)%varid,                           &
                buf3d,                                          &
                (/1,1,1,1/),(/nzonal_domain+(6*nzdoms),nmerid_domain+(6*nmdoms),nverti+1,1/) &
                )
        NCERR(__LINE__)
    end do
    !print*,' ACONCENTRATIONS O'
    do ispec=1,nspectot
        do ivert=1,nverti
            do ime=1,nmerid_domain
                do izo=1,nzonal_domain
                    buf3do(izo,ime,ivert)=aconco(ispec,izo,ime,ivert)
                end do
            end do
        end do
        ncstat=nf90_put_var(                                 &
                aend_ncid,                                       &
                species(ispec)%varido,                           &
                buf3do,                                          &
                (/1,1,1,1/),(/nzonal_domain,nmerid_domain,nverti,1/) &
                )
        NCERR(__LINE__)
    end do

    ! Synchronize disk to avoid data loss
    ncstat=nf90_sync(aend_ncid)
NCERR(__LINE__)
end subroutine awrite_aconcs
