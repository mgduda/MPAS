#define DM_BCAST_MACRO(A) call mpas_dmpar_bcast_reals(dminfo,size(A),A)

!=================================================================================================================
 module mpas_atmphys_gfs_grrad_init
 use mpas_dmpar
 use mpas_kind_types
 use mpas_grid_types
 use mpas_timekeeping

 use mpas_atmphys_gfs_vars

!GFS physics:
 use module_radiation_aerosols, only: aerinit
 use module_radiation_gases, only: gasinit
!use module_radiation_astronomy, only: solinit
 use module_radiation_clouds, only: cldinit
 use module_radiation_surface, only: sfcinit
 use module_radlw_main, only: rlwinit
 use module_radsw_main, only: rswinit

 implicit none
 private
 public:: radiation_init_forMPAS
 

!>\brief MPAS main initialization subroutine for the GFS radiation codes.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-01-27.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine radiation_init_forMPAS(dminfo,clock,mesh)
!=================================================================================================================

!input arguments:
 type(dm_info):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 type(MPAS_Clock_type),intent(in):: clock

!local pointers:
 integer,pointer:: nVertLevels

!local variables:
 type(MPAS_Time_Type):: currTime
 character(len=StrKind):: timeStamp
 integer:: me
 integer:: ierr,iyear,imonth,imonth2,iday,ihour

 integer:: nVertLevels
 real(kind=RKIND),dimension(:),allocatable:: sigma

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine radiation_init_forMPAS:'

 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

!... get the current year, month, day, and hour:
!100 format('     YEAR        =', i5   ,/, &
!           '     MONTH       =', i5   ,/, &
!           '     DAY         =', f16.9,/, &
!           '     UTC_H       =', f16.9,/, &
!           '     CURR_JULDAY =', f16.9,/, &
!           '     LEAP_YEAR   =', 1x,l1,/, &
!           '     TIME STAMP  =', 1x,a32,/)

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=iyear,MM=imonth, &
                    M=imonth2,DD=iday,H=ihour,ierr=ierr)
 write(0,*) 'timeStamp =', trim(timeStamp)
 write(0,*) 'iyear     =', iyear
 write(0,*) 'imonth    =', imonth
 write(0,*) 'imonth2   =', imonth2
 write(0,*) 'iday      =', iday
 write(0,*) 'ihour     =', ihour

!if(dminfo % my_proc_id == IO_NODE) me = 0
 me = dminfo % my_proc_id

!... initialization of surface albedo and surface emissivity:
 write(0,*) '--- enter subroutine sfcinit: me =',me
 call sfcinit(nVertLevels,iflip,ialb,iems,me)
 write(0,*) '--- end subroutine sfcinit:'
 write(0,*)

!... initialization of cloud-radiation calculations. set-up of boundary layer cloud-top.
 if(.not.allocated (sigma)) allocate(sigma(nVertLevels+1))
 sigma(1:nVertLevels) = 0._RKIND
 write(0,*) '--- enter subroutine cldinit: me =',me
 call cldinit(sigma,nVertLevels,iflip,np3d,icld,me)
 write(0,*) '--- end subroutine cldinit'
 write(0,*)
 if(allocated(sigma)) deallocate(sigma)

!... initialization of parameterization of longwave radiation:
 write(0,*) '--- enter subroutine rlwinit: me=',me
 call rlwinit(icld,me,nVertLevels,iovr_lw,isubclw)
 write(0,*) '--- end subroutine rlwinit'
 write(0,*)

!... initialization of parameterization of shortwave radiation:
 write(0,*) '--- enter subroutine rswinit: me=',me
 call rswinit(icld,me,nVertLevels,iovr_sw,isubcsw)
 write(0,*) '--- end subroutine rswinit'
 write(0,*)

!... initialization of solar constant:
!call solinit(isol,iyear,iyear,me)

!... initialization of aerosol effects:
!write(0,*) '--- enter subroutine aerinit: me=',me
!call aerinit(iyear,imon,iaer,me)

!... initialization of trace gases effects:
!write(0,*) '--- enter subroutine gasinit:'
!call gasinit(....)

 write(0,*) '--- end subroutine radiation_init_forMPAS'
 write(0,*)

 end subroutine radiation_init_forMPAS

!=================================================================================================================
 end module mpas_atmphys_gfs_grrad_init
!=================================================================================================================
