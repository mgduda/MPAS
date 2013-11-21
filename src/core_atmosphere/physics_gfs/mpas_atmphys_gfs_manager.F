! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_manager
 use mpas_configure
 use mpas_kind_types
 use mpas_grid_types
 use mpas_timekeeping

 use mpas_atmphys_utilities
 use mpas_atmphys_gfs_update_surface
 use mpas_atmphys_gfs_vars

 implicit none
 private
 public:: physics_timetracker,physics_run_init

 integer, public:: year                 !Current year.
 integer, public:: month                !Current month.
 integer, public:: day                  !Current day.
 integer, public:: hour                 !Current hour.
 integer, public:: julday               !Initial Julian day.
 real(kind=RKIND), public:: curr_julday !Current Julian day (= 0.0 at 0Z on January 1st).
 real(kind=RKIND), public:: gmt         !Greenwich mean time hour of model start (hr)
 real(kind=RKIND), public:: utc_h       !UTC time in hours.

!defines alarms for calling the long- and short-wave radiation codes:
 integer, parameter:: radtlwAlarmID = 11
 integer, parameter:: radtswAlarmID = 12

!defines alarm to update the background surface albedo and the greeness fraction:
 integer, parameter:: greenAlarmID  = 16

 integer :: h, d, dd, m, mm, s, s_n, s_d, DoY, yr

!>\brief main time manager for GFS physics parameterizations.
!>\author Laura D. Fowler (send comments to laura@mmm.ucar.edu).
!>\date 2013-11-20.
!>
!>\details
!>
!> subroutines in mpas_atmphys_gfs_manager:
!> ----------------------------------------
!> physics_timetracker: check alarms and update boundary conditions if needed.
!> physics_run_init   : create and initialize alarms used for GFS physics parameterizations.
!>


 contains

!==================================================================================================
 subroutine physics_timetracker(domain,dt,clock,itimestep,xtime_s)
!==================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(MPAS_Clock_type),intent(inout):: clock
 type(domain_type),intent(inout):: domain

!local variables:
 type(block_type),pointer:: block
 type(MPAS_Time_Type):: currTime

 logical:: LeapYear
 character(len=StrKIND)   :: timeStamp
 integer:: ierr
 real(kind=RKIND):: utc_s
 real(kind=RKIND):: xtime_m

!==================================================================================================
 write(0,*)
 write(0,*) '--- enter subroutine physics_timetracker: itimestep = ', itimestep

!update the current julian day and current year:
 100 format('     YEAR        =', i5   ,/, &
            '     MONTH       =', i5   ,/, &
            '     DAY         =', i5   ,/, &
            '     HOUR        =', i5   ,/, &
            '     JULDAY      =', i5   ,/, &
            '     GMT         =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     CURR_JULDAY =', f16.9,/, &
            '     LEAP_YEAR   =', 1x,l1,/, &
            '     TIME STAMP  =', 1x,a32,/)

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,MM=mm,DD=d,H=h, &
                    M=m,S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)

 utc_h = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 utc_s = real(s + s_n + s_d)
 year   = yr
 month  = mm
 day    = d
 hour   = h
 julday = DoY
 curr_julday = real(julday-1) + utc_h / 24.0
 LeapYear = isLeapYear(year)
 write(0,100) year,month,day,hour,julday,gmt,utc_h,curr_julday,LeapYear,timeStamp

 block => domain % blocklist
 do while(associated(block))

    !update the background surface albedo and greeness of vegetation: interpolation of input
    !monthly values to current day:
    if(mpas_is_alarm_ringing(clock,greenAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,greenAlarmID,ierr=ierr)
       write(0,*) '--- time to update background surface albedo, greeness fraction.'
       call physics_gfs_update_surface(timeStamp,block%mesh,block%sfc_input,block%diag_physics)
    endif

    block => block % next
 end do

!check to see if it is time to run the longwave and shortwave radiation codes:
 l_radtlw = .false.

 if(config_radtlw_interval /= "none") then
    if(mpas_is_alarm_ringing(clock,radtlwAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,radtlwAlarmID,ierr=ierr)
       l_radtlw = .true.
    endif
 elseif(config_radtlw_interval == "none") then
    l_radtlw = .true.
 endif
 write(0,*) '--- time to run the LW radiation scheme L_RADLW =',l_radtlw

 l_radtsw = .false.
 if(config_radtsw_interval /= "none") then
    if(mpas_is_alarm_ringing(clock,radtswAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,radtswAlarmID,ierr=ierr)
       l_radtsw = .true.
    endif
 elseif(config_radtsw_interval == "none") then
    l_radtsw = .true.
 endif
 write(0,*) '--- time to run the SW radiation scheme L_RADSW =',l_radtsw

 write(0,*) '--- end subroutine physics_timetracker.'

 end subroutine physics_timetracker

!==================================================================================================
 subroutine physics_run_init(mesh,state,clock)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh
 type(state_type),intent(in):: state
 type(MPAS_Clock_type):: clock

!local variables:
 type(MPAS_Time_Type):: startTime,alarmStartTime
 type(MPAS_TimeInterval_Type):: alarmTimeStep, alarmTimeStepHi
 integer:: ierr

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_run_init'

!initialization of gmt, julian day, and alarms:
 101 format('     YEAR        =', i5   ,/, &
            '     JULDAY      =', i5   ,/, &
            '     GMT         =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     CURR_JULDAY =', f16.9,/)

 102 format('     DT_RADTLW   =', f16.9,/, &
            '     DT_RADTSW   =', f16.9,/)

 startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
 call mpas_get_time(curr_time=startTime,YYYY=yr,H=h,M=m,S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)
 utc_h  = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 year   = yr
 gmt    = utc_h
 julday = DoY
 curr_julday = real(julday-1)
 write(0,101) year,julday,gmt,utc_h,curr_julday

!set alarms for calling the longwave and shortwave radiation schemes, the convection schemes,
!and the PBL schemes at intervals different (greater) than the dynamical time-step:

 if(trim(config_radtlw_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_radtlw_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,radtlwAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error creating radtlwAlarmID')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_radtlw,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining dt_radtlw')

 elseif(trim(config_radtlw_interval) == "none") then
    dt_radtlw = config_dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_radtlw is not defined')

 endif

 if(trim(config_radtsw_interval) /= "none") then
    call mpas_set_timeInterval(alarmTimeStep,timeString=config_radtsw_interval,ierr=ierr)
    alarmStartTime = startTime
    call mpas_add_clock_alarm(clock,radtswAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error creating alarm radtsw')

    call mpas_get_timeInterval(interval=alarmTimeStep,dt=dt_radtsw,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_run_init: error defining radtswAlarmID')

 elseif(trim(config_radtsw_interval) == "none") then
    dt_radtsw = config_dt

 else
    call physics_error_fatal('subroutine physics_run_init: dt_radtsw is not defined')

 endif

!set alarm for updating the background surface albedo and the greeness fraction:
 call mpas_set_timeInterval(alarmTimeStep,timeString=config_greeness_update,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,greenAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm greeness')

 write(0,102) dt_radtlw,dt_radtsw

!initialization of local physics time-steps:
!... dynamics:
 dt_dyn     = config_dt
!... radiation:
 l_radtlw = .false.
 l_radtsw = .false.

 write(0,*) '--- end subroutine physics_run_init.'

 end subroutine physics_run_init

!==================================================================================================
 end module mpas_atmphys_manager
!==================================================================================================
 
