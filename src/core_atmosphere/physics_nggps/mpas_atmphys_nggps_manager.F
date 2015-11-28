! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_manager
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_timekeeping

 use mpas_atmphys_utilities
 use mpas_atmphys_nggps_o3climatology
 use mpas_atmphys_update
 use mpas_atmphys_nggps_update_surface
 use mpas_atmphys_nggps_vars

 implicit none
 private
 public:: physics_timetracker,physics_run_init

 character(len=StrKIND),public:: physics_timeStamp
 integer,public:: year                 !Current year.
 integer,public:: month                !Current month.
 integer,public:: day                  !Current day.
 integer,public:: hour                 !Current hour.
 integer,public:: minute               !Current minute.
 integer,public:: second               !Current second.
 integer,public:: msecond              !Current milli-second.
 integer,public:: julday               !Initial Julian day.
 real(kind=RKIND),public:: curr_julday !Current Julian day (= 0.0 at 0Z on January 1st).
 real(kind=RKIND),public:: gmt         !Greenwich mean time hour of model start (hr)
 real(kind=RKIND),public:: utc_h       !UTC time in hours.

!defines alarms for calling the long- and short-wave radiation codes:
 character(len=*),parameter:: radtlwAlarmID = 'radtlw'
 character(len=*),parameter:: radtswAlarmID = 'radtsw'

!defines alarm to update the background surface albedo and the greeness fraction:
 character(len=*),parameter:: greenAlarmID  = 'green'

!defines alarm to reset to zero the radiation and other physics fields accumulated in subroutine
!gbphys. By default, these fields are initialized every 6 hours, as with the GFS.
 character(len=*), parameter:: acfieldsAlarmID = 'acfields'

 integer:: h, d, dd, m, mm, s, s_n, s_d, DoY, yr


!modifications to sourcecode:
!----------------------------
! * In subroutine physics_timetracker, renamed and made public the variable timeStamp to physics_timeStamp.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-28.


 contains


!=================================================================================================================
 subroutine physics_timetracker(domain,dt,clock,itimestep,xtime_s)
!=================================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(MPAS_Clock_type),intent(inout):: clock
 type(domain_type),intent(inout):: domain

!local pointers:
 character(len=StrKIND),pointer:: config_radtlw_interval, &
                                  config_radtsw_interval

 type(block_type),pointer    :: block
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: nggps_diagnostics
 type(mpas_pool_type),pointer:: nggps_sfc_properties
 type(mpas_pool_type),pointer:: nggps_atm_inputs

!local variables:
 type(MPAS_Time_Type):: currTime

 logical:: LeapYear
 integer:: ierr
 real(kind=RKIND):: utc_s
 real(kind=RKIND):: xtime_m

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_timetracker: itimestep = ', itimestep

 call mpas_pool_get_config(domain%blocklist%configs,'config_radtlw_interval',config_radtlw_interval)
 call mpas_pool_get_config(domain%blocklist%configs,'config_radtsw_interval',config_radtsw_interval)

!update the current julian day and current year:
 100 format('     YEAR        =', i5   ,/, &
            '     MONTH       =', i5   ,/, &
            '     DAY         =', i5   ,/, &
            '     HOUR        =', i5   ,/, &
            '     MINUTE      =', i5   ,/, &
            '     SECOND      =', i5   ,/, &
            '     MSECOND     =', i5   ,/, &
            '     JULDAY      =', i5   ,/, &
            '     GMT         =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     CURR_JULDAY =', f16.9,/, &
            '     LEAP_YEAR   =', 1x,l1,/, &
            '     TIME STAMP  =', 1x,a32,/)

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=physics_timeStamp,YYYY=yr,MM=mm,DD=d,H=h, &
                    M=m,S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)

 utc_h = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 utc_s = real(s + s_n + s_d)
 year    = yr
 month   = mm
 day     = d
 hour    = h
 minute  = m
 second  = s
 msecond = 0
 julday  = DoY
 curr_julday = real(julday-1) + utc_h / 24.0
 LeapYear = isLeapYear(year)
 write(0,100) year,month,day,hour,minute,second,msecond,julday,gmt,utc_h,curr_julday, &
              LeapYear,physics_timeStamp

 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'                ,mesh                )
    call mpas_pool_get_subpool(block%structs,'nggps_sfc_properties',nggps_sfc_properties)
    call mpas_pool_get_subpool(block%structs,'nggps_diagnostics'   ,nggps_diagnostics   )

    !initialize diagnostics accumulated in subroutine gbphys back to zero:
    if(mpas_is_alarm_ringing(clock,acfieldsAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,acfieldsAlarmID,ierr=ierr)
       write(0,*) '--- time to initialize diagnostics accumulated in gbphys.'
       call physics_nggps_update_diagnostics(mesh,nggps_diagnostics)
    endif

    !update the background surface albedo and greeness of vegetation: interpolation of input
    !monthly values to current day:
    !if(mpas_is_alarm_ringing(clock,greenAlarmID,ierr=ierr)) then
    !   call mpas_reset_clock_alarm(clock,greenAlarmID,ierr=ierr)
    !   write(0,*) '--- time to update background surface albedo, greeness fraction.'
    !   call physics_nggps_update_surface(timeStamp,mesh,nggps_sfc_properties)
    !endif

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

!update ozone to the current julian day for use in the parameterization of !progsnotic ozone:
 
 block => domain % blocklist
 do while(associated(block))
    call mpas_pool_get_subpool(block%structs,'mesh'            ,mesh            )
    call mpas_pool_get_subpool(block%structs,'nggps_atm_inputs',nggps_atm_inputs)

    call o3climatology_from_MPAS(mesh,nggps_atm_inputs)

    block => block % next
 end do

 write(0,*) '--- end subroutine physics_timetracker.'

 end subroutine physics_timetracker

!=================================================================================================================
 subroutine physics_run_init(configs,mesh,state,clock,stream_manager)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(MPAS_Clock_type):: clock
 type (MPAS_streamManager_type), intent(inout) :: stream_manager

!local pointers:
 character(len=StrKIND),pointer:: config_radtlw_interval, &
                                  config_radtsw_interval, &
                                  config_bucket_update,   &
                                  config_greeness_update

 real(kind=RKIND),pointer:: config_dt

!local variables:
 type(MPAS_Time_Type):: startTime,alarmStartTime
 type(MPAS_TimeInterval_Type):: alarmTimeStep, alarmTimeStepHi
 integer:: ierr

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_run_init NGGPS:'

 call mpas_pool_get_config(configs,'config_radtlw_interval',config_radtlw_interval)
 call mpas_pool_get_config(configs,'config_radtsw_interval',config_radtsw_interval)
 call mpas_pool_get_config(configs,'config_bucket_update'  ,config_bucket_update  )
 call mpas_pool_get_config(configs,'config_greeness_update',config_greeness_update)

 call mpas_pool_get_config(configs,'config_dt',config_dt)

!initialization of gmt, julian day, and alarms:
 101 format('     YEAR        =', i5   ,/, &
            '     JULDAY      =', i5   ,/, &
            '     GMT         =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     CURR_JULDAY =', f16.9,/)

 102 format('     DT_RADTLW   =', f16.9,/, &
            '     DT_RADTSW   =', f16.9  )

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
 write(0,102) dt_radtlw,dt_radtsw

!set alarm for updating the background surface albedo and the greeness fraction:
 call mpas_set_timeInterval(alarmTimeStep,timeString=config_greeness_update,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,greenAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm greeness')

!set alarm to reset to zero the radiation and other physics fields accumulated in subroutine
!gbphys.
 call mpas_set_timeInterval(alarmTimeStep,timeString=config_bucket_update,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,acfieldsAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine physics_init: error creating alarm acfields')

!initialization of local physics time-steps:
!... dynamics:
 dt_dyn     = config_dt
!... radiation:
 l_radtlw = .false.
 l_radtsw = .false.

 write(0,*) '--- end subroutine physics_run_init NGGPS.'

 end subroutine physics_run_init

!=================================================================================================================
 end module mpas_atmphys_manager
!=================================================================================================================
 
