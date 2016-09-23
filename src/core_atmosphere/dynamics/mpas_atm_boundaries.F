module mpas_atm_boundaries

    use mpas_derived_types, only : mpas_pool_type, mpas_clock_type, block_type, mpas_time_type, mpas_timeInterval_type, MPAS_NOW, &
                                   MPAS_STREAM_LATEST_BEFORE, MPAS_STREAM_EARLIEST_STRICTLY_AFTER, MPAS_STREAM_INPUT, &
                                   MPAS_streamManager_type
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_subpool, mpas_pool_shift_time_levels
    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_geometry_utils, only : mpas_arc_length
    use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_time, mpas_get_timeinterval, mpas_set_time, operator(-)
    use mpas_stream_manager, only : mpas_stream_mgr_read

    type (MPAS_Time_Type) :: LBC_intv_end


    contains


    subroutine mpas_atm_setup_bdy_tend(clock, block)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (block_type), intent(inout) :: block

        type (mpas_pool_type), pointer :: mesh
        type (mpas_pool_type), pointer :: state
        type (mpas_pool_type), pointer :: lbc

        real (kind=RKIND), dimension(:,:), pointer :: u
        real (kind=RKIND), dimension(:,:), pointer :: w
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_theta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_scalars

        integer :: ierr


write(0,*) 'Setting up next tendency'
        call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        call mpas_pool_get_array(state, 'u', u, 1)
        call mpas_pool_get_array(state, 'w', w, 1)
        call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
        call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
        call mpas_pool_get_array(state, 'scalars', scalars, 1)

        call mpas_pool_get_array(lbc, 'lbc_u', lbc_u, 1)
        call mpas_pool_get_array(lbc, 'lbc_w', lbc_w, 1)
        call mpas_pool_get_array(lbc, 'lbc_theta_m', lbc_theta_m, 1)
        call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_rho_zz, 1)
        call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_scalars, 1)

        lbc_u(:,:) = u(:,:)
        lbc_w(:,:) = w(:,:)
        lbc_theta_m(:,:) = theta_m(:,:)
        lbc_rho_zz(:,:) = rho_zz(:,:)
        lbc_scalars(:,:,:) = scalars(:,:,:)

    end subroutine mpas_atm_setup_bdy_tend


    subroutine mpas_atm_update_bdy_tend(clock, streamManager, block, firstCall)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (MPAS_streamManager_type), intent(inout) :: streamManager
        type (block_type), intent(inout) :: block
        logical, intent(in) :: firstCall

        type (mpas_pool_type), pointer :: mesh
        type (mpas_pool_type), pointer :: state
        type (mpas_pool_type), pointer :: lbc
        real (kind=RKIND) :: dt

        real (kind=RKIND), dimension(:,:), pointer :: u
        real (kind=RKIND), dimension(:,:), pointer :: w
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

        integer :: ierr
        integer :: dd_intv, s_intv, sn_intv, sd_intv
        type (MPAS_Time_Type) :: currTime
        type (MPAS_TimeInterval_Type) :: lbc_interval
        character(len=StrKIND) :: read_time


        call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        if (firstCall) then
            call MPAS_stream_mgr_read(streamManager, streamID='lbc_in', timeLevel=2, whence=MPAS_STREAM_LATEST_BEFORE, actualWhen=read_time, ierr=ierr)
        else
            call mpas_pool_shift_time_levels(lbc)
            call MPAS_stream_mgr_read(streamManager, streamID='lbc_in', timeLevel=2, whence=MPAS_STREAM_EARLIEST_STRICTLY_AFTER, actualWhen=read_time, ierr=ierr)
        end if
        call mpas_set_time(currTime, dateTimeString=trim(read_time))

        if (.not. firstCall) then
            lbc_interval = currTime - LBC_intv_end
            call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
            dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))

write(0,*) 'Computing tendency'
write(0,*) '    dt = ', dt

            call mpas_pool_get_array(lbc, 'lbc_u', u, 2)
            call mpas_pool_get_array(lbc, 'lbc_w', w, 2)
            call mpas_pool_get_array(lbc, 'lbc_theta_m', theta_m, 2)
            call mpas_pool_get_array(lbc, 'lbc_rho_zz', rho_zz, 2)
            call mpas_pool_get_array(lbc, 'lbc_scalars', scalars, 2)

            call mpas_pool_get_array(lbc, 'lbc_u', lbc_tend_u, 1)
            call mpas_pool_get_array(lbc, 'lbc_w', lbc_tend_w, 1)
            call mpas_pool_get_array(lbc, 'lbc_theta_m', lbc_tend_theta_m, 1)
            call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_tend_rho_zz, 1)
            call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_tend_scalars, 1)

            dt = 1.0_RKIND / dt
            lbc_tend_u(:,:) = (u(:,:) - lbc_tend_u(:,:)) * dt
            lbc_tend_w(:,:) = (w(:,:) - lbc_tend_w(:,:)) * dt
            lbc_tend_theta_m(:,:) = (theta_m(:,:) - lbc_tend_theta_m(:,:)) * dt
            lbc_tend_rho_zz(:,:) = (rho_zz(:,:) - lbc_tend_rho_zz(:,:)) * dt
            lbc_tend_scalars(:,:,:) = (scalars(:,:,:) - lbc_tend_scalars(:,:,:)) * dt
        end if

        LBC_intv_end = currTime

    end subroutine mpas_atm_update_bdy_tend


    subroutine mpas_atm_get_bdy_tend(clock, block)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (block_type), intent(inout) :: block

        type (mpas_pool_type), pointer :: state
        type (mpas_pool_type), pointer :: lbc


        real (kind=RKIND), dimension(:,:), pointer :: u
        real (kind=RKIND), dimension(:,:), pointer :: w
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_theta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

        type (MPAS_Time_Type) :: currTime
        type (MPAS_TimeInterval_Type) :: lbc_interval

        integer :: dd_intv, s_intv, sn_intv, sd_intv
        real (kind=RKIND) :: dt
        integer :: ierr


        currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
!        call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)

        lbc_interval = LBC_intv_end - currTime

        call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
        dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))


        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        call mpas_pool_get_array(state, 'u', u, 1)
        call mpas_pool_get_array(state, 'w', w, 1)
        call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
        call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
        call mpas_pool_get_array(state, 'scalars', scalars, 1)

        call mpas_pool_get_array(lbc, 'lbc_u', lbc_tend_u, 1)
        call mpas_pool_get_array(lbc, 'lbc_w', lbc_tend_w, 1)
        call mpas_pool_get_array(lbc, 'lbc_theta_m', lbc_tend_theta_m, 1)
        call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_tend_rho_zz, 1)
        call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_tend_scalars, 1)

        call mpas_pool_get_array(lbc, 'lbc_u', lbc_u, 2)
        call mpas_pool_get_array(lbc, 'lbc_w', lbc_w, 2)
        call mpas_pool_get_array(lbc, 'lbc_theta_m', lbc_theta_m, 2)
        call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_rho_zz, 2)
        call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_scalars, 2)

        u(:,:) = lbc_u(:,:) - dt * lbc_tend_u(:,:)
        w(:,:) = lbc_w(:,:) - dt * lbc_tend_w(:,:)
        theta_m(:,:) = lbc_theta_m(:,:) - dt * lbc_tend_theta_m(:,:)
        rho_zz(:,:) = lbc_rho_zz(:,:) - dt * lbc_tend_rho_zz(:,:)
        scalars(:,:,:) = lbc_scalars(:,:,:) - dt * lbc_tend_scalars(:,:,:)

    end subroutine mpas_atm_get_bdy_tend

end module mpas_atm_boundaries
