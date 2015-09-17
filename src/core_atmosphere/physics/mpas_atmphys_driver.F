! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_driver_cloudiness
 use mpas_atmphys_driver_convection
 use mpas_atmphys_driver_gwdo
 use mpas_atmphys_driver_lsm
 use mpas_atmphys_driver_pbl
 use mpas_atmphys_driver_radiation_lw
 use mpas_atmphys_driver_radiation_sw 
 use mpas_atmphys_driver_sfclayer
 use mpas_atmphys_constants
 use mpas_atmphys_interface
 use mpas_atmphys_update
 use mpas_atmphys_vars, only: l_camlw,l_conv,l_radtlw,l_radtsw
 use mpas_timer

 implicit none
 private
 public:: physics_driver


!>\brief MPAS top physics driver.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!> subroutine physics_driver is the top physics driver from which separate drivers for all physics
!> parameterizations, except cloud microphysics parameterizations are called.
!>
!> subroutines called in mpas_atmphys_driver:
!> ------------------------------------------
!> allocate_forall_physics     : allocate local arrays defining atmospheric soundings (pressure,..)
!> allocate_cloudiness         : allocate all local arrays used in driver_cloudiness. 
!> allocate_convection         : allocate all local arrays used in driver_convection.
!> allocate_gwdo               : allocate all local arrays used in driver_gwdo.
!> allocate_lsm                : allocate all local arrays used in driver_lsm.
!> allocate_pbl                : allocate all local arrays used in driver_pbl.
!> allocate_radiation_lw       : allocate all local arrays used in driver_radiation_lw.
!> allocate_radiation_sw       : allocate all local arrays used in driver_radiation_sw.
!> allocate_sfclayer           : allocate all local arrays used in driver_sfclayer.
!>
!> deallocate_forall_physics   : deallocate local arrays defining atmospheric soundings.
!> deallocate_cloudiness       : dedeallocate all local arrays used in driver_cloudiness.
!> deallocate_convection       : deallocate all local arrays used in driver_convection.
!> deallocate_gwdo             : deallocate all local arrays used in driver_gwdo.
!> deallocate_lsm              : deallocate all local arrays used in driver_lsm.
!> deallocate_pbl              : deallocate all local arrays used in driver_pbl.
!> deallocate_radiation_lw     : deallocate all local arrays used in driver_radiation_lw.
!> deallocate_radiation_sw     : deallocate all local arrays used in driver_radiation_sw.
!> deallocate_sfclayer         : deallocate all local arrays used in driver_sfclayer.
!>
!> MPAS_to_physics             :
!> driver_cloudiness           : driver for parameterization of fractional cloudiness.
!> driver_convection           : driver for parameterization of convection.
!> driver_gwdo                 : driver for parameterization of gravity wave drag over orography.
!> driver_lsm                  : driver for land-surface scheme.
!> driver_pbl                  : driver for planetary boundary layer scheme.
!> driver_radiation_sw         : driver for short wave radiation schemes.
!> driver_radiation_lw         : driver for long wave radiation schemes.
!> driver_sfclayer             : driver for surface layer scheme.
!> update_convection_step1     : updates lifetime of deep convective clouds in Kain-Fritsch scheme.
!> update_convection_step2     : updates accumulated precipitation output from convection schemes.
!> update_radiation_diagnostics: updates accumualted radiation diagnostics from radiation schemes.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * removed call to calculate atmospheric soundings for the hydrostatic dynamical core.
!>      Laura D. Fowler (2013-05-06).
!>    * removed the namelist option config_eddy_scheme and associated sourcecode.
!>    * removed the namelist option config_conv_shallow_scheme and associated sourcecode.
!>      Laura D. Fowler (birch.ucar.edu) / 2013-05-29.
!>    * added block%atm_input in calls to subroutines driver_radiation_lw amd driver_radiation_lw.
!>      Laura D. Fowler (laura@ucar.edu) / 2013-07-03.
!>    * modified sourcecode to use pools.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
!>    * renamed config_conv_deep_scheme to config_convection_scheme.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-09-18.


 contains


!==================================================================================================
 subroutine physics_driver(domain,itimestep,xtime_s)
!==================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(domain_type),intent(inout):: domain

!local pointers:
 type(mpas_pool_type),pointer::  configs,      &
                                 mesh,         &
                                 state,        &
                                 diag,         &
                                 diag_physics, &
                                 tend_physics, &
                                 atm_input,    &
                                 sfc_input

 real(kind=RKIND),pointer:: config_bucket_radt

 character(len=StrKIND),pointer:: config_bucket_update,     &
                                  config_convection_scheme, &
                                  config_gwdo_scheme,       &
                                  config_lsm_scheme,        &
                                  config_pbl_scheme,        &
                                  config_radt_lw_scheme,    &
                                  config_radt_sw_scheme,    &
                                  config_sfclayer_scheme

!local variables:
 type(block_type),pointer:: block

 integer:: time_lev
 integer:: thread

 double precision omp_get_wtime
 double precision wtime_start, wtime_end, wtime_diff
 double precision :: tarr(128), tmin, tmax
 integer threadid, omp_get_thread_num

!==================================================================================================

 call mpas_timer_start('physics driver')

 call mpas_pool_get_config(domain%configs,'config_convection_scheme',config_convection_scheme)
 call mpas_pool_get_config(domain%configs,'config_gwdo_scheme'      ,config_gwdo_scheme      )
 call mpas_pool_get_config(domain%configs,'config_lsm_scheme'       ,config_lsm_scheme       )
 call mpas_pool_get_config(domain%configs,'config_pbl_scheme'       ,config_pbl_scheme       )
 call mpas_pool_get_config(domain%configs,'config_radt_lw_scheme'   ,config_radt_lw_scheme   )
 call mpas_pool_get_config(domain%configs,'config_radt_sw_scheme'   ,config_radt_sw_scheme   )
 call mpas_pool_get_config(domain%configs,'config_sfclayer_scheme'  ,config_sfclayer_scheme  )
 call mpas_pool_get_config(domain%configs,'config_bucket_radt'      ,config_bucket_radt      )
 call mpas_pool_get_config(domain%configs,'config_bucket_update'    ,config_bucket_update    )

 if(config_convection_scheme .ne. 'off' .or. &
    config_lsm_scheme        .ne. 'off' .or. &
    config_pbl_scheme        .ne. 'off' .or. & 
    config_radt_lw_scheme    .ne. 'off' .or. &
    config_radt_sw_scheme    .ne. 'off' .or. &
    config_sfclayer_scheme   .ne. 'off') then

 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'        ,mesh        )
    call mpas_pool_get_subpool(block%structs,'state'       ,state       )
    call mpas_pool_get_subpool(block%structs,'diag'        ,diag        )
    call mpas_pool_get_subpool(block%structs,'diag_physics',diag_physics)
    call mpas_pool_get_subpool(block%structs,'atm_input'   ,atm_input   )
    call mpas_pool_get_subpool(block%structs,'sfc_input'   ,sfc_input   )
    call mpas_pool_get_subpool(block%structs,'tend_physics',tend_physics)

    !allocate arrays shared by all physics parameterizations:
    call allocate_forall_physics

    !physics prep step:
    time_lev = 1

    write(0,*) 'RSS START MPAS_to_physics'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
    do thread=1,block % nThreads
       wtime_start = omp_get_wtime()
       threadid = omp_get_thread_num()
       call MPAS_to_physics(mesh,state,time_lev,diag,diag_physics, &
                            block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
       wtime_end = omp_get_wtime()
       wtime_diff = wtime_end - wtime_start
       tarr(threadid+1) = wtime_diff
       write(0,'(A50,I3,A6,F13.7)') 'RSS MPAS_to_physics: threadid ', threadid, ' time ', wtime_diff
    end do
!$OMP END PARALLEL DO
    tmin = minval(tarr(1:block%nThreads))
    tmax = maxval(tarr(1:block%nThreads))
    write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS MPAS_to_physics: stats ', tmin, tmax, tmax-tmin

    !call to cloud scheme:
    if(l_radtlw .or. l_radtsw) then
       call allocate_cloudiness
       write(0,*) 'RSS START driver_cloudiness'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_cloudiness(diag_physics, &
                                 block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_cloudiness: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
    tmin = minval(tarr(1:block%nThreads))
    tmax = maxval(tarr(1:block%nThreads))
    write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_cloudiness: stats ', tmin, tmax, tmax-tmin
    endif

    !call to short wave radiation scheme:
    if(l_radtsw) then
       time_lev = 1
       call allocate_radiation_sw(xtime_s)
       write(0,*) 'RSS START driver_radiation_sw'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_radiation_sw(itimestep,block%configs,mesh,state,time_lev,diag_physics, &
                                   atm_input,sfc_input,tend_physics,xtime_s, &
                                   block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_radiation_sw: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_radiation_sw: stats ', tmin, tmax, tmax-tmin
    endif

    !call to long wave radiation scheme:
    if(l_radtlw) then
       time_lev = 1
       call allocate_radiation_lw(xtime_s)
       write(0,*) 'RSS START driver_radiation_lw'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_radiation_lw(xtime_s,block%configs,mesh,state,time_lev,diag_physics, &
                                   atm_input,sfc_input,tend_physics, &
                                   block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_radiation_lw: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_radiation_lw: stats ', tmin, tmax, tmax-tmin
    endif

    !call to accumulate long- and short-wave diagnostics if needed:
     if(config_bucket_update /= 'none' .and. config_bucket_radt .gt. 0._RKIND) then
        write(0,*) 'RSS START update_radiation_diagnostics'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
        do thread=1,block % nThreads
           wtime_start = omp_get_wtime()
           threadid = omp_get_thread_num()
           call update_radiation_diagnostics(block%configs,mesh,diag_physics, &
                                            block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS update_radiation_diagnostics: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS update_radiation_diagnostics: stats ', tmin, tmax, tmax-tmin
     endif

    !deallocate all radiation arrays:
    if(config_radt_sw_scheme.ne.'off' .or. config_radt_lw_scheme.ne.'off') &
       call deallocate_cloudiness
    if(config_radt_sw_scheme.ne.'off') call deallocate_radiation_sw
    if(config_radt_lw_scheme.ne.'off') call deallocate_radiation_lw

    !call to surface-layer scheme:
    if(config_sfclayer_scheme .ne. 'off') then
       call allocate_sfclayer
       write(0,*) 'RSS START driver_sfclayer'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_sfclayer(mesh,diag_physics,sfc_input, &
                               block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_sfclayer: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_sfclayer: stats ', tmin, tmax, tmax-tmin
       call deallocate_sfclayer
    endif

    !call to land-surface scheme:
    if(config_lsm_scheme .ne. 'off') then
       call allocate_lsm
       write(0,*) 'RSS START driver_lsm'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_lsm(itimestep,block%configs,mesh,diag_physics,sfc_input, &
                          block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_lsm: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_lsm: stats ', tmin, tmax, tmax-tmin
       call deallocate_lsm
    endif

    !call to pbl schemes:
    if(config_pbl_scheme .ne. 'off' .and. config_sfclayer_scheme .ne. 'off') then
       call allocate_pbl
       write(0,*) 'RSS START driver_pbl'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_pbl(sfc_input,diag_physics,tend_physics, &
                          block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_pbl: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_pbl: stats ', tmin, tmax, tmax-tmin
       call deallocate_pbl
    endif

    !call to gravity wave drag over orography scheme:
    if(config_gwdo_scheme .ne. 'off') then
       call allocate_gwdo
       write(0,*) 'RSS START driver_gwdo'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_gwdo(itimestep,mesh,sfc_input,diag_physics,tend_physics, &
                           block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_gwdo: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_gwdo: stats ', tmin, tmax, tmax-tmin
       call deallocate_gwdo
    endif

    !call to convection scheme:
    write(0,*) 'RSS START update_convection_step1'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
    do thread=1,block % nThreads
       wtime_start = omp_get_wtime()
       threadid = omp_get_thread_num()
       call update_convection_step1(mesh,diag_physics,tend_physics, &
                            block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
       wtime_end = omp_get_wtime()
       wtime_diff = wtime_end - wtime_start
       tarr(threadid+1) = wtime_diff
       write(0,'(A50,I3,A6,F13.7)') 'RSS update_convection_step1: threadid ', threadid, ' time ', wtime_diff
    end do
!$OMP END PARALLEL DO
    tmin = minval(tarr(1:block%nThreads))
    tmax = maxval(tarr(1:block%nThreads))
    write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS update_convection_step1: stats ', tmin, tmax, tmax-tmin
    if(l_conv) then
       call allocate_convection
       write(0,*) 'RSS START driver_convection'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
       do thread=1,block % nThreads
          wtime_start = omp_get_wtime()
          threadid = omp_get_thread_num()
          call driver_convection(itimestep,mesh,sfc_input,diag_physics,tend_physics, &
                                 block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
          wtime_end = omp_get_wtime()
          wtime_diff = wtime_end - wtime_start
          tarr(threadid+1) = wtime_diff
          write(0,'(A50,I3,A6,F13.7)') 'RSS driver_convection: threadid ', threadid, ' time ', wtime_diff
       end do
!$OMP END PARALLEL DO
       tmin = minval(tarr(1:block%nThreads))
       tmax = maxval(tarr(1:block%nThreads))
       write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS driver_convection: stats ', tmin, tmax, tmax-tmin
       call deallocate_convection
    endif
    !update diagnostics:
    write(0,*) 'RSS START update_convection_step2'
!$OMP PARALLEL DO PRIVATE(wtime_end, wtime_start, threadid)
    do thread=1,block % nThreads
       wtime_start = omp_get_wtime()
       threadid = omp_get_thread_num()
       call update_convection_step2(block%configs,mesh,diag_physics, &
                            block % cellSolveThreadStart(thread), block % cellSolveThreadEnd(thread))
       wtime_end = omp_get_wtime()
       wtime_diff = wtime_end - wtime_start
       tarr(threadid+1) = wtime_diff
       write(0,'(A50,I3,A6,F13.7)') 'RSS update_convection_step2: threadid ', threadid, ' time ', wtime_diff
    end do
!$OMP END PARALLEL DO
    tmin = minval(tarr(1:block%nThreads))
    tmax = maxval(tarr(1:block%nThreads))
    write(0,'(A50,1x,F13.7,1x,F13.7,1x,F13.7)') 'RSS update_convection_step2: stats ', tmin, tmax, tmax-tmin

    !deallocate arrays shared by all physics parameterizations:
    call deallocate_forall_physics

    block => block % next
 end do 

 endif
! write(0,*)
! write(0,*) '--- end physics_driver:'
! write(0,*)

 call mpas_timer_stop('physics driver')

 end subroutine physics_driver

!==================================================================================================
 end module mpas_atmphys_driver
!==================================================================================================
