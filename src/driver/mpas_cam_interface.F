module mpas_cam_interface

   use mpas_grid_types
   use mpas_configure
   use mpas_io_input
   use mpas_io_output
   use mpas_dmpar
   use mpas_constants
   use mpas_core


#ifdef CAM_DYCORE
   type (dm_info), pointer :: dminfo
   type (domain_type), pointer :: domain
   real (kind=RKIND) :: dt_dynamics, dt_physics, p0
   integer :: itimestep, n_subcycle_steps
   type (io_output_object), save :: output_obj
   logical :: restart_run
   integer :: output_frame
#endif


   contains


#ifdef CAM_DYCORE
   subroutine mpas_init1(mpi_comm, phys_dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_INIT1
   !
   ! Initializes the MPAS software infrastructure, reads run-time configuration
   !    information, reads grid information from an MPAS grid file, and allocates
   !    storage for fields to be provided by CAM through either 
   !    cam_inidat_to_mpas() or cam_restart_to_mpas().
   !
   ! Input: mpi_comm - an MPI communicator supplied by CAM to be used by MPAS
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      integer, intent(in) :: mpi_comm
      real (kind=RKIND), intent(in) :: phys_dt

      integer :: iCell, iEdge, j
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), pointer :: east, north

      write(0,*) 'Called MPAS_INIT1'

      call mpas_framework_init(dminfo, domain, mpi_comm)

      itimestep = 1
      restart_run = .false.    ! If a restart run, this will be set to true in cam_restart_to_mpas

      !
      ! Set physics timestep and verify that it is evenly divided by dynamics timestep
      !
      dt_physics = phys_dt
      dt_dynamics = dt_physics / ceiling(dt_physics / config_dt)
      n_subcycle_steps = ceiling(dt_physics / config_dt)



      !
      ! Read initial conditions
      !
      call mpas_input_state_for_domain(domain)

      latCell => domain % blocklist % mesh % latCell % array
      lonCell => domain % blocklist % mesh % lonCell % array
      east => domain % blocklist % mesh % east % array
      north => domain % blocklist % mesh % north % array


      !
      ! Compute unit vectors in east and north directions for each cell
      !
      do iCell = 1,domain % blocklist % mesh % nCellsSolve

         east(1,iCell) = -sin(lonCell(iCell))
         east(2,iCell) =  cos(lonCell(iCell))
         east(3,iCell) =  0.0
         call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

         north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
         north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
         north(3,iCell) =  cos(latCell(iCell))
         call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

      end do

      p0 = 1.0E5    ! This should be consistent with the value used in MPAS to
                    !    originally compute theta from temperature
   
   end subroutine mpas_init1
   
   
   subroutine mpas_get_pref_profile(pref_edge, pref_mid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_GET_PREF_PROFILE
   !
   ! Returns reference profiles for pressure at layer interfaces and layer
   !    midpoints. Pressures should be in ascending order, i.e., from the top of
   !    the atmosphere downward.
   !
   ! Output: pref_edge - reference pressure at layer interfaces (vertical edges)
   !         pref_mid - reference pressure at layer midpoints
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: pref_edge
      real (kind=RKIND), dimension(:), intent(out) :: pref_mid


      integer :: k, nVertLevels, nCellsSolve
      real (kind=RKIND), dimension(:), pointer :: temp
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp

      nVertLevels = domain % blocklist % mesh % nVertLevels
      nCellsSolve = domain % blocklist % mesh % nCellsSolve

      fzm => domain % blocklist % mesh % fzm % array
      fzp => domain % blocklist % mesh % fzp % array

      allocate(temp(nVertLevels))

      !
      ! Take layer midpoint pressures directly from pressure_base
      !
      do k=1,nVertLevels
         temp(nVertLevels-k+1) = maxval(domain % blocklist % diag % pressure_base % array(k,1:nCellsSolve))
      end do

      call mpas_dmpar_max_real_array(domain % dminfo, nVertLevels, temp, pref_mid)

      deallocate(temp)

      !
      ! Interpolate interface pressure for interior edges, and extrapolate for
      !    top and bottom pressure
      !
      do k=2,nVertLevels
         pref_edge(k) = fzm(k)*pref_mid(k) + fzp(k)*pref_mid(k-1)
      end do
      pref_edge(1) = pref_mid(1) + pref_mid(1) - pref_edge(2)
      pref_edge(nVertLevels+1) = pref_mid(nVertLevels) + pref_mid(nVertLevels) - pref_edge(nVertLevels)

   end subroutine mpas_get_pref_profile
   
   
   subroutine mpas_init2()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_INIT2
   !
   ! Calls the core-specific (most likely the hydrostatic atmosphere core) 
   !    initialization routine mpas_core_init() after initial fields have been provided
   !    through a call to either cam_inidat_to_mpas() or cam_restart_to_mpas().
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      integer :: i, iCell, k, iScalar
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      type (block_type), pointer :: block

      integer :: index_qv
      character (len=StrKIND) :: initialTimeStamp
   
      write(0,*) 'Called MPAS_INIT2' 


      ! Compute diagnostic fields needed in solve loop, and initialize
      !    simulation time to 0 for all blocks
      block => domain % blocklist
      do while (associated(block))
         call mpas_core_init(domain, initialTimeStamp)
         block => block % next
      end do

      output_frame = 1
      if(config_frames_per_outfile > 0) then
         call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(initialTimeStamp))
      else
         call mpas_output_state_init(output_obj, domain, "OUTPUT")
      end if
      call atm_write_output_frame(output_obj, output_frame, domain)

   end subroutine mpas_init2
   
   
   subroutine mpas_to_cam(Numcols, Plev, Pcnst, Psd, Phis, Pint, Pmid, Zint, Zmid, &
                          T, Ux, Uy, Omega, Tracer)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_TO_CAM
   !
   ! Input: Numcols - number of columns/cells
   !        Plev - number of vertical levels
   !        Pcnst - number of tracers
   !
   ! Output: Psd - dry surface pressure (Pa)
   !         Phis - surface geopotential (m^2/s^2)
   !         Pint - dry pressure at vertical layer interfaces
   !         Pmid - dry pressure at vertical layer mid-points
   !         Zint - geopotential height at vertical layer interfaces
   !         Zmid - geopotential height at vertical layer mid-points
   !         T - temperature (K)
   !         Ux - longitudinal velocity at cell centers (m/s)
   !         Uy - latitudinal velocity at cell centers (m/s)
   !         Omega - omega (Pa/s)
   !         Tracer - tracer mixing ratios (kg/kg)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      integer, intent(in) :: Numcols
      integer, intent(in) :: Plev
      integer, intent(in) :: Pcnst
      real (kind=RKIND), dimension(Numcols), intent(out) :: Psd
      real (kind=RKIND), dimension(Numcols), intent(out) :: Phis
      real (kind=RKIND), dimension(Numcols,Plev+1), intent(out) :: Pint
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Pmid
      real (kind=RKIND), dimension(Numcols,Plev+1), intent(out) :: Zint
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Zmid
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: T
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Ux
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Uy
      real (kind=RKIND), dimension(Numcols,Plev), intent(out) :: Omega
      real (kind=RKIND), dimension(Numcols,Plev,Pcnst), intent(out) :: Tracer

      integer :: iCell, k, iScalar
      integer :: index_qv
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:,:), pointer :: theta, pressure_base, pressure_p, ww, uReconstZonal, uReconstMeridional, east, north
      real (kind=RKIND), dimension(:,:), pointer :: w, rho_p, rho_base, zz
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      type (block_type), pointer :: block

      write(0,*) 'Called MPAS_TO_CAM'


      block => domain % blocklist

      index_qv = block % state % time_levs(1) % state % index_qv


      !
      !  Perform basic sanity check on expected and available field dimensions
      !
      if (Numcols /= block % mesh % nCellsSolve) then
         write(0,*) 'Error: mismatch between Numcols and nCellsSolve: ', Numcols, block % mesh % nCellsSolve
         return
      end if
      if (Plev /= block % mesh % nVertLevels) then
         write(0,*) 'Error: mismatch between Plev and nVertLevels: ', Plev, block % mesh % nVertLevels
         return
      end if
      if (Pcnst /= block % state % time_levs(1) % state % num_scalars) then
         write(0,*) 'Error: mismatch between Pcnst and num_scalars: ', Pcnst, block % state % time_levs(1) % state % num_scalars
         return
      end if

      uReconstZonal => block % diag % uReconstructZonal % array
      uReconstMeridional => block % diag % uReconstructMeridional % array
      scalars => block % state % time_levs(1) % state % scalars % array
      theta => block % diag % theta % array
      pressure_base => block % diag % pressure_base % array
      pressure_p => block % diag % pressure_p % array
      w => block % state % time_levs(1) % state % w % array
      zz => block % mesh % zz % array
      rho_p => block % diag % rho_p % array
      rho_base => block % diag % rho_base % array

      !
      ! Fill in CAM arrays from block % state % time_levs(1) % state arrays
      !
      do iCell=1,block % mesh % nCellsSolve
         Psd(iCell) = block % diag % surface_pressure % array(iCell)
         Phis(iCell) = block % mesh % zgrid % array(1,iCell) * 9.81
         do k=1,block % mesh % nVertLevels
            theta(k,iCell) = block % state % time_levs(1) % state % theta_m % array(k,iCell) &
                             / (1.0 + 1.61 * scalars(index_qv,k,iCell)) 
            T(iCell,k) = theta(k,iCell) * ((pressure_base(k,iCell)+pressure_p(k,iCell)) / p0) ** (rgas/cp)
            Pmid(iCell,k) = pressure_base(k,iCell)+pressure_p(k,iCell)
            if (k == 1) then
               Pint(iCell,k) = Psd(iCell)
            else
               Pint(iCell,k) = block % mesh % fzm % array(k) * Pmid(iCell,k-1) + &
                               block % mesh % fzp % array(k) * Pmid(iCell,k)
            end if
            Zint(iCell,k) = block % mesh % zgrid % array(k,iCell)
            Zmid(iCell,k) = 0.5*(block % mesh % zgrid % array(k,iCell) + block % mesh % zgrid % array(k+1,iCell))

!            Omega(iCell,k) = 0.5*(ww(k,iCell) + ww(k+1,iCell))
            Omega(iCell,k) = -1.0 * gravity * (rho_p(k,iCell)+rho_base(k,iCell)) * zz(k,iCell) * (w(k,iCell)+w(k+1,iCell))
            !
            ! NOTE: Eventually, we need to ensure that the moisture variables we return to CAM
            !       are what CAM expects, rather than simply what we have in our scalars array
            !
            do iScalar=1,block % state % time_levs(1) % state % num_scalars
               Tracer(iCell,k,iScalar) = scalars(iScalar,k,iCell) / (1.0 + scalars(index_qv,k,iCell))
            end do
         end do
         Zint(iCell,Plev+1) = block % mesh % zgrid % array(Plev+1,iCell)
         Pint(iCell,Plev+1) = Pmid(iCell,Plev) + (Pmid(iCell,Plev) - Pint(iCell,Plev))
      end do

      do iCell=1,block % mesh % nCellsSolve
         do k=1,block % mesh % nVertLevels
            Ux(iCell,k) =   uReconstZonal(k,iCell)
            Uy(iCell,k) =   uReconstMeridional(k,iCell)
         end do
      end do 

write(0,*) '<- MPAS DEBUGGING: MIN/MAX uReconstructZonal=', minval(Ux), maxval(Ux)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX uReconstructMeridional=', minval(Uy), maxval(Uy)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Psd=', minval(Psd), maxval(Psd)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Phis=', minval(Phis), maxval(Phis)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Pint=', minval(Pint), maxval(Pint)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Pmid=', minval(Pmid), maxval(Pmid)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Zint=', minval(Zint), maxval(Zint)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Zmid=', minval(Zmid), maxval(Zmid)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX T=', minval(T), maxval(T)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Omega=', minval(Omega), maxval(Omega)
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Tracer1=', minval(Tracer(:,:,1)), maxval(Tracer(:,:,1))
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Tracer2=', minval(Tracer(:,:,2)), maxval(Tracer(:,:,2))
write(0,*) '<- MPAS DEBUGGING: MIN/MAX Tracer3=', minval(Tracer(:,:,3)), maxval(Tracer(:,:,3))

   end subroutine mpas_to_cam
   
   
   subroutine cam_to_mpas(Numcols, Plev, Pcnst, T_tend, Ux_tend, Uy_tend, Tracer)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CAM_TO_MPAS
   !
   ! Input: Numcols - number of columns/cells
   !        Plev - number of vertical levels
   !        Pcnst - number of tracers
   !        T_tend - temperature tendency (K/s)
   !        Ux_tend - longitudinal velocity tendency at cell centers (m/s^2)
   !        Uy_tend - latitudinal velocity tendency at cell centers (m/s^2)
   !        Tracer - tracer mixing ratios (kg/kg)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      integer, intent(in) :: Numcols
      integer, intent(in) :: Plev
      integer, intent(in) :: Pcnst
      real (kind=RKIND), dimension(Numcols,Plev), intent(in) :: T_tend
      real (kind=RKIND), dimension(Numcols,Plev), intent(in) :: Ux_tend
      real (kind=RKIND), dimension(Numcols,Plev), intent(in) :: Uy_tend
      real (kind=RKIND), dimension(Numcols,Plev,Pcnst), intent(in) :: Tracer

      integer :: iCell, iEdge, iScalar, k, j
      integer :: index_qv
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), dimension(:,:), pointer :: pressure_base, pressure_p
      real (kind=RKIND), dimension(:,:), pointer :: east, north, edge_normal, pressure
      real (kind=RKIND), dimension(:,:), pointer :: theta_tend, u_tend
      real (kind=RKIND), dimension(:,:), pointer :: zonal_tend, meridional_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars, scalars_tend
      type (block_type), pointer :: block

      write(0,*) 'Called CAM_TO_MPAS'

      block => domain % blocklist

      index_qv = block % state % time_levs(1) % state % index_qv


      !
      !  Perform basic sanity check on expected and available field dimensions
      !
      if (Numcols /= block % mesh % nCellsSolve) then
         write(0,*) 'Error: mismatch between Numcols and nCellsSolve: ', Numcols, block % mesh % nCellsSolve
         return
      end if
      if (Plev /= block % mesh % nVertLevels) then
         write(0,*) 'Error: mismatch between Plev and nVertLevels: ', Plev, block % mesh % nVertLevels
         return
      end if
      if (Pcnst /= block % state % time_levs(1) % state % num_scalars) then
         write(0,*) 'Error: mismatch between Pcnst and num_scalars: ', Pcnst, block % state % time_levs(1) % state % num_scalars
         return
      end if

      nEdgesOnCell => block % mesh % nEdgesOnCell % array
      edgesOnCell => block % mesh % edgesOnCell % array
      east => block % mesh % east % array
      north => block % mesh % north % array
      edge_normal => block % mesh % edgeNormalVectors % array
      scalars => block % state % time_levs(1) % state % scalars % array
      pressure_base => block % diag % pressure_base % array
      pressure_p => block % diag % pressure_p % array
      scalars_tend => block % tend_physics % scalars % array
      theta_tend => block % tend_physics % theta % array
      u_tend => block % tend_physics % u % array
      zonal_tend => block % tend_physics % ux % array
      meridional_tend => block % tend_physics % uy % array


      !
      ! Fill in MPAS tendency arrays from arguments
      !
      do iCell=1,block % mesh % nCellsSolve
         do k=1,block % mesh % nVertLevels
            zonal_tend(k,iCell) = Ux_tend(iCell,k)
            meridional_tend(k,iCell) = Uy_tend(iCell,k)
         end do
      end do

      call mpas_dmpar_exch_halo_field(block % tend_physics % ux)
      call mpas_dmpar_exch_halo_field(block % tend_physics % uy)

      u_tend(:,:) = 0.0
      do iCell=1,block % mesh % nCells
         do j=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(j,iCell)
            do k=1,block % mesh % nVertLevels
               u_tend(k,iEdge) = u_tend(k,iEdge) + 0.5 * zonal_tend(k,iCell) * (edge_normal(1,iEdge) * east(1,iCell) + &
                                                                                edge_normal(2,iEdge) * east(2,iCell) + &
                                                                                edge_normal(3,iEdge) * east(3,iCell))  &
                                                 + 0.5 * meridional_tend(k,iCell) * (edge_normal(1,iEdge) * north(1,iCell) + &
                                                                                     edge_normal(2,iEdge) * north(2,iCell) + &
                                                                                     edge_normal(3,iEdge) * north(3,iCell))
            end do
         end do

      end do

      do iCell=1,block % mesh % nCellsSolve
         do k=1,block % mesh % nVertLevels
            theta_tend(k,iCell) = T_tend(iCell,k) * (p0 / (pressure_base(k,iCell)+pressure_p(k,iCell))) ** (rgas/cp)
            
            ! Couple theta tendency with h valid at time when tendency was computed
!            theta_tend(k,iCell) = theta_tend(k,iCell) * h(k,iCell) 

            !
            ! NOTE: Once we begin to use more than just qv in MPAS, we need to make sure to
            !       properly assign variables from Tracer array provided from CAM physics
            !
            do iScalar=1,block % state % time_levs(1) % state % num_scalars
               scalars_tend(iScalar,k,iCell) = (Tracer(iCell,k,iScalar) * (1.0 + Tracer(iCell,k,index_qv)) &
                                              - scalars(iScalar,k,iCell)&
                                               ) / dt_physics

               ! Couple scalar tendencies with h valid at time when tendencies were computed
!               scalars_tend(iScalar,k,iCell) = scalars_tend(iScalar,k,iCell) * h(k,iCell) 
            end do
         end do
      end do


      !
      ! Do a ghost-cell update on tendency arrays
      !
      call mpas_dmpar_exch_halo_field(block % tend_physics % theta)
      call mpas_dmpar_exch_halo_field(block % tend_physics % u)
      call mpas_dmpar_exch_halo_field(block % tend_physics % scalars)

write(0,*) '-> MPAS DEBUGGING: MIN/MAX u_tend=', minval(u_tend(:,1:block % mesh % nEdges)),maxval(u_tend(:,1:block % mesh % nEdges))
write(0,*) '-> MPAS DEBUGGING: MIN/MAX theta_tend=', minval(theta_tend(:,1:block % mesh % nCells)),maxval(theta_tend(:,1:block % mesh % nCells))
write(0,*) '-> MPAS DEBUGGING: MIN/MAX qv_tend=', minval(scalars_tend(1,:,1:block % mesh % nCells)),maxval(scalars_tend(1,:,1:block % mesh % nCells))
write(0,*) '-> MPAS DEBUGGING: MIN/MAX qc_tend=', minval(scalars_tend(2,:,1:block % mesh % nCells)),maxval(scalars_tend(2,:,1:block % mesh % nCells))
write(0,*) '-> MPAS DEBUGGING: MIN/MAX qr_tend=', minval(scalars_tend(3,:,1:block % mesh % nCells)),maxval(scalars_tend(3,:,1:block % mesh % nCells))
   
   end subroutine cam_to_mpas
   
   
   subroutine mpas_dyn_run()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_DYN_RUN
   !
   ! Advance MPAS dynamics by one timestep.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr
      
      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: itimestep
      integer :: ierr
      integer :: idynstep

      write(0,*) 'Called MPAS_DYN_RUN'

      
      ! Eventually, dt should be domain specific
      dt = config_dt
      
      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt_dynamics in time by timestep(...)
write(0,*) 'MGD running dynamics for n_subcycle_steps=',n_subcycle_steps
      do idynstep=1,n_subcycle_steps
      
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(0,*) 'Begin timestep ', trim(timeStamp)
         
         ! Input external updates (i.e. surface)
         if (mpas_is_alarm_ringing(clock, sfcAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, sfcAlarmID, ierr=ierr)
            
            call mpas_read_and_distribute_fields(sfc_update_obj)
            sfc_update_obj % time = sfc_update_obj % time + 1
         end if
         
         call mpas_timer_start("time integration")
         call atm_do_timestep(domain, dt_dynamics, itimestep)
         call mpas_timer_stop("time integration")   
         
         ! Move time level 2 fields back into time level 1 for next time step
         call mpas_shift_time_levels_state(domain % blocklist % state)
         

         ! Advance clock before writing output
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         
         !TODO: MPAS_getClockRingingAlarms is probably faster than multiple MPAS_isAlarmRinging...
         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=ierr)
            ! output_frame will always be > 1 here unless it was reset after the maximum number of frames per outfile was reached
            if(output_frame == 1) then 
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if
            call atm_write_output_frame(output_obj, output_frame, domain)
         end if
         
         if (mpas_is_alarm_ringing(clock, restartAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=ierr)
            
            block_ptr => domain % blocklist
            do while (associated(block_ptr))
               call atm_compute_restart_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % diag, block_ptr % mesh)
               block_ptr => block_ptr % next
            end do

            ! Write one restart time per file
            call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
            call mpas_output_state_for_domain(restart_obj, domain, 1)
            call mpas_output_state_finalize(restart_obj, domain % dminfo)
         end if

      end do

   end subroutine mpas_dyn_run
   
   
   subroutine mpas_final()
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE MPAS_FINAL
   !
   ! Deallocates storage allocated by MPAS and shuts down the MPAS software
   !    infrastructure.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      implicit none
   
      write(0,*) 'Called MPAS_FINAL'

      call mpas_output_state_finalize(output_obj, domain % dminfo)

      call mpas_core_finalize(domain)

      call mpas_framework_finalize(dminfo, domain)

   
   end subroutine mpas_final


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE ROTATE_ABOUT_VECTOR
   !
   ! Rotates the point (x,y,z) through an angle theta about the vector
   !   originating at (a, b, c) and having direction (u, v, w).
   !
   ! Reference: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine rotate_about_vector(x, y, z, theta, a, b, c, u, v, w, xp, yp, zp)
   
      implicit none
   
      real (kind=RKIND), intent(in) :: x, y, z, theta, a, b, c, u, v, w
      real (kind=RKIND), intent(out) :: xp, yp, zp
   
      real (kind=RKIND) :: vw2, uw2, uv2
      real (kind=RKIND) :: m
   
      vw2 = v**2.0 + w**2.0
      uw2 = u**2.0 + w**2.0
      uv2 = u**2.0 + v**2.0
      m = sqrt(u**2.0 + v**2.0 + w**2.0)
   
      xp = (a*vw2 + u*(-b*v-c*w+u*x+v*y+w*z) + ((x-a)*vw2+u*(b*v+c*w-v*y-w*z))*cos(theta) + m*(-c*v+b*w-w*y+v*z)*sin(theta))/m**2.0
      yp = (b*uw2 + v*(-a*u-c*w+u*x+v*y+w*z) + ((y-b)*uw2+v*(a*u+c*w-u*x-w*z))*cos(theta) + m*( c*u-a*w+w*x-u*z)*sin(theta))/m**2.0
      zp = (c*uv2 + w*(-a*u-b*v+u*x+v*y+w*z) + ((z-c)*uv2+w*(a*u+b*v-u*x-v*y))*cos(theta) + m*(-b*u+a*v-v*x+u*y)*sin(theta))/m**2.0
   
   end subroutine rotate_about_vector


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE R3_CROSS
   !
   ! Computes the cross product of (ax, ay, az) and (bx, by, bz).
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine r3_cross(ax, ay, az, bx, by, bz, cx, cy, cz)

      implicit none

      real (kind=RKIND), intent(in) :: ax, ay, az
      real (kind=RKIND), intent(in) :: bx, by, bz
      real (kind=RKIND), intent(out) :: cx, cy, cz

      cx = ay * bz - az * by
      cy = az * bx - ax * bz
      cz = ax * by - ay * bx

   end subroutine r3_cross 


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE R3_NORMALIZE
   !
   ! Normalizes the vector (ax, ay, az)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine r3_normalize(ax, ay, az)

      implicit none

      real (kind=RKIND), intent(inout) :: ax, ay, az

      real (kind=RKIND) :: mi

      mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
      ax = ax * mi
      ay = ay * mi
      az = az * mi

   end subroutine r3_normalize 
#endif

end module mpas_cam_interface
