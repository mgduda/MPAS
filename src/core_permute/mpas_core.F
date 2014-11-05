! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_core

   use mpas_framework
   use mpas_timekeeping

   type (MPAS_Clock_type), pointer :: clock


   contains


   subroutine mpas_core_init(domain, stream_manager, startTimeStamp)
   
      use mpas_configure
      use mpas_grid_types
      use mpas_stream_manager
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      character(len=*), intent(out) :: startTimeStamp
   
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool

      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime
      integer :: ierr


      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock


      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_array(statePool, 'xtime', xtime, 1)

         xtime = startTimeStamp
         block => block % next
      end do

      call MPAS_stream_mgr_read(stream_manager, ierr=ierr)

   end subroutine mpas_core_init


   subroutine mpas_core_run(domain, stream_manager)
   
      use mpas_grid_types
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager

      type (block_type), pointer :: block_ptr
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool

      integer :: i
      integer, pointer :: nCells, nEdges, nVertices
      integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve
      integer :: cellStart, edgeStart, vertexStart
      integer, dimension(:), pointer :: indexToCellID, indexToEdgeID, indexToVertexID
      type (field1DInteger), pointer :: indexToCellIDField, indexToEdgeIDField, indexToVertexIDField
      integer :: ierr


      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
      call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)

      call mpas_pool_get_field(meshPool, 'indexToCellID', indexToCellIDField)
      call mpas_pool_get_field(meshPool, 'indexToEdgeID', indexToEdgeIDField)
      call mpas_pool_get_field(meshPool, 'indexToVertexID', indexToVertexIDField)

      write(0,*) 'nCellsSolve = ', nCellsSolve
      write(0,*) 'nEdgesSolve = ', nEdgesSolve
      write(0,*) 'nVerticesSolve = ', nVerticesSolve

      call mpas_dmpar_get_start(domain % dminfo, nCellsSolve, cellStart, write_partition=.true.)
      call mpas_dmpar_get_start(domain % dminfo, nEdgesSolve, edgeStart)
      call mpas_dmpar_get_start(domain % dminfo, nVerticesSolve, vertexStart)

      write(0,*) 'cellStart=', cellStart
      write(0,*) 'edgeStart=', edgeStart
      write(0,*) 'vertexStart=', vertexStart

      do i=1,nCellsSolve
         indexToCellID(i) = cellStart + i - 1
      end do
      call mpas_dmpar_exch_halo_field(indexToCellIDField)

      do i=1,nEdgesSolve
         indexToEdgeID(i) = edgeStart + i - 1
      end do
      call mpas_dmpar_exch_halo_field(indexToEdgeIDField)

      do i=1,nVerticesSolve
         indexToVertexID(i) = vertexStart + i - 1
      end do
      call mpas_dmpar_exch_halo_field(indexToVertexIDField)

   end subroutine mpas_core_run

   
   subroutine mpas_core_finalize(domain, stream_manager)
   
      use mpas_grid_types
      use mpas_stream_manager
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      integer :: ierr

      call MPAS_stream_mgr_write(stream_manager, ierr=ierr)
 
      call mpas_destroy_clock(clock, ierr)

   end subroutine mpas_core_finalize


   !***********************************************************************
   !
   !  routine mpas_core_setup_packages
   !
   !> \brief   Pacakge setup routine
   !> \author  Doug Jacobsen
   !> \date    September 2011
   !> \details 
   !>  This routine is intended to correctly configure the packages for this MPAS
   !>   core. It can use any Fortran logic to properly configure packages, and it
   !>   can also make use of any namelist options. All variables in the model are
   !>   *not* allocated until after this routine is called.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_setup_packages(configPool, packagePool, ierr)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool
      integer, intent(out) :: ierr

      ierr = 0

   end subroutine mpas_core_setup_packages!}}}


   !***********************************************************************
   !
   !  routine mpas_core_setup_clock
   !
   !> \brief   Pacakge setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_setup_clock(core_clock, configs, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (MPAS_time_type) :: startTime, stopTime
      type (MPAS_timeInterval_type) :: timeStep
      character (len=StrKIND) :: config_start_time, config_stop_time 
      real (kind=RKIND) :: config_dt


      ierr = 0

      config_start_time = '0000-01-01_00:00:00'
      config_stop_time = '0000-01-01_00:00:00'
      config_dt = 86400

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=ierr)
      call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
      call mpas_set_timeInterval(timeStep, dt=config_dt, ierr=ierr)

      call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=ierr)

   end subroutine mpas_core_setup_clock


   !***********************************************************************
   !
   !  routine mpas_core_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_get_mesh_stream(configs, stream, ierr)

      implicit none

      type (mpas_pool_type), intent(in) :: configs
      character(len=*), intent(out) :: stream
      integer, intent(out) :: ierr


      ierr = 0

      write(stream,'(a)') 'input'

   end subroutine mpas_core_get_mesh_stream

end module mpas_core
