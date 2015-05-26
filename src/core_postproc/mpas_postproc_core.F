! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module postproc_core

   use mpas_framework
   use mpas_timekeeping

   type (MPAS_Clock_type), pointer :: clock

   contains

   !***********************************************************************
   !
   !  function postproc_core_init
   !
   !> \brief   MPAS Core Initialization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs the necessary initialization of a core. This can
   !>   involve things like setting up coeffiecients for remapping and advection,
   !>   but more generally can include anything the core needs to initialize.
   !
   !-----------------------------------------------------------------------
   function postproc_core_init(domain, startTimeStamp) result(iErr)!{{{
   
      use mpas_configure
      use mpas_derived_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      type (MPAS_Time_Type) :: startTime
      integer :: iErr


      iErr = 0

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, iErr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

   end function postproc_core_init!}}}


   !***********************************************************************
   !
   !  function postproc_core_run
   !
   !> \brief   MPAS Core Run Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs a run of the MPAS core. This can involve time
   !>  stepping if the core requires time stepping, but more generally includes
   !>  anything a core would typically do after it was initialized.
   !
   !-----------------------------------------------------------------------
   function postproc_core_run(domain) result(iErr)!{{{
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
      use mpas_stream_manager
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: iErr

      type (mpas_pool_type), pointer :: pool
      type (mpas_pool_iterator_type) :: itr
      type (MPAS_time_type) :: currTime
      character (len=StrKIND) :: timeStamp

      type (mpas_pool_type), pointer :: mesh
      integer, pointer :: nCells
      real (kind=RKIND), dimension(:), pointer :: xCell

      do while (.not. mpas_is_clock_stop_time(clock))
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         

         write(0,*) 'Processing '//trim(timeStamp)

         call mpas_stream_mgr_read(domain % streamManager, ierr=ierr)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
         call mpas_pool_get_dimension(mesh, 'nCells', nCells)
         call mpas_pool_get_array(mesh, 'xCell', xCell)

         write(0,*) 'nCells = ', nCells

         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)

         call mpas_advance_clock(clock)
      end do

      iErr = 0

   end function postproc_core_run!}}}

   !***********************************************************************
   !
   !  function postproc_core_finalize
   !
   !> \brief   MPAS Core Finalization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function finalizes the MPAS core. It should at a minimum destroy the
   !>  simulation clock, but can perform any functions necessary to clean up the
   !>  MPAS core.
   !
   !-----------------------------------------------------------------------
   function postproc_core_finalize(domain) result(iErr)!{{{
   
      use mpas_derived_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: iErr
 

      iErr = 0

      call mpas_destroy_clock(clock, iErr)

   end function postproc_core_finalize!}}}

end module postproc_core
