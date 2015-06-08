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

   real (kind=RKIND), parameter :: Re = 6371229.0
   real (kind=RKIND), parameter :: Pi = 3.14159265359
   real (kind=RKIND), parameter :: rad2deg = 57.2957795131

   integer, parameter :: topo_x = 43200
   integer, parameter :: topo_y = 21600
   real (kind=RKIND), parameter :: pts_per_degree = real(topo_x,RKIND) / 360.0
   real (kind=RKIND), parameter :: start_lat = -90.0
   real (kind=RKIND), parameter :: start_lon = -180.0

   real (kind=RKIND), dimension(:,:), pointer :: topo
   real (kind=RKIND), dimension(:,:), pointer :: box
   real (kind=RKIND) :: box_mean
   integer :: nx, ny

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

      type (mpas_pool_type), pointer :: model, mesh, state
      integer :: iCell, i
      real (kind=RKIND) :: dc
      integer, pointer :: nCells
      character(len=StrKIND), pointer :: xtime
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, dcEdge
      real (kind=RKIND), dimension(:), pointer :: var2d, con, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4

      do while (.not. mpas_is_clock_stop_time(clock))
         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         

         write(0,*) 'Processing '//trim(timeStamp)

         call mpas_stream_mgr_read(domain % streamManager, whence=MPAS_STREAM_NEAREST, ierr=ierr)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', state)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'model', model)
         call mpas_pool_get_dimension(mesh, 'nCells', nCells)
         call mpas_pool_get_array(state, 'topo', topo)
         call mpas_pool_get_array(mesh, 'latCell', latCell)
         call mpas_pool_get_array(mesh, 'lonCell', lonCell)
         call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
         call mpas_pool_get_array(mesh, 'var2d', var2d)
         call mpas_pool_get_array(mesh, 'con', con)
         call mpas_pool_get_array(mesh, 'ol1', ol1)
         call mpas_pool_get_array(mesh, 'ol2', ol2)
         call mpas_pool_get_array(mesh, 'ol3', ol3)
         call mpas_pool_get_array(mesh, 'ol4', ol4)
         call mpas_pool_get_array(mesh, 'oa1', oa1)
         call mpas_pool_get_array(mesh, 'oa2', oa2)
         call mpas_pool_get_array(mesh, 'oa3', oa3)
         call mpas_pool_get_array(mesh, 'oa4', oa4)
         call mpas_pool_get_array(model, 'xtime', xtime)

         write(xtime,'(a)') trim(timeStamp)

         write(0,*) 'nCells = ', nCells
         write(0,*) 'min(topo)=', minval(topo)
         write(0,*) 'max(topo)=', maxval(topo)
         write(0,*) 'pointval=', topo(5,5)

         do iCell=1,nCells

            dc = 0.0
            do i=1,nEdgesOnCell(iCell)
               dc = dc + dcEdge(edgesOnCell(i,iCell))
            end do
            dc = dc / real(nEdgesOnCell(iCell),RKIND)
dc = 2.0 * dc

            call get_box(latCell(iCell)*rad2deg, lonCell(iCell)*rad2deg, dc)

            var2d(iCell) = get_var()
            con(iCell) = get_con()
            oa1(iCell) = get_oa1()
            oa2(iCell) = get_oa2()
            oa3(iCell) = get_oa3()
            oa4(iCell) = get_oa4()
            ol1(iCell) = get_ol1()
            ol2(iCell) = get_ol2()
            ol3(iCell) = get_ol3()
            ol4(iCell) = get_ol4()
         end do

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


   subroutine get_box(lat, lon, dx)

      implicit none

      real (kind=RKIND), intent(in) :: lat, lon, dx

      integer :: i, j, ii, jj, ic, jc

      if (cos(lat/rad2deg) > (2.0 * pts_per_degree * dx * 180.0) / (real(topo_x,RKIND) * Pi * Re)) then
         nx = ceiling((180.0 * dx * pts_per_degree) / (Pi * Re * cos(lat/rad2deg)))
      else
         nx = topo_x / 2
      end if
      ny = ceiling((180.0 * dx * pts_per_degree) / (Pi * Re))

      ic = nint((lon - start_lon) * pts_per_degree) + 1
      jc = nint((lat - start_lat) * pts_per_degree) + 1

      if (ic <= 0) ic = ic + topo_x
      if (ic > topo_x) ic = ic - topo_x

      if (associated(box)) deallocate(box)
      allocate(box(nx,ny))

      box_mean = 0.0
      do j=1,ny
      do i=1,nx

         ii = i - nx/2 + ic
         jj = j - ny/2 + jc

         if (jj <= 0) then
            jj = -jj + 1
            ii = ii + topo_y
         end if
         if (jj > topo_y) then
            jj = topo_y - (jj - topo_y - 1)
            ii = ii + topo_y
         end if
         do while (ii <= 0)
            ii = ii + topo_x
         end do
         do while (ii > topo_x)
            ii = ii - topo_x
         end do
 
         box(i,j) = topo(ii,jj)
         box_mean = box_mean + box(i,j)

      end do
      end do
      box_mean = box_mean / real(nx*ny, RKIND)

   end subroutine get_box


   real (kind=RKIND) function get_var()

      implicit none

      integer :: i, j
      real (kind=RKIND) :: s2

      s2 = 0.0

      do j=1,ny
         do i=1,nx
            s2 = s2 + (box(i,j) - box_mean)**2
         end do
      end do

      get_var = sqrt(s2 / real(nx*ny,RKIND))

   end function get_var


   real (kind=RKIND) function get_con()

      implicit none

      integer :: i, j
      real (kind=RKIND) :: s2, s4, var

      s2 = 0.0
      s4 = 0.0

      do j=1,ny
         do i=1,nx
            s2 = s2 + (box(i,j) - box_mean)**2
            s4 = s4 + (box(i,j) - box_mean)**4
         end do
      end do

      var = s2 / real(nx*ny,RKIND)
 
      if (abs(var) < 1.0e-5) then
         get_con = 0.0
      else
         get_con = s4 / (var**2 * real(nx*ny,RKIND))
      end if

   end function get_con


   real (kind=RKIND) function get_oa1()

      implicit none

      integer :: i, j
      integer :: nu, nd

      nu = 0
      nd = 0
      do j=1,ny
         do i=1,nx/2
            if (box(i,j) > box_mean) nu = nu + 1
         end do
         do i=nx/2+1,nx
            if (box(i,j) > box_mean) nd = nd + 1
         end do
      end do

      if (nu + nd > 0) then
         get_oa1 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa1 = 0.0
      end if

   end function get_oa1


   real (kind=RKIND) function get_oa2()

      implicit none

      integer :: i, j
      integer :: nu, nd

      nu = 0
      nd = 0
      do j=1,ny/2
         do i=1,nx
            if (box(i,j) > box_mean) nu = nu + 1
         end do
      end do
      do j=ny/2+1,ny
         do i=1,nx
            if (box(i,j) > box_mean) nd = nd + 1
         end do
      end do

      if (nu + nd > 0) then
         get_oa2 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa2 = 0.0
      end if

   end function get_oa2


   real (kind=RKIND) function get_oa3()

      implicit none

      integer :: i, j
      integer :: nu, nd
      real (kind=RKIND) :: ratio

      nu = 0
      nd = 0
      ratio = real(ny,RKIND)/real(nx,RKIND)
      do j=1,ny
         do i=1,nx
            if (nint(real(i,RKIND) * ratio) < (ny - j)) then
               if (box(i,j) > box_mean) nu = nu + 1
            else
               if (box(i,j) > box_mean) nd = nd + 1
            end if
         end do
      end do

      if (nu + nd > 0) then
         get_oa3 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa3 = 0.0
      end if

   end function get_oa3


   real (kind=RKIND) function get_oa4()

      implicit none

      integer :: i, j
      integer :: nu, nd
      real (kind=RKIND) :: ratio

      nu = 0
      nd = 0
      ratio = real(ny,RKIND)/real(nx,RKIND)
      do j=1,ny
         do i=1,nx
            if (nint(real(i,RKIND) * ratio) < j) then
               if (box(i,j) > box_mean) nu = nu + 1
            else
               if (box(i,j) > box_mean) nd = nd + 1
            end if
         end do
      end do

      if (nu + nd > 0) then
         get_oa4 = real((nu - nd),RKIND) / real((nu + nd),RKIND)
      else
         get_oa4 = 0.0
      end if

   end function get_oa4


   real (kind=RKIND) function get_ol1()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      do j=ny/4,3*ny/4
         do i=1,nx
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol1 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol1


   real (kind=RKIND) function get_ol2()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      do j=1,ny
         do i=nx/4,3*nx/4
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol2 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol2


   real (kind=RKIND) function get_ol3()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      do j=1,ny/2
         do i=1,nx/2
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do
      do j=ny/2+1,ny
         do i=nx/2+1,nx
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol3 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol3


   real (kind=RKIND) function get_ol4()

      implicit none

      integer :: i, j
      integer :: nw
      integer :: nt

      nw = 0
      nt = 0

      get_ol4 = real(nw,RKIND) / real(nx*ny,RKIND)
      do j=ny/2+1,ny
         do i=1,nx/2
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do
      do j=1,ny/2
         do i=nx/2+1,nx
            if (box(i,j) > box_mean) nw = nw + 1
            nt = nt + 1
         end do
      end do

      get_ol4 = real(nw,RKIND) / real(nt,RKIND)

   end function get_ol4

end module postproc_core
