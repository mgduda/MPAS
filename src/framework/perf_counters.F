module perf_counters

   use iso_c_binding

   integer :: first_call

!#include "fpapi.h"

   contains


   subroutine perf_counters_init()

      implicit none

      integer(c_int) :: ierr

!      call PAPIF_library_init(ierr)
!write(0,*) 'PERF: PAPIF_library_init ', ierr

      first_call = 1

   end subroutine perf_counters_init


   subroutine perf_counters_finalize()

      implicit none

!      call PAPIF_shutdown()
!write(0,*) 'PERF: PAPIF_shutdown '

   end subroutine perf_counters_finalize


   subroutine start_counters()

      implicit none

      integer(c_int) :: ierr
      integer(c_long_long) :: flpops
      real(c_float) :: real_time, proc_time, mflops
      integer(c_long_long), dimension(1) :: values
      integer(c_int) :: len

      len = 1
!      call PAPIF_stop_counters(values, len, ierr)
!write(0,*) 'PERF: ierr=', ierr
!      call PAPIF_flops(real_time, proc_time, flpops, mflops, ierr)
!write(0,*) 'PERF: ierr=', ierr

   end subroutine start_counters


   subroutine stop_counters(arg)
  
      implicit none

      character (len=*), intent(in) :: arg

      integer(c_int) :: ierr
      integer(c_long_long) :: flpops
      real(c_float) :: real_time, proc_time, mflops

!      call PAPIF_flops(real_time, proc_time, flpops, mflops, ierr)
!write(0,'(a,f9.4,a,i10,a,f9.4,a,i3)') 'PERF '//trim(arg)//' : rt/flpops/mflops=', real_time, ' ', flpops, ' ', mflops, ' ', ierr

   end subroutine stop_counters

end module perf_counters
