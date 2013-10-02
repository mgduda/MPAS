      module tracer_const
      use machine , only : kind_phys
      implicit none
      SAVE

      real(kind=kind_phys) ri(0:20),cpi(0:20)
      integer, parameter :: num_tracer=3

      contains
! -------------------------------------------------------------------   
#if defined(mpas)
      subroutine set_tracer_const(me,ntrac)
      use machine , only : kind_phys
      use physcons , only : rd => con_rd , cpd => con_cp
      implicit none
      integer,intent(in):: ntrac,me
#else
      subroutine set_tracer_const (ntrac,me,nlunit)
      use machine , only : kind_phys
      use physcons , only : rd => con_rd , cpd => con_cp
      implicit none
      integer ntrac,me,nlunit
#endif
      namelist /tracer_constant/ ri,cpi

!
      if( ntrac.ne.num_tracer ) then
        if( me.eq.0 ) then
          write(*,*) ' Error ; inconsistent number of tracer '
          write(*,*) ' ntrac=',ntrac,' num_tracer=',num_tracer
        endif
        call abort
      endif

      ri=0.0
      cpi=0.0
      ri(0)=rd
      cpi(0)=cpd

#if defined(mpas)
      return
      end subroutine set_tracer_const
#else
      rewind(nlunit)
      read(nlunit, tracer_constant)

      return
      end subroutine set_tracer_const
#endif

      end module tracer_const
