#if defined(mpas)

! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module machine

 implicit none
 save

 integer, parameter :: kind_io4  = selected_real_kind(6),      &
                       kind_io8  = selected_real_kind(12),     &
                       kind_ior  = selected_real_kind(12),     &
                       kind_evod = selected_real_kind(12),     &
                       kind_dbl_prec = selected_real_kind(12), &
                       kind_rad  = selected_real_kind(12),     & ! maps to 64-bit real
                       kind_phys = selected_real_kind(12),     & ! maps to 64-bit real
                       kind_REAL = selected_real_kind(12),     & ! used in cmp_comm
                       kind_INTEGER = selected_int_kind(9)

!==================================================================================================
 end module machine
!==================================================================================================

#else

      MODULE MACHINE

      IMPLICIT NONE

#ifndef SINGLE_PREC
      integer, parameter :: kind_io4  = 4, kind_io8  = 8 , kind_ior = 8 &
     &,                     kind_evod = 8, kind_dbl_prec = 8            &
     &,                     kind_qdt_prec = 16                          &
     &,                     kind_rad  = 8                               &
     &,                     kind_phys = 8     ,kind_taum=8              &
     &,                     kind_grid = 8                               &
     &,                     kind_REAL = 8                               &! used in cmp_comm
     &,                     kind_INTEGER = 4                             ! -,,-

#else
      integer, parameter :: kind_io4  = 4, kind_io8  = 8 , kind_ior = 8 &
     &,                     kind_evod = 4, kind_dbl_prec = 8            &
     &,                     kind_qdt_prec = 16                          &
     &,                     kind_rad  = 4                               &
     &,                     kind_phys = 4     ,kind_taum=4              &
     &,                     kind_grid = 4                               &
     &,                     kind_REAL = 4                               &! used in cmp_comm
     &,                     kind_INTEGER = 4                             ! -,,-

#endif

!
      real(kind=kind_evod), parameter :: mprec = 1.e-12           ! machine precision to restrict dep
      real(kind=kind_evod), parameter :: grib_undef = 9.99e20     ! grib undefine value
!
      END MODULE MACHINE

#endif

