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
