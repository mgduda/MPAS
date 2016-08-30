! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_io

#ifndef NO_PIO
#include "mpas_io_pio.inc"
#else
#include "mpas_io_nc4.inc"
#endif
 
end module mpas_io
