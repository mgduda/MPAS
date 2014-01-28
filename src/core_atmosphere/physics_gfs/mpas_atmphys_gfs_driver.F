! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_driver
 use mpas_dmpar
 use mpas_configure
 use mpas_grid_types

 use mpas_atmphys_gfs_driver_gbphys
 use mpas_atmphys_gfs_driver_grrad
 use mpas_atmphys_gfs_interface
 use mpas_atmphys_gfs_vars, only: l_radtlw,l_radtsw

 implicit none
 private
 public:: physics_driver


!>\brief MPAS top GFS physics driver.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-01-27.
!>
!>\details
!> subroutine physics_driver is the top physics driver from which separate drivers for the GFS
!> long- and short-wave radiation codes (grrad) and all other physics parameterizations are 
!> are called.
!>
!> subroutines called in mpas_atmphys_driver:
!> ------------------------------------------
!> allocate_forall_physics     : allocate local arrays defining atmospheric soundings (pressure,..)
!> allocate_radiation          : allocate all local arrays used in mpas_atmphys_grrad.
!> allocate_gbphys             : allocate all local arrays used in mpas_atmphys_gbphys.
!>
!> deallocate_forall_physics   : deallocate local arrays defining atmospheric soundings.
!> deallocate_radiation        : deallocate all local arrays used in mpas_atmphys_grrad.
!> deallocate_gbphys           : deallocate all local arrays used in mpas_atmphys_gbphys.
!>
!> MPAS_to_physics             : interface between the dynamical core and GFS physics.
!> mpas_atmphys_grrad          : driver for GFS long-wave and short-wave radiation codes. 
!> mpas_atmphys_gbphys         : driver for all other GFS physics, except radiation.


 contains

!==================================================================================================
 subroutine physics_driver(domain,itimestep,xtime_s)
!==================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables:
 type(block_type),pointer:: block

!==================================================================================================
 write(0,*)
 write(0,*) '--- enter subroutine physics_gfs_driver:'

 block => domain % blocklist
 do while(associated(block))

    !allocate arrays shared by all physics parameterizations:
    call allocate_forall_physics(block%mesh)
    call MPAS_to_physics(block%mesh,block%state%time_levs(1)%state,block%diag)

    !call to longwave and shortwave radiation:
    if(l_radtlw .or. l_radtsw) then
       call allocate_radiation(block%mesh)
       call mpas_atmphys_grrad(itimestep,xtime_s,domain%dminfo,block%mesh,block%sfc_input, &
                               block%diag_physics,block%tend_physics)
       call deallocate_radiation
    endif

    !call to all physics, except radiation:
    call allocate_gbphys(block%mesh) 
    call mpas_atmphys_gbphys(itimestep,xtime_s,domain%dminfo,block%mesh,block%atm_input, &
                             block%sfc_input,block%diag_physics,block%tend_physics)
    call deallocate_gbphys

    !deallocate arrays shared by all physics parameterizations:
    call deallocate_forall_physics

    block => block % next

 end do

 write(0,*) '--- end subroutine physics_gfs_driver'
 write(0,*)

 end subroutine physics_driver

!==================================================================================================
 end module mpas_atmphys_driver
!==================================================================================================
