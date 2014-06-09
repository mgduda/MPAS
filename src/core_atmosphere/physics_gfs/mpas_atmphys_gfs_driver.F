! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver
 use mpas_dmpar
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
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine physics_driver(domain,itimestep,xtime_s)
!=================================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(domain_type),intent(inout):: domain

!local pointers:
 type(mpas_pool_type),pointer::  mesh,         &
                                 state,        &
                                 diag,         &
                                 diag_physics, &
                                 tend_physics, &
                                 atm_input,    &
                                 sfc_input

!local variables:
 type(block_type),pointer:: block

 integer:: time_lev

!=================================================================================================================
 write(0,*)
 write(0,*) '--- enter subroutine physics_gfs_driver:'

 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'        ,mesh        )
    call mpas_pool_get_subpool(block%structs,'state'       ,state       )
    call mpas_pool_get_subpool(block%structs,'diag'        ,diag        )
    call mpas_pool_get_subpool(block%structs,'diag_physics',diag_physics)
    call mpas_pool_get_subpool(block%structs,'atm_input'   ,atm_input   )
    call mpas_pool_get_subpool(block%structs,'sfc_input'   ,sfc_input   )
    call mpas_pool_get_subpool(block%structs,'tend_physics',tend_physics)

    !allocate arrays shared by all physics parameterizations:
    call allocate_forall_physics(mesh)

    !physics prep step:
    time_lev = 1
    call MPAS_to_physics(mesh,state,time_lev,diag)

    !call to longwave and shortwave radiation:
    if(l_radtlw .or. l_radtsw) then
       call allocate_radiation(mesh)
       call mpas_atmphys_grrad(itimestep,xtime_s,domain%dminfo,mesh,sfc_input, &
                               diag_physics,tend_physics)
       call deallocate_radiation
    endif

    !call to all physics, except radiation:
    call allocate_gbphys(mesh) 
    call mpas_atmphys_gbphys(itimestep,xtime_s,domain%dminfo,mesh,atm_input, &
                             sfc_input,diag_physics,tend_physics)
    call deallocate_gbphys

    !deallocate arrays shared by all physics parameterizations:
    call deallocate_forall_physics

    block => block % next

 end do

 write(0,*) '--- end subroutine physics_gfs_driver'
 write(0,*)

 end subroutine physics_driver

!=================================================================================================================
 end module mpas_atmphys_driver
!=================================================================================================================
