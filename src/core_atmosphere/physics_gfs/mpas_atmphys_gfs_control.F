! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_control
 use mpas_grid_types

 implicit none
 private
 public:: physics_namelist_check, &
          physics_registry_init

 logical,public:: moist_physics


!>\brief MPAS control and initialization routines. Note that the length of the module is strongly reduced relative
!> to that in the physics directory. It is only used to follow the implementation of the WRF physics, and use the
!> same set of calls in the module mpas_atm_mpas_core.F.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-01-27.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine physics_namelist_check(mesh,configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine physics_namelist_check:'
 moist_physics = .true.
 write(0,*) '--- end subroutine physics_namelist_check.'

 end subroutine physics_namelist_check

!=================================================================================================================
 subroutine physics_registry_init(mesh,configs,sfc_input)
!=================================================================================================================

!input and inout arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(inout):: sfc_input
 
!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine physics_registry_init:'
 write(0,*) '--- end subroutine physics_registry_init.'

 end subroutine physics_registry_init

!=================================================================================================================
 end module mpas_atmphys_control
!=================================================================================================================


