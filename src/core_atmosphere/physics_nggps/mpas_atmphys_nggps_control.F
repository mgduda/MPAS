! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_control
 use mpas_derived_types
 use mpas_pool_routines

 implicit none
 private
 public:: physics_namelist_check, &
          physics_registry_init

 logical,public:: moist_physics


!>\brief MPAS control and initialization routines.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2015-08-01.
!>


 contains

!==================================================================================================
 subroutine physics_namelist_check(mesh,configs)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine physics_namelist_check NGGPS:'
 moist_physics = .true.
 write(0,*) '--- end subroutine physics_namelist_check NGGPS.'

 end subroutine physics_namelist_check

!==================================================================================================
 subroutine physics_registry_init(mesh,configs,sfc_input)
!==================================================================================================

!input and inout arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(inout):: sfc_input
 
!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine physics_registry_init NGGPS:'
 write(0,*) '--- end subroutine physics_registry_init NGGPS.'

 end subroutine physics_registry_init

!==================================================================================================
 end module mpas_atmphys_control
!==================================================================================================


