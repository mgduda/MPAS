! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_state_fields_out
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_vars
 use nuopc_physics

 implicit none
 private

 public:: mpas_state_fields_out

 type,extends(state_fields_out):: mpas_state_fields_out

    contains
       procedure:: mpas_setphys   => mpas_state_fld_out_setphys
       procedure:: allocate_mem   => mpas_state_fld_out_allocate
       procedure:: deallocate_mem => mpas_state_fld_out_deallocate

 end type

!modifications to sourcecode:
!----------------------------
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.

 contains


!=================================================================================================================
 subroutine mpas_state_fld_out_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_state_fields_out),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve
 integer,pointer:: nVertLevels
 integer:: ntrac

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_state_fld_out_allocate:'

 call mpas_pool_get_dimension(mesh,"nCellsSolve",nCellsSolve)
 call mpas_pool_get_dimension(mesh,"nVertLevels",nVertLevels)

 ntrac = mdl_parameters%ntrac

 if(.not.allocated(gt0)) allocate(gt0(1:nCellsSolve,1:nVertLevels))
 if(.not.allocated(gu0)) allocate(gu0(1:nCellsSolve,1:nVertLevels))
 if(.not.allocated(gv0)) allocate(gv0(1:nCellsSolve,1:nVertLevels))

 if(.not.allocated(gq0)) allocate(gq0(1:nCellsSolve,1:nVertLevels,ntrac))

!write(0,*) '--- end subroutine mpas_state_fld_out_allocate'

 end subroutine mpas_state_fld_out_allocate

!=================================================================================================================
 subroutine mpas_state_fld_out_deallocate(this)
!=================================================================================================================

 class(mpas_state_fields_out),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_state_fld_out_deallocate:'

 if(allocated(gt0)) deallocate(gt0)
 if(allocated(gu0)) deallocate(gu0)
 if(allocated(gv0)) deallocate(gv0)
 if(allocated(gq0)) deallocate(gq0)

!write(0,*) '--- end subroutine mpas_state_fld_out_allocate'

 end subroutine mpas_state_fld_out_deallocate

!=================================================================================================================
 subroutine mpas_state_fld_out_setphys(this,mesh)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh

!inout variables:
 class(mpas_state_fields_out),intent(inout):: this

!local variables and pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer:: iCell,k,n,ntrac

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_state_fld_setphys_out:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 ntrac = mdl_parameters%ntrac

 do k = 1, nVertLevels
    do iCell = 1, nCellsSolve
       gt0(iCell,k) = tgrs(iCell,k)
       gu0(iCell,k) = ugrs(iCell,k)
       gv0(iCell,k) = vgrs(iCell,k)
       do n = 1, ntrac
          gq0(iCell,k,n) = qgrs(iCell,k,n)
       enddo
    enddo
 enddo

 call this%setphys(gt0,gq0,gu0,gv0)

!write(0,*) '--- end subroutine mpas_state_fld_setphys_out'

 end subroutine mpas_state_fld_out_setphys

!=================================================================================================================
 end module mpas_atmphys_nggps_state_fields_out
!=================================================================================================================
