! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_update
 use mpas_configure
 use mpas_grid_types

 implicit none
 private
 public:: physics_update, &
          physics_gfs_update_diagnostics


!>\brief update diagnostics.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-01-28.
!>
!>\details
!>
!> subroutines in mpas_atmphys_update:
!> -----------------------------------
!> physics_update                : not used.
!> physics_gfs_update_diagnostics: reset all accumulated diagnostics to zero.


 contains
 
!==================================================================================================
 subroutine physics_update(domain,dt)
!==================================================================================================

!input arguments:
 type(domain_type),intent(in):: domain
 real(kind=RKIND),intent(in) :: dt 

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- begin physics_update:'
 write(0,*) '--- end physics_update.'

 end subroutine physics_update

!==================================================================================================
 subroutine physics_gfs_update_diagnostics(mesh,diag_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout)::diag_physics

!local variables:
 integer:: iCell,k

!--------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine physics_gfs_update_diagnostics:'
!do iCell = 1,5
!   write(0,201) iCell,diag_physics%suntim%array(iCell),diag_physics%dlwsfc%array(iCell), &
!                diag_physics%dusfc%array(iCell)
!enddo

 do iCell = 1, mesh%nCellsSolve
    diag_physics % psmean  % array(iCell) = 0._RKIND
    diag_physics % suntim  % array(iCell) = 0._RKIND
    diag_physics % dlwsfc  % array(iCell) = 0._RKIND
    diag_physics % ulwsfc  % array(iCell) = 0._RKIND
    diag_physics % gflux   % array(iCell) = 0._RKIND
    diag_physics % evbsa   % array(iCell) = 0._RKIND
    diag_physics % evcwa   % array(iCell) = 0._RKIND
    diag_physics % transa  % array(iCell) = 0._RKIND
    diag_physics % sbsnoa  % array(iCell) = 0._RKIND
    diag_physics % snowca  % array(iCell) = 0._RKIND
    diag_physics % snohfa  % array(iCell) = 0._RKIND
    diag_physics % ep      % array(iCell) = 0._RKIND
    diag_physics % gtmp2m  % array(iCell) = 0._RKIND
    diag_physics % gu10m   % array(iCell) = 0._RKIND
    diag_physics % gv10m   % array(iCell) = 0._RKIND
    diag_physics % gustar  % array(iCell) = 0._RKIND
    diag_physics % gzorl   % array(iCell) = 0._RKIND
    diag_physics % goro    % array(iCell) = 0._RKIND
    diag_physics % dusfc   % array(iCell) = 0._RKIND
    diag_physics % dvsfc   % array(iCell) = 0._RKIND
    diag_physics % dtsfc   % array(iCell) = 0._RKIND
    diag_physics % dqsfc   % array(iCell) = 0._RKIND
    diag_physics % gpblh   % array(iCell) = 0._RKIND
    diag_physics % dugwd   % array(iCell) = 0._RKIND
    diag_physics % dvgwd   % array(iCell) = 0._RKIND
    diag_physics % cldwrk  % array(iCell) = 0._RKIND
    diag_physics % cnvprcp % array(iCell) = 0._RKIND
    diag_physics % totprcp % array(iCell) = 0._RKIND
    diag_physics % runoff  % array(iCell) = 0._RKIND
    diag_physics % srunoff % array(iCell) = 0._RKIND
    diag_physics % gsoil   % array(iCell) = 0._RKIND

    diag_physics % tmpmax  % array(iCell) = 0._RKIND
    diag_physics % tmpmin  % array(iCell) = 0._RKIND
    diag_physics % spfhmax % array(iCell) = 0._RKIND
    diag_physics % spfhmin % array(iCell) = 0._RKIND

    do k = 1, mesh%nVertLevels
       diag_physics % dkh % array(k,iCell) = 0._RKIND
    enddo

 enddo

 write(0,*) '--- end subroutine physics_gfs_update_diagnostics:'
!do iCell = 1,5
!   write(0,201) iCell,diag_physics%suntim%array(iCell),diag_physics%dlwsfc%array(iCell), &
!                diag_physics%dusfc%array(iCell)
!enddo
!write(0,*)
 201 format(i9,10(1x,e15.8))

 end subroutine physics_gfs_update_diagnostics

!==================================================================================================
 end module mpas_atmphys_update
!==================================================================================================
