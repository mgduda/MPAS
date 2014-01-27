! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_gfs_update_surface
 use mpas_grid_types
 use mpas_kind_types

 use mpas_atmphys_date_time
 use mpas_atmphys_gfs_date_time

 implicit none
 private
 public:: physics_gfs_update_surface

!>\brief update surface boundary conditions.
!>\author Laura D. Fowler (send comments to laura@mmm.ucar.edu).
!>\date 2013-11-20.
!>
!>\details
!>
!> subroutines in mpas_atmphys_update_surface:
!> -------------------------------------------
!> physics_gfs_update_surface : update climatological spectral surface albedo,greeness fraction,
!>                              and mixing length.
!>


 contains

!==================================================================================================
 subroutine physics_gfs_update_surface(current_date,mesh,sfc_input,diag_physics)
!==================================================================================================

!input variables:
 type(mesh_type),intent(in):: mesh
 type(sfc_input_type),intent(in):: sfc_input
 character(len=*),intent(in):: current_date

!inout variables:
 type(diag_physics_type),intent(inout):: diag_physics
 
!local variables:
 integer:: iCell,nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer:: xice
 real(kind=RKIND),dimension(:),pointer:: alnsf,alnwf,alvsf,alvwf,vegfra,zorl
 real(kind=RKIND),dimension(:,:),pointer:: alnsf4s,alnwf4s,alvsf4s,alvwf4s,greenfrac,zorl12m

!--------------------------------------------------------------------------------------------------
 write(0,*) '--- enter physics_gfs_update_surface:'

 nCellsSolve = mesh % nCellsSolve

 landmask  => sfc_input % landmask  % array
 alnsf4s   => sfc_input % alnsf4s   % array
 alnwf4s   => sfc_input % alnwf4s   % array
 alvsf4s   => sfc_input % alvsf4s   % array
 alvwf4s   => sfc_input % alvwf4s   % array
 greenfrac => sfc_input % greenfrac % array
 xice      => sfc_input % xice      % array
 zorl12m   => sfc_input % zorl12m   % array

 alnsf  => diag_physics % alnsf  % array
 alnwf  => diag_physics % alnwf  % array
 alvsf  => diag_physics % alvsf  % array
 alvwf  => diag_physics % alvwf  % array
 vegfra => diag_physics % vegfra % array
 zorl   => diag_physics % zorl   % array

!initialization of the green-ness (vegetation) fraction: interpolation of the monthly mean values
!to the initial date:
 call monthly_interp_to_date(nCellsSolve,current_date,greenfrac,vegfra)
 write(0,*) '--- end monthly_interp_to_date'
 write(0,*) '--- max vegfra = ',maxval(vegfra(1:nCellsSolve))
 write(0,*) '--- min vegfra = ',minval(vegfra(1:nCellsSolve))

!initialization of the surface roughness length (here zorl is in cm) interpolation of the monthly
!mean values to the current date. set the minimum values for zorl over the oceans.
!call monthly_interp_to_date(nCellsSolve,current_date,zorl12m,zorl)
!do iCell = 1, nCellsSolve
!   if(landmask(iCell) == 1 .and. zorl(iCell) == 0._RKIND) then 
!      zorl(iCell) = 2.0_RKIND
!   endif    
!enddo

!initialization of the seasonal mean near-infrared and visible spectral surface albedos with strong
!and weak dependence on the solar zenith angle
 call season_interp_to_date(nCellsSolve,current_date,alnsf4s,alnsf)
 call season_interp_to_date(nCellsSolve,current_date,alnwf4s,alnwf)
 call season_interp_to_date(nCellsSolve,current_date,alvsf4s,alvsf)
 call season_interp_to_date(nCellsSolve,current_date,alvwf4s,alvwf)
 write(0,*) '--- end season_interp_to_date'
 write(0,*) '--- max alnsf = ',maxval(alnsf(1:nCellsSolve))
 write(0,*) '--- min alnwf = ',minval(alnwf(1:nCellsSolve))
 write(0,*) '--- max alvsf = ',maxval(alvsf(1:nCellsSolve))
 write(0,*) '--- min alvwf = ',minval(alvwf(1:nCellsSolve))

 write(0,*) '--- end physics_gfs_update_surface:'
 write(0,*)

 end subroutine physics_gfs_update_surface

!==================================================================================================
 end module mpas_atmphys_gfs_update_surface
!==================================================================================================


