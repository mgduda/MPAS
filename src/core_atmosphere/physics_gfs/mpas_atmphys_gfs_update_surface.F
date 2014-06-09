! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
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
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine physics_gfs_update_surface(current_date,mesh,sfc_input,diag_physics)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 character(len=*),intent(in):: current_date

!inout variables:
 type(mpas_pool_type),intent(inout):: diag_physics
 
!local pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer:: xice
 real(kind=RKIND),dimension(:),pointer:: alnsf,alnwf,alvsf,alvwf,vegfra,zorl
 real(kind=RKIND),dimension(:,:),pointer:: alnsf4s,alnwf4s,alvsf4s,alvwf4s,greenfrac,zorl12m

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter physics_gfs_update_surface:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(sfc_input,'landmask' ,landmask )
 call mpas_pool_get_array(sfc_input,'alnsf4s'  ,alnsf4s  )
 call mpas_pool_get_array(sfc_input,'alnwf4s'  ,alnwf4s  )
 call mpas_pool_get_array(sfc_input,'alvsf4s'  ,alvsf4s  )
 call mpas_pool_get_array(sfc_input,'alvwf4s'  ,alvwf4s  )
 call mpas_pool_get_array(sfc_input,'greenfrac',greenfrac)
 call mpas_pool_get_array(sfc_input,'xice'     ,xice     )
 call mpas_pool_get_array(sfc_input,'zorl12m'  ,zorl12m  )

 call mpas_pool_get_array(diag_physics,'alnsf' ,alnsf )
 call mpas_pool_get_array(diag_physics,'alnwf' ,alnwf )
 call mpas_pool_get_array(diag_physics,'alvsf' ,alvsf )
 call mpas_pool_get_array(diag_physics,'alvwf' ,alvwf )
 call mpas_pool_get_array(diag_physics,'vegfra',vegfra)
 call mpas_pool_get_array(diag_physics,'zorl'  ,zorl  )

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

!=================================================================================================================
 end module mpas_atmphys_gfs_update_surface
!=================================================================================================================


