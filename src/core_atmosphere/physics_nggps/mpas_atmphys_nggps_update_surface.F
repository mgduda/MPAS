! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_update_surface
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_date_time
 use mpas_atmphys_nggps_date_time

 implicit none
 private
 public:: physics_nggps_update_surface


 contains


!=================================================================================================================
 subroutine physics_nggps_update_surface(current_date,mesh,nggps_sfc_properties)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh
 character(len=*),intent(in):: current_date

!inout variables:
 type(mpas_pool_type),intent(inout):: nggps_sfc_properties
 
!local variables:
 integer:: iCell
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=R8KIND),dimension(:),pointer:: alnsf,alnwf,alvsf,alvwf,vfrac,zorl
 real(kind=R8KIND),dimension(:,:),pointer:: alnsf4s,alnwf4s,alvsf4s,alvwf4s,greenfrac,zorl12m

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter physics_nggps_update_surface:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(nggps_sfc_properties,'landmask' ,landmask )
 call mpas_pool_get_array(nggps_sfc_properties,'alnsf4s'  ,alnsf4s  )
 call mpas_pool_get_array(nggps_sfc_properties,'alnwf4s'  ,alnwf4s  )
 call mpas_pool_get_array(nggps_sfc_properties,'alvsf4s'  ,alvsf4s  )
 call mpas_pool_get_array(nggps_sfc_properties,'alvwf4s'  ,alvwf4s  )
 call mpas_pool_get_array(nggps_sfc_properties,'greenfrac',greenfrac)
 call mpas_pool_get_array(nggps_sfc_properties,'zorl12m'  ,zorl12m  )

 call mpas_pool_get_array(nggps_sfc_properties,'alnsf',alnsf)
 call mpas_pool_get_array(nggps_sfc_properties,'alnwf',alnwf)
 call mpas_pool_get_array(nggps_sfc_properties,'alvsf',alvsf)
 call mpas_pool_get_array(nggps_sfc_properties,'alvwf',alvwf)
 call mpas_pool_get_array(nggps_sfc_properties,'vfrac',vfrac)
 call mpas_pool_get_array(nggps_sfc_properties,'zorl' ,zorl )

!initialization of the green-ness (vegetation) fraction: interpolation of the monthly mean values
!to the initial date:
 call monthly_interp_to_date(nCellsSolve,current_date,greenfrac,vfrac)
 write(0,*) '--- end monthly_interp_to_date'
 write(0,*) '--- max vfrac = ',maxval(vfrac(1:nCellsSolve))
 write(0,*) '--- min vfrac = ',minval(vfrac(1:nCellsSolve))

!initialization of the surface roughness length (here zorl is in cm) interpolation of the monthly
!mean values to the current date. set the minimum values for zorl over the oceans.
!call monthly_interp_to_date(nCellsSolve,current_date,zorl12m,zorl)
!do iCell = 1, nCellsSolve
!   if(landmask(iCell) == 1 .and. zorl(iCell) == 0._R8KIND) then 
!      zorl(iCell) = 2.0_R8KIND
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
 write(0,*) '--- max alnwf = ',maxval(alnwf(1:nCellsSolve))
 write(0,*) '--- max alvsf = ',maxval(alvsf(1:nCellsSolve))
 write(0,*) '--- max alvwf = ',maxval(alvwf(1:nCellsSolve))

 write(0,*) '--- end physics_nggps_update_surface.'
 write(0,*)

 end subroutine physics_nggps_update_surface

!=================================================================================================================
 end module mpas_atmphys_nggps_update_surface
!=================================================================================================================


