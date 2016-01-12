! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_cloud_properties
 use mpas_atmphys_nggps_diagnostics
 use mpas_atmphys_nggps_driver_out
 use mpas_atmphys_nggps_dynamic_parameters
 use mpas_atmphys_nggps_interface_fields
 use mpas_atmphys_nggps_radiation_tendencies
 use mpas_atmphys_nggps_sfc_properties
 use mpas_atmphys_nggps_state_fields_in
 use mpas_atmphys_nggps_state_fields_out
 use mpas_atmphys_nggps_tbd_ddt
 use mpas_atmphys_nggps_vars
 use nuopc_physics

 implicit none
 private
 public:: physics_driver

!modifications to sourcecode:
!----------------------------
! * added the mpas_pool_type nggps_dynamic_parameters to the argument list of subroutines mpas_dyn_param_setrad
!   and mpas_dyn_param_setphys.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-28.


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
 type(mpas_pool_type),pointer:: configs,                    &
                                mesh,                       &
                                state,                      &
                                diag,                       &
                                nggps_atm_inputs,           & 
                                nggps_cloud_properties,     &
                                nggps_diagnostics,          &
                                nggps_dynamic_parameters,   &
                                nggps_interface_fields,     &
                                nggps_physics_tendencies,   &
                                nggps_radiation_tendencies, &
                                nggps_sfc_properties,       &
                                nggps_tbd_ddt
 type(block_type),pointer:: block

!local variables:
 integer:: time_lev,iCell,k

 type(mpas_state_fields_in)     :: mpas_statein
 type(mpas_state_fields_out)    :: mpas_stateout
 type(mpas_cloud_properties)    :: mpas_cld_props
 type(mpas_diagnostics)         :: mpas_diags
 type(mpas_dynamic_parameters)  :: mpas_dyn_parameters 
 type(mpas_interface_fields)    :: mpas_intrfc_flds
 type(mpas_radiation_tendencies):: mpas_rad_tendencies
 type(mpas_sfc_properties)      :: mpas_sfc_prop
 type(mpas_tbd_ddt)             :: mpas_tbds

!=================================================================================================================
 write(0,*)
 write(0,*) '--- enter subroutine physics_ngpps_driver:'

 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh' ,mesh   )
    call mpas_pool_get_subpool(block%structs,'state',state  )
    call mpas_pool_get_subpool(block%structs,'diag' ,diag   )
    call mpas_pool_get_subpool(block%structs,'nggps_atm_inputs'          ,nggps_atm_inputs          )
    call mpas_pool_get_subpool(block%structs,'nggps_cloud_properties'    ,nggps_cloud_properties    )
    call mpas_pool_get_subpool(block%structs,'nggps_diagnostics'         ,nggps_diagnostics         )
    call mpas_pool_get_subpool(block%structs,'nggps_dynamic_parameters'  ,nggps_dynamic_parameters  )
    call mpas_pool_get_subpool(block%structs,'nggps_interface_fields'    ,nggps_interface_fields    )
    call mpas_pool_get_subpool(block%structs,'nggps_physics_tendencies'  ,nggps_physics_tendencies  )
    call mpas_pool_get_subpool(block%structs,'nggps_radiation_tendencies',nggps_radiation_tendencies)
    call mpas_pool_get_subpool(block%structs,'nggps_sfc_properties'      ,nggps_sfc_properties      )
    call mpas_pool_get_subpool(block%structs,'nggps_tbd_ddt'             ,nggps_tbd_ddt             )

    !allocate arrays shared by the nggps radiation and physics parameterizations:
    time_lev = 1
    call mpas_statein%allocate_mem(mesh)
    call mpas_statein%mpas_set(mesh,state,time_lev,diag)

    call mpas_cld_props%allocate_mem(mesh)
    call mpas_diags%allocate_mem(mesh)
    call mpas_intrfc_flds%allocate_mem(mesh)
    call mpas_rad_tendencies%allocate_mem(mesh)
    call mpas_sfc_prop%allocate_mem(mesh)

    !... call to the NGGPS interface between MPAS and the GFS radiation:
    call mpas_cld_props%mpas_setrad(mesh,nggps_cloud_properties)
    call mpas_diags%mpas_setrad(mesh,nggps_diagnostics)
    call mpas_intrfc_flds%mpas_setrad(mesh,nggps_interface_fields)
    call mpas_rad_tendencies%mpas_setrad(mesh,nggps_radiation_tendencies)
    call mpas_sfc_prop%mpas_setrad(mesh,nggps_sfc_properties)

    if(l_radtlw .or. l_radtsw) then

!      call mpas_cld_props%mpas_setrad(mesh,nggps_cloud_properties)
!      call mpas_diags%mpas_setrad(mesh,nggps_diagnostics)
       call mpas_dyn_parameters%setrad_allocate_mem(mesh)
       call mpas_dyn_parameters%mpas_setrad(block%configs,mesh,nggps_dynamic_parameters,itimestep,xtime_s)
!      call mpas_intrfc_flds%mpas_setrad(mesh,nggps_interface_fields)
!      call mpas_rad_tendencies%mpas_setrad(mesh,nggps_radiation_tendencies)
!      call mpas_sfc_prop%mpas_setrad(mesh,nggps_sfc_properties)

       call nuopc_rad_run( &
                      mpas_statein%state_fields_in,             &
                      mpas_sfc_prop%sfc_properties,             &
                      mpas_diags%diagnostics,                   &
                      mpas_intrfc_flds%interface_fields,        &
                      mpas_cld_props%cloud_properties,          &
                      mpas_rad_tendencies%radiation_tendencies, &
                      mdl_parameters,                           &
                      mpas_dyn_parameters%dynamic_parameters    &
                         )

       call nuopc_rad_write( &
                      mpas_cld_props%cloud_properties,          &
                      mpas_diags%diagnostics,                   &
                      mpas_intrfc_flds%interface_fields,        &
                      mpas_rad_tendencies%radiation_tendencies, &
                      mesh,                                     &
                      nggps_cloud_properties,                   &
                      nggps_diagnostics,                        &
                      nggps_interface_fields,                   &
                      nggps_radiation_tendencies                &
                           )

       call mpas_dyn_parameters%setrad_deallocate_mem

    endif

    call mpas_stateout%allocate_mem(mesh)
    call mpas_stateout%mpas_setphys(mesh)
    call mpas_cld_props%mpas_setphys(mesh,nggps_cloud_properties)
    call mpas_diags%mpas_setphys(mesh,nggps_diagnostics)
    call mpas_dyn_parameters%mpas_setphys(block%configs,mesh,nggps_dynamic_parameters,itimestep,xtime_s)
    call mpas_intrfc_flds%mpas_setphys(mesh,nggps_interface_fields)
    call mpas_sfc_prop%mpas_setphys(mesh,nggps_sfc_properties)
    call mpas_tbds%allocate_mem(mesh)
    call mpas_tbds%mpas_set(mesh,nggps_atm_inputs,nggps_tbd_ddt,mpas_statein%state_fields_in)

    !... call to the NGGPS interface between MPAS and the GFS physics:
    call nuopc_phys_run( &
                    mpas_statein%state_fields_in,             &
                    mpas_stateout%state_fields_out,           & 
                    mpas_sfc_prop%sfc_properties,             &
                    mpas_diags%diagnostics,                   &
                    mpas_intrfc_flds%interface_fields,        & 
                    mpas_cld_props%cloud_properties,          & 
                    mpas_rad_tendencies%radiation_tendencies, &
                    mdl_parameters,                           &
                    mpas_tbds%tbd_ddt,                        & 
                    mpas_dyn_parameters%dynamic_parameters)

    call nuopc_phys_write( &
                    mpas_statein%state_fields_in,   &
                    mpas_stateout%state_fields_out, &
                    mpas_diags%diagnostics,         &
                    mesh,                           &
                    nggps_physics_tendencies)

    !... back to Registry (to be moved below nuopc_phys_run:
    call mpas_cld_props%mpas_write(mesh,nggps_cloud_properties)
    call mpas_diags%mpas_write(mesh,nggps_diagnostics)
    call mpas_rad_tendencies%mpas_write(mesh,nggps_radiation_tendencies)
    call mpas_sfc_prop%mpas_write(nggps_sfc_properties)
    call mpas_tbds%mpas_write(mesh,nggps_tbd_ddt)

   !deallocate arrays shared by all physics parameterizations:
    call mpas_statein%deallocate_mem
    call mpas_stateout%deallocate_mem
    call mpas_cld_props%deallocate_mem
    call mpas_diags%deallocate_mem
    call mpas_intrfc_flds%deallocate_mem
    call mpas_rad_tendencies%deallocate_mem
    call mpas_sfc_prop%deallocate_mem
    call mpas_tbds%deallocate_mem

    block => block % next

 end do

 write(0,*)
 write(0,*) '--- end subroutine physics_nggps_driver'
 write(0,*)

 end subroutine physics_driver

!=================================================================================================================
 end module mpas_atmphys_driver
!=================================================================================================================
