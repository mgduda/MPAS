! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_driver_out
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_vars
 use nuopc_physics

 use machine
 use module_radlw_parameters, only: sfcflw_type,topflw_type
 use module_radsw_parameters, only: sfcfsw_type,topfsw_type

 implicit none
 private
 public:: nuopc_rad_write, &
          nuopc_phys_write 

!modifications:
!--------------
! * in subroutine nuopc_phys_write, corrected the physics tendency for water vapor so that it is expressed in term
!   of a mixing ratio instead of in term of a specific humidity.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-19. 
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.
! * corrected the calculation of the tendencies for the ozone and cloud condensate "specific" mixing ratios to
!   mixing ratios.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-24.

 contains


!=================================================================================================================
 subroutine nuopc_rad_write(cld_properties,diags,intr_fields,rad_tendencies,mesh,nggps_cloud_properties, &
                            nggps_diagnostics,nggps_interface_fields,nggps_radiation_tendencies)
!=================================================================================================================

!input arguments:
 type(cloud_properties),intent(in)    :: cld_properties
 type(diagnostics),intent(in)         :: diags
 type(interface_fields),intent(in)    :: intr_fields
 type(radiation_tendencies),intent(in):: rad_tendencies
 
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_cloud_properties
 type(mpas_pool_type),intent(inout):: nggps_diagnostics
 type(mpas_pool_type),intent(inout):: nggps_interface_fields
 type(mpas_pool_type),intent(inout):: nggps_radiation_tendencies

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nFxr

 real(kind=PHYSKIND),dimension(:),pointer  :: sfcnirbmd,sfcnirdfd,sfcnirbmu,sfcnirdfu
 real(kind=PHYSKIND),dimension(:),pointer  :: sfcvisbmd,sfcvisdfd,sfcvisbmu,sfcvisdfu
 real(kind=PHYSKIND),dimension(:),pointer  :: coszdg,coszen,semis,sfalb,tsflw
 real(kind=PHYSKIND),dimension(:),pointer  :: sfclwup,sfclwupc,sfclwdn,sfclwdnc
 real(kind=PHYSKIND),dimension(:),pointer  :: sfcswup,sfcswupc,sfcswdn,sfcswdnc
 real(kind=PHYSKIND),dimension(:),pointer  :: toplwup,toplwupc,topswup,topswupc,topswdn
 real(kind=PHYSKIND),dimension(:,:),pointer:: cldcov,htrlw,htrsw
 real(kind=PHYSKIND),dimension(:,:),pointer:: fluxr,dswcmp,uswcmp

 real(kind=kind_phys),dimension(:,:),pointer:: fluxr_l
 real(kind=kind_phys),dimension(:,:),pointer:: cldcov_l,htrlw_l,htrsw_l
 real(kind=kind_phys),dimension(:,:),pointer:: dswcmp_l,uswcmp_l

 type(sfcflw_type),dimension(:),pointer:: sfcflw_l
 type(sfcfsw_type),dimension(:),pointer:: sfcfsw_l
 type(topflw_type),dimension(:),pointer:: topflw_l
 type(topfsw_type),dimension(:),pointer:: topfsw_l

!local variables:
 integer:: iCell,k

!diags%fluxr


!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine nuopc_rad_out:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nFxr'       ,nFxr       )

!--- cloud_properties:
 call mpas_pool_get_array(nggps_cloud_properties,'cldcov',cldcov)

!--- diagnostics:
 call mpas_pool_get_array(nggps_diagnostics,'toplwup' ,toplwup )
 call mpas_pool_get_array(nggps_diagnostics,'toplwupc',toplwupc)
 call mpas_pool_get_array(nggps_diagnostics,'topswup' ,topswup )
 call mpas_pool_get_array(nggps_diagnostics,'topswupc',topswupc)
 call mpas_pool_get_array(nggps_diagnostics,'topswdn' ,topswdn )
 call mpas_pool_get_array(nggps_diagnostics,'fluxr'   ,fluxr   )
 call mpas_pool_get_array(nggps_diagnostics,'dswcmp'  ,dswcmp  )
 call mpas_pool_get_array(nggps_diagnostics,'uswcmp'  ,uswcmp  )

!--- interface_fields:
 call mpas_pool_get_array(nggps_interface_fields,'sfclwup'  ,sfclwup  )
 call mpas_pool_get_array(nggps_interface_fields,'sfclwdn'  ,sfclwdn  )
 call mpas_pool_get_array(nggps_interface_fields,'sfclwupc' ,sfclwupc )
 call mpas_pool_get_array(nggps_interface_fields,'sfclwdnc' ,sfclwdnc )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswup'  ,sfcswup  )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswdn'  ,sfcswdn  )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswupc' ,sfcswupc )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswdnc' ,sfcswdnc )

 call mpas_pool_get_array(nggps_interface_fields,'sfcnirbmd',sfcnirbmd)
 call mpas_pool_get_array(nggps_interface_fields,'sfcnirdfd',sfcnirdfd)
 call mpas_pool_get_array(nggps_interface_fields,'sfcnirbmu',sfcnirbmu)
 call mpas_pool_get_array(nggps_interface_fields,'sfcnirdfu',sfcnirdfu)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisbmd',sfcvisbmd)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisdfd',sfcvisdfd)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisbmu',sfcvisbmu)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisdfu',sfcvisdfu)

!--- radiation_tendencies:
 call mpas_pool_get_array(nggps_radiation_tendencies,'coszen',coszen)
 call mpas_pool_get_array(nggps_radiation_tendencies,'coszdg',coszdg)
 call mpas_pool_get_array(nggps_radiation_tendencies,'semis' ,semis )
 call mpas_pool_get_array(nggps_radiation_tendencies,'sfalb' ,sfalb )
 call mpas_pool_get_array(nggps_radiation_tendencies,'tsflw' ,tsflw )
 call mpas_pool_get_array(nggps_radiation_tendencies,'htrlw' ,htrlw )
 call mpas_pool_get_array(nggps_radiation_tendencies,'htrsw' ,htrsw )

!... back to Registry:
 cldcov_l  => cld_properties%cldcov

 fluxr_l   => diags%fluxr
 dswcmp_l  => diags%dswcmp
 uswcmp_l  => diags%uswcmp
 topflw_l  => diags%topflw
 topfsw_l  => diags%topfsw

 sfcflw_l  => intr_fields%sfcflw
 sfcfsw_l  => intr_fields%sfcfsw 

 coszen    => rad_tendencies%coszen
 coszdg    => rad_tendencies%coszdg
 semis     => rad_tendencies%semis
 sfalb     => rad_tendencies%sfalb
 tsflw     => rad_tendencies%tsflw
 htrsw_l   => rad_tendencies%htrsw
 htrlw_l   => rad_tendencies%htrlw

 do iCell = 1, nCellsSolve
    sfclwup(iCell)  = sfcflw_l(iCell)%upfxc
    sfclwupc(iCell) = sfcflw_l(iCell)%upfx0
    sfclwdn(iCell)  = sfcflw_l(iCell)%dnfxc
    sfclwdnc(iCell) = sfcflw_l(iCell)%dnfx0
    sfcswup(iCell)  = sfcfsw_l(iCell)%upfxc
    sfcswupc(iCell) = sfcfsw_l(iCell)%upfx0
    sfcswdn(iCell)  = sfcfsw_l(iCell)%dnfxc
    sfcswdnc(iCell) = sfcfsw_l(iCell)%dnfx0
    toplwup(iCell)  = topflw_l(iCell)%upfxc
    toplwupc(iCell) = topflw_l(iCell)%upfx0
    topswup(iCell)  = topfsw_l(iCell)%upfxc
    topswupc(iCell) = topfsw_l(iCell)%upfx0
    topswdn(iCell)  = topfsw_l(iCell)%dnfxc
    do k = 1, 4
       dswcmp(k,iCell) = dswcmp_l(iCell,k)
       uswcmp(k,iCell) = uswcmp_l(iCell,k)
    enddo
    do k = 1, nFxr
       fluxr(k,iCell)  = fluxr_l(iCell,k)
    enddo 
    do k = 1, nVertLevels
       cldcov(k,iCell) = cldcov_l(iCell,k)
       htrsw(k,iCell)  = htrsw_l(iCell,k) 
       htrlw(k,iCell)  = htrlw_l(iCell,k) 
    enddo
 enddo

 do iCell = 1, nCellsSolve
    sfcnirbmd(iCell)=  dswcmp(1,iCell)
    sfcnirdfd(iCell) = dswcmp(2,iCell)
    sfcvisbmd(iCell) = dswcmp(3,iCell)
    sfcvisdfd(iCell) = dswcmp(4,iCell)
    sfcnirbmu(iCell) = uswcmp(1,iCell)
    sfcnirdfu(iCell) = uswcmp(2,iCell)
    sfcvisbmu(iCell) = uswcmp(3,iCell)
    sfcvisdfu(iCell) = uswcmp(4,iCell)
 enddo

!do iCell = 1,20
!   write(0,201) iCell,sfclwdn(iCell),sfcswdn(iCell),sfcswdn(iCell)-sfcswup(iCell),sfcnirbmd(iCell),          &
!      sfcnirdfd(iCell),sfcvisbmd(iCell),sfcvisdfd(iCell),sfcnirbmu(iCell),sfcnirdfu(iCell),sfcvisbmu(iCell), &
!      sfcvisdfu(iCell)
!enddo
 201 format(i7,12(1x,e15.8))

!write(0,*) '--- end subroutine nuopc_nggps_rad_out:'

 end subroutine nuopc_rad_write

!=================================================================================================================
 subroutine nuopc_phys_write(statein,stateout,diags,mesh,nggps_physics_tendencies)
!=================================================================================================================

!input arguments:
 type(diagnostics),intent(in)     :: diags
 type(state_fields_in),intent(in) :: statein
 type(state_fields_out),intent(in):: stateout
 type(mpas_pool_type),intent(in)  :: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_physics_tendencies

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

 real(kind=PHYSKIND),dimension(:,:),pointer:: tend_thradlw,tend_thradsw,tend_thdiffusion,tend_thcudeep, &
                                           tend_thcushallow,tend_thmicrop
 real(kind=PHYSKIND),dimension(:,:),pointer:: tend_qvdiffusion,tend_qvcudeep,tend_qvcushallow,tend_qvmicrop
 real(kind=PHYSKIND),dimension(:,:),pointer:: tend_qo3diffusion,tend_qo3prodct
 real(kind=PHYSKIND),dimension(:,:),pointer:: tend_udiffusion,tend_ucudeep,tend_ugwdc,tend_ugwdo
 real(kind=PHYSKIND),dimension(:,:),pointer:: tend_vdiffusion,tend_vcudeep,tend_vgwdc,tend_vgwdo

 real(kind=PHYSKIND),dimension(:,:),pointer:: tend_thphysics,tend_qvphysics,tend_qcphysics,tend_qo3physics, &
                                           tend_uphysics,tend_vphysics

 real(kind=kind_phys),dimension(:,:),pointer  :: prslk,tgrs,ugrs,vgrs
 real(kind=kind_phys),dimension(:,:),pointer  :: gt0,gu0,gv0
 real(kind=kind_phys),dimension(:,:,:),pointer:: qgrs,gq0
 real(kind=kind_phys),dimension(:,:,:),pointer:: dq3dt,dt3dt,du3dt,dv3dt

!local variables:
 integer:: iCell,k
 real(kind=PHYSKIND):: dti
 real(kind=PHYSKIND):: sphm1,sphp1

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine nuopc_write_out:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_physics_tendencies,'thradlw'     ,tend_thradlw     )
 call mpas_pool_get_array(nggps_physics_tendencies,'thradsw'     ,tend_thradsw     )
 call mpas_pool_get_array(nggps_physics_tendencies,'thdiffusion' ,tend_thdiffusion )
 call mpas_pool_get_array(nggps_physics_tendencies,'thcudeep'    ,tend_thcudeep    )
 call mpas_pool_get_array(nggps_physics_tendencies,'thcushallow' ,tend_thcushallow )
 call mpas_pool_get_array(nggps_physics_tendencies,'thmicrop'    ,tend_thmicrop    )

 call mpas_pool_get_array(nggps_physics_tendencies,'qvdiffusion' ,tend_qvdiffusion )
 call mpas_pool_get_array(nggps_physics_tendencies,'qvcudeep'    ,tend_qvcudeep    )
 call mpas_pool_get_array(nggps_physics_tendencies,'qvcushallow' ,tend_qvcushallow )
 call mpas_pool_get_array(nggps_physics_tendencies,'qvmicrop'    ,tend_qvmicrop    )

 call mpas_pool_get_array(nggps_physics_tendencies,'qo3diffusion',tend_qo3diffusion)
 call mpas_pool_get_array(nggps_physics_tendencies,'qo3prodct'   ,tend_qo3prodct   )

 call mpas_pool_get_array(nggps_physics_tendencies,'udiffusion'  ,tend_udiffusion  )
 call mpas_pool_get_array(nggps_physics_tendencies,'ucudeep'     ,tend_ucudeep     )
 call mpas_pool_get_array(nggps_physics_tendencies,'ugwdc'       ,tend_ugwdc       )
 call mpas_pool_get_array(nggps_physics_tendencies,'ugwdo'       ,tend_ugwdo       )

 call mpas_pool_get_array(nggps_physics_tendencies,'vdiffusion'  ,tend_vdiffusion  )
 call mpas_pool_get_array(nggps_physics_tendencies,'vcudeep'     ,tend_vcudeep     )
 call mpas_pool_get_array(nggps_physics_tendencies,'vgwdc'       ,tend_vgwdc       )
 call mpas_pool_get_array(nggps_physics_tendencies,'vgwdo'       ,tend_vgwdo       )

 call mpas_pool_get_array(nggps_physics_tendencies,'thphysics'   ,tend_thphysics   )
 call mpas_pool_get_array(nggps_physics_tendencies,'qvphysics'   ,tend_qvphysics   )
 call mpas_pool_get_array(nggps_physics_tendencies,'qcphysics'   ,tend_qcphysics   )
 call mpas_pool_get_array(nggps_physics_tendencies,'qo3physics'  ,tend_qo3physics  )
 call mpas_pool_get_array(nggps_physics_tendencies,'uphysics'    ,tend_uphysics    )
 call mpas_pool_get_array(nggps_physics_tendencies,'vphysics'    ,tend_vphysics    )

 prslk => statein%prslk
 tgrs  => statein%tgrs
 qgrs  => statein%qgrs
 ugrs  => statein%ugrs
 vgrs  => statein%vgrs

 gt0   => stateout%gt0
 gq0   => stateout%gq0
 gu0   => stateout%gu0
 gv0   => stateout%gv0

 dq3dt => diags%dq3dt
 dt3dt => diags%dt3dt
 du3dt => diags%du3dt
 dv3dt => diags%dv3dt

 dti = 1._PHYSKIND/dt_dyn
!potential temperature:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_thradlw(k,iCell)     = dt3dt(iCell,k,1)*dti/prslk(iCell,k)
       tend_thradsw(k,iCell)     = dt3dt(iCell,k,2)*dti/prslk(iCell,k)
       tend_thdiffusion(k,iCell) = dt3dt(iCell,k,3)*dti/prslk(iCell,k)
       tend_thcudeep(k,iCell)    = dt3dt(iCell,k,4)*dti/prslk(iCell,k)
       tend_thcushallow(k,iCell) = dt3dt(iCell,k,5)*dti/prslk(iCell,k)
       tend_thmicrop(k,iCell)    = dt3dt(iCell,k,6)*dti/prslk(iCell,k)
    enddo
 enddo

!water vapor and ozone mixing ratios:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_qvdiffusion(k,iCell)  = dq3dt(iCell,k,1)*dti
       tend_qvcudeep(k,iCell)     = dq3dt(iCell,k,2)*dti
       tend_qvcushallow(k,iCell)  = dq3dt(iCell,k,3)*dti
       tend_qvmicrop(k,iCell)     = dq3dt(iCell,k,4)*dti

       tend_qo3diffusion(k,iCell) = dq3dt(iCell,k,5)*dti
       tend_qo3prodct(k,iCell)    = dq3dt(iCell,k,6)*dti
    enddo
 enddo

!zonal and meridional winds:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_udiffusion(k,iCell) = du3dt(iCell,k,1)*dti
       tend_ugwdo(k,iCell)      = du3dt(iCell,k,2)*dti
       tend_ucudeep(k,iCell)    = du3dt(iCell,k,3)*dti
       tend_ugwdc(k,iCell)      = du3dt(iCell,k,4)*dti

       tend_vdiffusion(k,iCell) = dv3dt(iCell,k,1)*dti
       tend_vgwdo(k,iCell)      = dv3dt(iCell,k,2)*dti
       tend_vcudeep(k,iCell)    = dv3dt(iCell,k,3)*dti
       tend_vgwdc(k,iCell)      = dv3dt(iCell,k,4)*dti
    enddo
 enddo

!total tendencies using updated temperature,updated water vapor mixing ratio,updated total cloud
!condensate,updated ozone mixing ratio, and updated u- and v-momentum:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_thphysics(k,iCell)  = (gt0(iCell,k)-tgrs(iCell,k))*dti/prslk(iCell,k)
       tend_uphysics(k,iCell)   = (gu0(iCell,k)-ugrs(iCell,k))*dti
       tend_vphysics(k,iCell)   = (gv0(iCell,k)-vgrs(iCell,k))*dti
       !conversion from specific humidity to water vapor mixing ratio:
       tend_qvphysics(k,iCell) = gq0(iCell,k,1)/(1._PHYSKIND-gq0(iCell,k,1))
       tend_qvphysics(k,iCell) = (tend_qvphysics(k,iCell)-qgrs(iCell,k,1)/(1._PHYSKIND-qgrs(iCell,k,1)))*dti

       sphp1 = 1._PHYSKIND + gq0(iCell,k,1)
       sphm1 = 1._PHYSKIND + qgrs(iCell,k,1)
       !conversion from "specific" ozone/cloud condensate mixing ratios to mixing ratios:
       tend_qo3physics(k,iCell) = (gq0(iCell,k,2)*sphp1 - qgrs(iCell,k,2)*sphm1) * dti
       tend_qcphysics(k,iCell)  = (gq0(iCell,k,3)*sphp1 - qgrs(iCell,k,3)*sphm1) * dti

    enddo
 enddo

!write(0,*) 'max tend_uphysics   = ',maxval(tend_uphysics(:,1:nCellsSolve)  )
!write(0,*) 'min tend_uphysics   = ',minval(tend_uphysics(:,1:nCellsSolve)  )
!write(0,*) 'max tend_vphysics   = ',maxval(tend_vphysics(:,1:nCellsSolve)  )
!write(0,*) 'min tend_vphysics   = ',minval(tend_vphysics(:,1:nCellsSolve)  )
!write(0,*) 'max tend_thphysics  = ',maxval(tend_thphysics(:,1:nCellsSolve) )
!write(0,*) 'min tend_thphysics  = ',minval(tend_thphysics(:,1:nCellsSolve) )
!write(0,*) 'max tend_qvphysics  = ',maxval(tend_qvphysics(:,1:nCellsSolve) )
!write(0,*) 'min tend_qvphysics  = ',minval(tend_qvphysics(:,1:nCellsSolve) )
!write(0,*) 'max tend_qcphysics  = ',maxval(tend_qcphysics(:,1:nCellsSolve) )
!write(0,*) 'min tend_qcphysics  = ',minval(tend_qcphysics(:,1:nCellsSolve) )
!write(0,*) 'max tend_qo3physics = ',maxval(tend_qo3physics(:,1:nCellsSolve))
!write(0,*) 'min tend_qo3physics = ',minval(tend_qo3physics(:,1:nCellsSolve))

!write(0,*) '--- end subroutine nuopc_nggps_physics_write:'

 end subroutine nuopc_phys_write

!=================================================================================================================
 end module mpas_atmphys_nggps_driver_out
!=================================================================================================================

