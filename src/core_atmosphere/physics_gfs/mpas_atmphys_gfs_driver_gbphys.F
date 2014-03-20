!==================================================================================================
 module mpas_atmphys_gfs_driver_gbphys
 use mpas_dmpar
 use mpas_grid_types
 use mpas_kind_types

 use mpas_atmphys_gfs_interface
 use mpas_atmphys_manager, only: gmt,year,month,day,hour,utc_h,julday,curr_julday 
 use mpas_atmphys_gfs_vars

!GFS physics:
 use mersenne_twister
 use namelist_def

 implicit none
 private
 public:: allocate_gbphys,    &
          deallocate_gbphys,  &
          mpas_atmphys_gbphys


!>\brief MPAS top driver for the GFS physics, except the radiation codes.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-01-27.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Added the calculation of logical lscca for outputs of the diagnostics cv, cvb, cvt, acv,
!>      acvb, and acvt which are instantaneous and accumulated convection cloud fraction, and
!>      cloud base and cloud-top of convective clouds.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-19.
!>    * Replaced the variables xkzm_m, xkzm_h, and xkzm_s, with the actual background vertical
!>      diffusivity coefficients provided by the GFS physics.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04.20. 
!>    * Corrected the argument list in the call to subroutine gbphys for the water equivalent snow
!>      snow depth (snow), the actual snow depth (snowh), and the fractional area of snow (snowc).
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-20.
!>    * Removed inputs/outputs in the argument list to subroutine gbphys needed in subroutine
!>      sfc_nst since the variable nst_fcst is set to 0 and that subroutine is not called. Doing
!>      so helps reduce the size of the restart and output files.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-20.
!>    * Cleaned-up sourcecode for non-needed variables and diagnostics.
!>    * Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-20.


 contains

!==================================================================================================
 subroutine allocate_gbphys(mesh)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh

!--------------------------------------------------------------------------------------------------

 if(.not.allocated(dt3dt)  ) allocate(dt3dt(mesh%nCellsSolve,mesh%nVertLevels,6)            )
 if(.not.allocated(dq3dt)  ) allocate(dq3dt(mesh%nCellsSolve,mesh%nVertLevels,5+pl_coeff)   )
 if(.not.allocated(du3dt)  ) allocate(du3dt(mesh%nCellsSolve,mesh%nVertLevels,4)            )
 if(.not.allocated(dv3dt)  ) allocate(dv3dt(mesh%nCellsSolve,mesh%nVertLevels,4)            )

!allocation of arrays for miscellaneous:
 if(.not.allocated(fscav)  ) allocate(fscav(ntrac-ncld-1)                                   )

 if(.not.allocated(slmsk)  ) allocate(slmsk(mesh%nCellsSolve)                               )
 if(.not.allocated(nlons)  ) allocate(nlons(mesh%nCellsSolve)                               )
 if(.not.allocated(sinlat) ) allocate(sinlat(mesh%nCellsSolve)                              )
 if(.not.allocated(coslat) ) allocate(coslat(mesh%nCellsSolve)                              )
 if(.not.allocated(rann)   ) allocate(rann(mesh%nCellsSolve,nrcm)                           )
 if(.not.allocated(rnp)    ) allocate(rnp(mesh%nCellsSolve,mesh%nVertLevels)                )

!allocation of arrays for ozone parameterization:
 if(.not.allocated(poz)    ) allocate(poz(levozp)                                           )
 if(.not.allocated(prdout) ) &
     allocate(prdout(1:mesh%nCellsSolve,1:mesh%nOznLevels,1:mesh%nOznPlCoeff))

!allocation of arrays for parameterization of vertical diffusion (gocart):
 if(.not.allocated(dkh)    ) allocate(dkh(mesh%nCellsSolve,mesh%nVertLevels)                )

!allocation of arrays for parameterized convection:
 if(.not.allocated(dpshc)  ) allocate(dpshc(1:mesh%nCellsSolve)                             )
 if(.not.allocated(upd_mf) ) allocate(upd_mf(mesh%nCellsSolve,mesh%nVertLevels)             )
 if(.not.allocated(dwn_mf) ) allocate(dwn_mf(mesh%nCellsSolve,mesh%nVertLevels)             )
 if(.not.allocated(det_mf) ) allocate(det_mf(mesh%nCellsSolve,mesh%nVertLevels)             )

!allocation of arrays for 2D and 3D restarts:
 if(.not.allocated(phy_f2d)) allocate(phy_f2d(mesh%nCellsSolve,mesh%numP2d)                 )
 if(.not.allocated(phy_f3d)) allocate(phy_f3d(mesh%nCellsSolve,mesh%nVertLevels,mesh%numP3d))

!allocation of arrays for parameterization of long- and short-wave radiation:
 if(.not.allocated(htrlw)  ) allocate(htrlw(mesh%nCellsSolve,mesh%nVertLevels)              )
 if(.not.allocated(htrsw)  ) allocate(htrsw(mesh%nCellsSolve,mesh%nVertLevels)              )

!allocation of arrays for parameterization of gravity wave drag over orography:
 if(.not.allocated(hprime) ) allocate(hprime(mesh%nCellsSolve,mesh%nMtvr)                   )

!allocation of arrays for parameterization of land surface scheme:
 if(.not.allocated(slc)    ) allocate(slc(mesh%nCellsSolve,mesh%nSoilLevels)                )
 if(.not.allocated(smc)    ) allocate(smc(mesh%nCellsSolve,mesh%nSoilLevels)                )
 if(.not.allocated(stc)    ) allocate(stc(mesh%nCellsSolve,mesh%nSoilLevels)                )

 end subroutine allocate_gbphys

!==================================================================================================
 subroutine deallocate_gbphys
!==================================================================================================

 if(allocated(dt3dt)  ) deallocate(dt3dt  )
 if(allocated(dq3dt)  ) deallocate(dq3dt  )
 if(allocated(du3dt)  ) deallocate(du3dt  )
 if(allocated(dv3dt)  ) deallocate(dv3dt  )

!deallocation of arrays for miscellaneous:
 if(allocated(fscav)  ) deallocate(fscav  )
 if(allocated(slmsk)  ) deallocate(slmsk  )
 if(allocated(nlons)  ) deallocate(nlons  )
 if(allocated(sinlat) ) deallocate(sinlat )
 if(allocated(coslat) ) deallocate(coslat )
 if(allocated(rann)   ) deallocate(rann   )
 if(allocated(rnp)    ) deallocate(rnp    )

!deallocation of arrays for ozone parameterization:
 if(allocated(poz)    ) deallocate(poz    )
 if(allocated(prdout) ) deallocate(prdout )

!deallocation of arrays for parameterization of vertical diffusion (gocart):
 if(allocated(dkh)    ) deallocate(dkh    )

!deallocation of arrays for parameterized convection:
 if(allocated(dpshc)  ) deallocate(dpshc  )
 if(allocated(upd_mf) ) deallocate(upd_mf )
 if(allocated(dwn_mf) ) deallocate(dwn_mf )
 if(allocated(det_mf) ) deallocate(det_mf )

!deallocation of arrays for 2D and 3D restarts:
 if(allocated(phy_f2d)) deallocate(phy_f2d)
 if(allocated(phy_f3d)) deallocate(phy_f3d)

!deallocation of arrays for parameterization of long- and short-wave radiation:
 if(allocated(htrlw)  ) deallocate(htrlw  )
 if(allocated(htrsw)  ) deallocate(htrsw  )

!deallocation of arrays for parameterization of gravity wave drag over orography:
 if(allocated(hprime) ) deallocate(hprime )

!deallocation of arrays for parameterization of land surface scheme:
 if(allocated(slc)    ) deallocate(slc    )
 if(allocated(smc)    ) deallocate(smc    )
 if(allocated(stc)    ) deallocate(stc    )  

 end subroutine deallocate_gbphys

!==================================================================================================
 subroutine gbphys_from_MPAS(mesh,atm_input,sfc_input,diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(atm_input_type),intent(in):: atm_input

!inout arguments:
 type(sfc_input_type),intent(inout):: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables:
 integer:: iCell,k,n

!--------------------------------------------------------------------------------------------------

 do iCell = 1, mesh%nCellsSolve
    slmsk(iCell) = sfc_input%landmask%array(iCell)
    if(sfc_input%xice%array(iCell) .gt. 0._RKIND) then
       slmsk(iCell) = 2._RKIND
       diag_physics%zorl%array(iCell) = z0_sib(13)
    endif

!temporary:
    sfc_input%slope%array(iCell) = min(sfc_input%slope%array(iCell),9._RKIND)

    nlons(iCell)  = mesh%nCellsSolve
    sinlat(iCell) = sin(mesh%latCell%array(iCell))
    coslat(iCell) = cos(mesh%latCell%array(iCell))
 enddo

!initialization of local physics tendencies:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       dt3dt(iCell,k,:) = 0._RKIND
       dq3dt(iCell,k,:) = 0._RKIND
       du3dt(iCell,k,:) = 0._RKIND
       dv3dt(iCell,k,:) = 0._RKIND
    enddo
 enddo

 do iCell = 1, mesh%nCellsSolve
    !initialization of 2D arrays for restarts:
    do n = 1, mesh%numP2d
       phy_f2d(iCell,n) = diag_physics % phy_f2d % array(n,iCell)
    enddo
    do k = 1, mesh%nVertLevels
       !initialization of 3D arrays for restarts:
       do n = 1, mesh%numP3d
          phy_f3d(iCell,k,n) = diag_physics % phy_f3d % array(n,k,iCell)
       enddo
       !initialization of in cloud precipitation rate:
       rnp(iCell,k) = diag_physics % rnp % array(k,iCell)
       !initialization of arrays for parameterization of vertical diffusion (gocart):
       dkh(iCell,k) = diag_physics % dkh % array(k,iCell)
    enddo
 enddo

!initialization of the maximum pressure depth for shallow convection, others...:
 write(0,*) '--- gen_coord_hybrid =', gen_coord_hybrid
 write(0,*) '--- thermodyn_id     =', thermodyn_id
 do iCell = 1, mesh%nCellsSolve
    phil(iCell,mesh%nVertLevels) = 0._RKIND !forces the calculation of geopotential in gbphys.
    dpshc(iCell) = 0.3 * pgr_hyd(iCell)

    if(gen_coord_hybrid .and. thermodyn_id == 3) then
       prslk(iCell,1) = 0._RKIND !forces the calculation of geopotential in gbphys.
       prsik(iCell,1) = 0._RKIND !forces the calculation of geopotential in gbphys.
    endif
 enddo

!initialization for parameterization of gravity wave drag over orography:
 do iCell = 1, mesh%nCellsSolve
    do n = 1, mesh%nMtvr
       hprime(iCell,n)  = sfc_input % orogwd % array(n,iCell)
    enddo
    diag_physics % dugwd % array(iCell) = 0._RKIND
    diag_physics % dvgwd % array(iCell) = 0._RKIND
 enddo

!initialization for parameterized convection: CHECK IF INITIALIZATION IS NEEDED FOR ACV,ACVB,ACVT ??
 do iCell = 1, mesh%nCellsSolve
    diag_physics % acv  % array(iCell) = 0._RKIND
    diag_physics % acvb % array(iCell) = 0._RKIND
    diag_physics % acvt % array(iCell) = 0._RKIND
    do k = 1, mesh%nVertLevels
       upd_mf(iCell,k) = 0._RKIND
       dwn_mf(iCell,k) = 0._RKIND
       det_mf(iCell,k) = 0._RKIND
    enddo
 enddo

!initialization for parameterization of long-wave and short-wave radiation:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       htrlw(iCell,k) = tend_physics % rthratenlw % array(k,iCell)
       htrsw(iCell,k) = tend_physics % rthratensw % array(k,iCell)
    enddo
 enddo

!initialization of local liquid soil moisture, soil moisture, and soil temperature for land
!surface scheme:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nSoilLevels
       slc(iCell,k) = sfc_input % smcrel % array(k,iCell)
       smc(iCell,k) = sfc_input % smois  % array(k,iCell)
       stc(iCell,k) = sfc_input % tslb   % array(k,iCell)
    enddo
 enddo

!initialization for ozone:
 do k = 1, mesh%nOznLevels
    poz(k) = atm_input%pl_pres%array(k)
 enddo
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nOznLevels
       do n = 1, mesh%nOznPlCoeff
          prdout(iCell,k,n) = diag_physics%o3clim%array(k,n,iCell)
       enddo
    enddo
 enddo

!initialization of temperature, momentum, water vapor, water/ice condensates, and ozone mixing
!ratios updated in gbphys:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       gt0(iCell,k)   = 0._RKIND
       gu0(iCell,k)   = 0._RKIND
       gv0(iCell,k)   = 0._RKIND
       gq0(iCell,k,:) = 0._RKIND       
    enddo
 enddo

!outputs to subroutine to gbphys:
 do iCell = 1, mesh%nCellsSolve
    diag_physics%dlwsfci%array(iCell) = 0._RKIND
    diag_physics%ulwsfci%array(iCell) = 0._RKIND
    diag_physics%dswsfci%array(iCell) = 0._RKIND
    diag_physics%uswsfci%array(iCell) = 0._RKIND
    diag_physics%dqsfci%array(iCell)  = 0._RKIND
    diag_physics%dtsfci%array(iCell)  = 0._RKIND
    diag_physics%epi%array(iCell)     = 0._RKIND
    diag_physics%gfluxi%array(iCell)  = 0._RKIND
 enddo

 end subroutine gbphys_from_MPAS

!==================================================================================================
 subroutine gbphys_to_MPAS(dminfo,mesh,diag_physics,tend_physics,sfc_input)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics
 type(sfc_input_type),intent(inout):: sfc_input

!local variables:
 integer:: iCell,k,n
 real(kind=RKIND):: dti

 real(kind=RKIND):: tend_min,tend_max
 real(kind=RKIND):: global_tend_min,global_tend_max

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine gbphys_to_MPAS:'

!fill MPAS physics tendencies:
 dti = 1._RKIND/dt_dyn
!potential temperature:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_physics % thradlw     % array(k,iCell) = dt3dt(iCell,k,1)*dti/prslk(iCell,k)
       tend_physics % thradsw     % array(k,iCell) = dt3dt(iCell,k,2)*dti/prslk(iCell,k)
       tend_physics % thdiffusion % array(k,iCell) = dt3dt(iCell,k,3)*dti/prslk(iCell,k)
       tend_physics % thcudeep    % array(k,iCell) = dt3dt(iCell,k,4)*dti/prslk(iCell,k)
       tend_physics % thcushallow % array(k,iCell) = dt3dt(iCell,k,5)*dti/prslk(iCell,k)
       tend_physics % thmicrop    % array(k,iCell) = dt3dt(iCell,k,6)*dti/prslk(iCell,k)
    enddo
 enddo

!water vapor and ozone mixing ratios:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_physics % qvdiffusion  % array(k,iCell) = dq3dt(iCell,k,1)*dti
       tend_physics % qvcudeep     % array(k,iCell) = dq3dt(iCell,k,2)*dti
       tend_physics % qvcushallow  % array(k,iCell) = dq3dt(iCell,k,3)*dti
       tend_physics % qvmicrop     % array(k,iCell) = dq3dt(iCell,k,4)*dti

       tend_physics % qo3diffusion % array(k,iCell) = dq3dt(iCell,k,5)*dti
       tend_physics % qo3prodct    % array(k,iCell) = dq3dt(iCell,k,6)*dti
    enddo
 enddo
 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_min = min(tend_min,tend_physics%qvdiffusion%array(k,iCell))
       tend_max = max(tend_max,tend_physics%qvdiffusion%array(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvdiffusion',global_tend_min,global_tend_max

 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_min = min(tend_min,tend_physics%qvcudeep%array(k,iCell))
       tend_max = max(tend_max,tend_physics%qvcudeep%array(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvcudeep',global_tend_min,global_tend_max

 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_min = min(tend_min,tend_physics%qvcushallow%array(k,iCell))
       tend_max = max(tend_max,tend_physics%qvcushallow%array(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvcushallow',global_tend_min,global_tend_max

 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_min = min(tend_min,tend_physics%qvmicrop%array(k,iCell))
       tend_max = max(tend_max,tend_physics%qvmicrop%array(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvmicrop',global_tend_min,global_tend_max

!zonal and meridional winds:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_physics % udiffusion % array(k,iCell) = du3dt(iCell,k,1)*dti
       tend_physics % ugwdo      % array(k,iCell) = du3dt(iCell,k,2)*dti 
       tend_physics % ucudeep    % array(k,iCell) = du3dt(iCell,k,3)*dti
       tend_physics % ugwdc      % array(k,iCell) = du3dt(iCell,k,4)*dti 

       tend_physics % vdiffusion % array(k,iCell) = dv3dt(iCell,k,1)*dti
       tend_physics % vgwdo      % array(k,iCell) = dv3dt(iCell,k,2)*dti
       tend_physics % vcudeep    % array(k,iCell) = dv3dt(iCell,k,3)*dti 
       tend_physics % vgwdc      % array(k,iCell) = dv3dt(iCell,k,4)*dti
    enddo
 enddo

!total tendencies using updated temperature,updated water vapor mixing ratio,updated total cloud
!condensate,updated ozone mixing ratio, and updated u- and v-momentum:
 write(0,*) "max tend thradlw     = ",maxval(tend_physics%thradlw%array(:,1:mesh%nCellsSolve)    )
 write(0,*) "min tend thradlw     = ",minval(tend_physics%thradlw%array(:,1:mesh%nCellsSolve)    )
 write(0,*) "max tend thradsw     = ",maxval(tend_physics%thradsw%array(:,1:mesh%nCellsSolve)    )
 write(0,*) "min tend thradsw     = ",minval(tend_physics%thradsw%array(:,1:mesh%nCellsSolve)    )
 write(0,*) "max tend thdiffusion = ",maxval(tend_physics%thdiffusion%array(:,1:mesh%nCellsSolve))
 write(0,*) "min tend thdiffusion = ",minval(tend_physics%thdiffusion%array(:,1:mesh%nCellsSolve))
 write(0,*) "max tend thcudeep    = ",maxval(tend_physics%thcudeep%array(:,1:mesh%nCellsSolve)   )
 write(0,*) "min tend thcudeep    = ",minval(tend_physics%thcudeep%array(:,1:mesh%nCellsSolve)   )
 write(0,*) "max tend thcushallow = ",maxval(tend_physics%thcushallow%array(:,1:mesh%nCellsSolve))
 write(0,*) "min tend thcushallow = ",minval(tend_physics%thcushallow%array(:,1:mesh%nCellsSolve))
 write(0,*) "max tend thmicrop    = ",maxval(tend_physics%thmicrop%array(:,1:mesh%nCellsSolve)   )
 write(0,*) "min tend thmicrop    = ",minval(tend_physics%thmicrop%array(:,1:mesh%nCellsSolve)   )


 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_physics % thphysics  % array(k,iCell) = (gt0(iCell,k)-tgrs(iCell,k))*dti/prslk(iCell,k)
       tend_physics % qvphysics  % array(k,iCell) = (gq0(iCell,k,1)-qgrs(iCell,k,1))*dti 
       tend_physics % qo3physics % array(k,iCell) = (gq0(iCell,k,2)-qgrs(iCell,k,2))*dti
       tend_physics % qcphysics  % array(k,iCell) = (gq0(iCell,k,3)-qgrs(iCell,k,3))*dti
       tend_physics % uphysics   % array(k,iCell) = (gu0(iCell,k)-ugrs(iCell,k))*dti
       tend_physics % vphysics   % array(k,iCell) = (gv0(iCell,k)-vgrs(iCell,k))*dti
    enddo
 enddo
 write(0,*) "max tend thphysics   = ",maxval(tend_physics%thphysics%array(:,1:mesh%nCellsSolve)  )
 write(0,*) "min tend thphysics   = ",minval(tend_physics%thphysics%array(:,1:mesh%nCellsSolve)  )
 write(0,*)

 do iCell = 1, mesh%nCellsSolve
    !update 2D arrays for restarts:
    do n = 1, mesh%numP2d
       diag_physics % phy_f2d % array(n,iCell) = phy_f2d(iCell,n)
    enddo
    do k = 1, mesh%nVertLevels
       !update 3D arrays for restarts:
       do n = 1, mesh%numP3d
          diag_physics % phy_f3d % array(n,k,iCell) = phy_f3d(iCell,k,n)
       enddo
       !update in-cloud precipitation rate:
       diag_physics % rnp % array(k,iCell) = rnp(iCell,k)
       !update vertical diffusion coefficients (gocart):
       diag_physics % dkh % array(k,iCell) = dkh(iCell,k)
    enddo
 enddo

!fill MPAS arrays from parameterized deep convection:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       diag_physics % upd_mf % array(k,iCell) = upd_mf(iCell,k)
       diag_physics % dwn_mf % array(k,iCell) = dwn_mf(iCell,k)
       diag_physics % det_mf % array(k,iCell) = det_mf(iCell,k)
    enddo
 enddo

!update MPAS arrays from land surface scheme:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nSoilLevels
       sfc_input % smcrel % array(k,iCell) = slc(iCell,k)
       sfc_input % smois  % array(k,iCell) = smc(iCell,k)
       sfc_input % tslb   % array(k,iCell) = stc(iCell,k)
    enddo
 enddo

 write(0,*) '--- end subroutine gbphys_to_MPAS'

 end subroutine gbphys_to_MPAS

!==================================================================================================
 subroutine mpas_atmphys_gbphys(itimestep,xtime_s,dminfo,mesh,atm_input,sfc_input, &
                                diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh

 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(atm_input_type),intent(inout):: atm_input
 type(sfc_input_type),intent(inout):: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables:
 logical:: flipv

 integer:: iCell,k,n
 integer:: nCellsSolve,nVertLevels,nSoilLevels
 integer:: numP2d,numP3d
 integer:: me
 integer:: lat,lonf,latg,kdt,ko3

 real(kind=RKIND),parameter:: cons_24 = 24._RKIND
 real(kind=RKIND),parameter:: cons_99 = 99._RKIND

 real(kind=RKIND):: xtime_h
 real(kind=RKIND):: dtp,dtf,solhr
 real(kind=RKIND):: clstp


!--------------------------------------------------------------------------------------------------
 write(0,100) itimestep

 !formats:
 100 format(/,' --- enter subroutine mpas_atmphys_gbphys ',i6)
 101 format('     ITIMESTEP   =', i6   ,/, &
            '     YEAR        =', i6   ,/, &
            '     MONTH       =', i6   ,/, &
            '     DAY         =', i6   ,/, &
            '     HOUR        =', i6   ,/, &
            '     JULDAY      =', i6   ,/, &
            '     CURR_JULDAY =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     XTIME_H     =', f16.9,/)
 xtime_h = xtime_s/3600._RKIND
 write(0,101) itimestep,year,month,day,hour,julday,curr_julday,utc_h,xtime_h

!... copy all MPAS arrays to GFS arrays:
 call gbphys_from_MPAS(mesh,atm_input,sfc_input,diag_physics,tend_physics)

!abort MPAS if the new SAS convection scheme is not used. If RAS or old SAS are used, we need
!to compute a random number seed (see additional sourcecode in gloopb.f):
 if(.not.newsas) then
    if(random_clds) then
       write(0,*) 'Error in choice of convection scheme, set newsas to true for now.'
       call mpas_dmpar_abort(dminfo)
    endif
 endif

!... date conversion:
 fhour = utc_h
 solhr = mod(fhour+idate(1),cons_24)

!... initialization of variable clstp needed in the calculation of bottom and top pressures for
!convective clouds (Ken Campana diagnostics). set switch for saving convective clouds:
 lscca = mod(itimestep,nsswr) == 0
 if(lscca.and.lsswr) then
    clstp = 1100.+min(fhswr,fhour,cons_99)  !initialize,accumulate,convert
 elseif(lscca) then
    clstp = 0100.+min(fhswr,fhour,cons_99)  !accumulate,convert
 elseif(lsswr) then
    clstp = 1100.                           !initialize,accumulate
 else
    clstp = 0100.                           !accumulate
 endif
 lslwr = .false.
 lsswr = .false.

!... initialization of local variables:
 nCellsSolve = mesh % nCellsSolve
 nVertLevels = mesh % nVertLevels
 nSoilLevels = mesh % nSoilLevels
 numP2d = mesh % numP2d
 numP3d = mesh % numP3d
!if(dminfo % my_proc_id == IO_NODE) me = 0
 me = dminfo % my_proc_id

 lonf = nCellsSolve
 lat  = 1
 latg = 1

 kdt  = itimestep
 ko3  = levozp
 dtp  = dt_dyn
 dtf  = dt_dyn

 flipv = .false.
 if(iflip == 1) flipv = .true.

!write(0,*) 'nCellsSolve = ', nCellsSolve
!write(0,*) 'nVertLevels = ', nVertLevels
!write(0,*) 'nSoilLevels = ', nSoilLevels
!write(0,*) 'lsm         = ', lsm         !
!write(0,*) 'ntrac       = ', ntrac       ! number of tracers, including ozone.
!write(0,*) 'ncld        = ', ncld        ! number of cloud condensates.
!write(0,*) 'ntoz        = ', ntoz        ! index number for ozone in ntrac.
!write(0,*) 'ntcw        = ', ntcw        ! index number for first cloud condenstates in ntrac.
!write(0,*)
!write(0,*) 'nmtvr       = ', nmtvr       ! number of input topography statistics.
!write(0,*) 'nrcm        = ', nrcm        ! 2nd dimension for random number rann.
!write(0,*) 'ko3         = ', ko3         ! number of ozone layers
!write(0,*) 'lonf        = ', lonf
!write(0,*) 'latg        = ', latg
!write(0,*) 'jcap        = ', jcap
!write(0,*) 'num_p3d     = ', num_p3d
!write(0,*) 'num_p2d     = ', num_p2d
!write(0,*)
!write(0,*) 'kdt         = ', kdt
!write(0,*) 'lat         = ', lat
!write(0,*) 'me          = ', me
!write(0,*) 'pl_coeff    = ', pl_coeff

!write(0,*) 'nlons       = ', nlons(1)
!write(0,*) 'ncw         = ', ncw
!write(0,*) 'flgmin      = ', flgmin
!write(0,*) 'crtrh       = ', crtrh
!write(0,*) 'cdmbgwd     = ', cdmbgwd
!write(0,*)
!write(0,*) 'ccwf        = ', ccwf
!write(0,*) 'dlqf        = ', dlqf
!write(0,*) 'ctei_rm     = ', ctei_rm
!write(0,*) 'clstp       = ', clstp
!write(0,*) 'dtp         = ', dtp
!write(0,*) 'fhour       = ', fhour
!write(0,*) 'solhr       = ', solhr
!write(0,*)
!write(0,*) 'ras              = ',ras
!write(0,*) 'pre_rad          = ',pre_rad
!write(0,*) 'ldiag3d          = ',ldiag3d
!write(0,*) 'lggfs3d          = ',lggfs3d
!write(0,*) 'lssav            = ',lssav
!write(0,*) 'lssav_cc         = ',lssav_cc
!write(0,*) 'psautco          = ',psautco
!write(0,*) 'prautco          = ',prautco
!write(0,*) 'evpco            = ',evpco
!write(0,*) 'wminco           = ',wminco
!write(0,*) 'flipv            = ',flipv
!write(0,*) 'old_monin        = ',old_monin
!write(0,*) 'cnvgwd           = ',cnvgwd
!write(0,*) 'shal_cnv         = ',shal_cnv
!write(0,*) 'sashal           = ',sashal
!write(0,*) 'newsas           = ',newsas
!write(0,*) 'cal_pre          =',cal_pre
!write(0,*) 'mom4ice          =',mom4ice
!write(0,*) 'mstrat           =',mstrat
!write(0,*) 'trans_trac       =',trans_trac
!write(0,*) 'nst_fcst         =',nst_fcst
!write(0,*) 'moist_adj        =',moist_adj
!write(0,*) 'fscav            =',fscav
!write(0,*) 'thermodyn_id     =',thermodyn_id
!write(0,*) 'sfcpress_id      =',sfcpress_id
!write(0,*) 'gen_coord_hybrid =',gen_coord_hybrid
!stop

 201 format(i9,10(1x,e15.8))
 202 format(2i9,10(1x,e15.8))

!... CALL GFS PHYSICS DRIVER FOR ALL PHYSICS, EXCEPT RADIATION:
 call gbphys( &
    !... inputs:
    nCellsSolve, nCellsSolve, nVertLevels, nSoilLevels, lsm, ntrac, ncld, ntoz, ntcw, nmtvr,    & 
    nrcm, ko3, lonf, latg, jcap, num_p3d, num_p2d, kdt, lat, me, pl_coeff, nlons, ncw, flgmin,  & 
    crtrh, cdmbgwd, ccwf, dlqf, ctei_rm, clstp, dtp, dtf, fhour, solhr, slag, sdec, cdec,       & 
    sinlat, coslat, pgr, ugrs, vgrs, tgrs, qgrs, vvel, prsi, prsl, prslk, prsik, phii, phil,    &
    rann, prdout, poz, dpshc, hprime,                                                           &
    mesh%lonCell%array(1:nCellsSolve)           , mesh%latCell%array(1:nCellsSolve),            &
    sfc_input%slope%array(1:nCellsSolve)        , sfc_input%shdmin%array(1:nCellsSolve),        &
    sfc_input%shdmax%array(1:nCellsSolve)       , sfc_input%snoalb%array(1:nCellsSolve),        &
    sfc_input%tmn%array(1:nCellsSolve),  slmsk  , diag_physics%vegfra%array(1:nCellsSolve),     & 
    DBLE(sfc_input%ivgtyp%array(1:nCellsSolve)) , DBLE(sfc_input%isltyp%array(1:nCellsSolve)),  &
    diag_physics%uustar%array(1:nCellsSolve)    , sfc_input%ter%array(1:nCellsSolve),           &
    diag_physics%coszen%array(1:nCellsSolve)    , diag_physics%sfcdsw%array(1:nCellsSolve),     &
    diag_physics%sfcnsw%array(1:nCellsSolve)    , diag_physics%sfcdlw% array(1:nCellsSolve),    &
    diag_physics%tsflw%array(1:nCellsSolve)     , diag_physics%sfcemis%array(1:nCellsSolve),    &
    diag_physics%sfcalb%array(1:nCellsSolve)    ,                                               &
    htrsw, htrlw, ras, pre_rad, ldiag3d, lggfs3d, lssav, lssav_cc, bkgd_vdif_m, bkgd_vdif_h,    &
    bkgd_vdif_s, psautco, prautco, evpco, wminco, flipv, old_monin, cnvgwd, shal_cnv, sashal,   &
    newsas, cal_pre, mom4ice, mstrat, trans_trac, nst_fcst, moist_adj, fscav, thermodyn_id,     & 
    sfcpress_id, gen_coord_hybrid,                                                              &
    !... inputs/outputs:
    diag_physics%hice%array(1:nCellsSolve)    , sfc_input%xice%array(1:nCellsSolve),            &
    diag_physics%tisfc%array(1:nCellsSolve)   , sfc_input%skintemp%array(1:nCellsSolve),        &
    diag_physics%tprcp%array(1:nCellsSolve)   , diag_physics%cv%array(1:nCellsSolve),           & 
    diag_physics%cvb%array(1:nCellsSolve)     , diag_physics%cvt%array(1:nCellsSolve),          &
    diag_physics%srflag%array(1:nCellsSolve)  , sfc_input%snowh%array(1:nCellsSolve),           &
    sfc_input%snow%array(1:nCellsSolve)       , sfc_input%snowc%array(1:nCellsSolve),           &
    diag_physics%zorl%array(1:nCellsSolve)    , diag_physics%canopy%array(1:nCellsSolve),       &
    diag_physics%ffmm%array(1:nCellsSolve)    , diag_physics%ffhh%array(1:nCellsSolve),         &
    diag_physics%f10m%array(1:nCellsSolve)    , diag_physics%srunoff%array(1:nCellsSolve),      &
    diag_physics%evbsa%array(1:nCellsSolve)   , diag_physics%evcwa%array(1:nCellsSolve),        &
    diag_physics%snohfa%array(1:nCellsSolve)  , diag_physics%transa%array(1:nCellsSolve),       &
    diag_physics%sbsnoa%array(1:nCellsSolve)  , diag_physics%snowca%array(1:nCellsSolve),       &
    diag_physics%soilm%array(1:nCellsSolve)   , diag_physics%tmpmin%array(1:nCellsSolve),       &
    diag_physics%tmpmax%array(1:nCellsSolve)  , diag_physics%dusfc%array(1:nCellsSolve),        &
    diag_physics%dvsfc%array(1:nCellsSolve)   , diag_physics%dtsfc%array(1:nCellsSolve),        &
    diag_physics%dqsfc%array(1:nCellsSolve)   , diag_physics%totprcp%array(1:nCellsSolve),      &
    diag_physics%gflux%array(1:nCellsSolve)   , diag_physics%dlwsfc%array(1:nCellsSolve),       &
    diag_physics%ulwsfc%array(1:nCellsSolve)  , diag_physics%suntim%array(1:nCellsSolve),       & 
    diag_physics%runoff%array(1:nCellsSolve)  , diag_physics%ep%array(1:nCellsSolve),           &
    diag_physics%cldwrk%array(1:nCellsSolve)  , diag_physics%dugwd%array(1:nCellsSolve),        &
    diag_physics%dvgwd%array(1:nCellsSolve)   , diag_physics%psmean%array(1:nCellsSolve),       &    
    diag_physics%cnvprcp%array(1:nCellsSolve) , diag_physics%spfhmin%array(1:nCellsSolve),      &
    diag_physics%spfhmax%array(1:nCellsSolve) , dt3dt , dq3dt , du3dt , dv3dt ,                 &
    diag_physics%acv%array(1:nCellsSolve)     , diag_physics%acvb%array(1:nCellsSolve),         &
    diag_physics%acvt%array(1:nCellsSolve)    ,                                                 &
    slc, smc, stc, upd_mf, dwn_mf, det_mf, dkh, rnp, phy_f3d, phy_f2d,                          &
    dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys,                     &
    dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys,                     &
    dumphys, dumphys, dumphys, dumphys,                                                         &
    !... outputs:
    gt0, gq0, gu0, gv0,                                                                         &
    diag_physics%t2m%array(1:nCellsSolve)     , diag_physics%q2m%array(1:nCellsSolve),          &
    diag_physics%u10m%array(1:nCellsSolve)    , diag_physics%v10m%array(1:nCellsSolve),         &
    dumphys, dumphys,                                                                           &
    diag_physics%hpbl%array(1:nCellsSolve)    , diag_physics%pwat%array(1:nCellsSolve),         &
    diag_physics%t1%array(1:nCellsSolve)      , diag_physics%q1%array(1:nCellsSolve),           &
    diag_physics%u1%array(1:nCellsSolve)      , diag_physics%v1%array(1:nCellsSolve),           &
    diag_physics%chh%array(1:nCellsSolve)     , diag_physics%cmm%array(1:nCellsSolve),          &
    diag_physics%dlwsfci%array(1:nCellsSolve) , diag_physics%ulwsfci%array(1:nCellsSolve),      &
    diag_physics%dswsfci%array(1:nCellsSolve) , diag_physics%uswsfci%array(1:nCellsSolve),      & 
    diag_physics%dtsfci%array(1:nCellsSolve)  , diag_physics%dqsfci%array(1:nCellsSolve),       &
    diag_physics%gfluxi%array(1:nCellsSolve)  , diag_physics%epi%array(1:nCellsSolve),          & 
    diag_physics%smcwlt2%array(1:nCellsSolve) , diag_physics%smcref2%array(1:nCellsSolve),      &
    diag_physics%gsoil%array(1:nCellsSolve)   , diag_physics%gtmp2m%array(1:nCellsSolve),       &
    diag_physics%gustar%array(1:nCellsSolve)  , diag_physics%gpblh%array(1:nCellsSolve),        &
    diag_physics%gu10m%array(1:nCellsSolve)   , diag_physics%gv10m%array(1:nCellsSolve),        &
    diag_physics%gzorl%array(1:nCellsSolve)   , diag_physics%goro%array(1:nCellsSolve),         &
    dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys,                     &
    dumphys, dumphys, dumphys,dumphys                                                           &
            )
 write(0,*) '--- end subroutine gbphys:'

!... copy all GFS arrays to MPAS arrays:
 call gbphys_to_MPAS(dminfo,mesh,diag_physics,tend_physics,sfc_input)

 write(0,*) '--- end subroutine mpas_atmphys_gbphys:'

 end subroutine mpas_atmphys_gbphys

!==================================================================================================
 end module mpas_atmphys_gfs_driver_gbphys
!==================================================================================================
