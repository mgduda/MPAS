!=================================================================================================================
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
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine allocate_gbphys(mesh)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nSoilLevels
 integer,pointer:: nOznLevels,nOznPlCoeff
 integer,pointer:: numP2d,numP3d
 integer,pointer:: nGwd

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)
 call mpas_pool_get_dimension(mesh,'numP2d'     ,numP2d     )
 call mpas_pool_get_dimension(mesh,'numP3d'     ,numP3d     )
 call mpas_pool_get_dimension(mesh,'nGwd'       ,nGwd       )

 if(.not.allocated(dt3dt)    ) allocate(dt3dt(nCellsSolve,nVertLevels,6)         )
 if(.not.allocated(dq3dt)    ) allocate(dq3dt(nCellsSolve,nVertLevels,5+pl_coeff))
 if(.not.allocated(du3dt)    ) allocate(du3dt(nCellsSolve,nVertLevels,4)         )
 if(.not.allocated(dv3dt)    ) allocate(dv3dt(nCellsSolve,nVertLevels,4)         )

!allocation of arrays for miscellaneous:
 if(.not.allocated(fscav)    ) allocate(fscav(ntrac-ncld-1)                      )

 if(.not.allocated(slmsk)    ) allocate(slmsk(nCellsSolve)                       )
 if(.not.allocated(nlons)    ) allocate(nlons(nCellsSolve)                       )
 if(.not.allocated(sinlat)   ) allocate(sinlat(nCellsSolve)                      )
 if(.not.allocated(coslat)   ) allocate(coslat(nCellsSolve)                      )
 if(.not.allocated(rann)     ) allocate(rann(nCellsSolve,nrcm)                   )
 if(.not.allocated(rnp_p)    ) allocate(rnp_p(nCellsSolve,nVertLevels)           )

!allocation of arrays for ozone parameterization:
 if(.not.allocated(poz)      ) allocate(poz(levozp)                              )
 if(.not.allocated(prdout)   ) &
     allocate(prdout(1:nCellsSolve,1:nOznLevels,1:nOznPlCoeff)                   )

!allocation of arrays for parameterization of vertical diffusion (gocart):
 if(.not.allocated(dkh_p)    ) allocate(dkh_p(nCellsSolve,nVertLevels)           )

!allocation of arrays for parameterized convection:
 if(.not.allocated(dpshc)    ) allocate(dpshc(1:nCellsSolve)                     )
 if(.not.allocated(upd_mf_p) ) allocate(upd_mf_p(nCellsSolve,nVertLevels)        )
 if(.not.allocated(dwn_mf_p) ) allocate(dwn_mf_p(nCellsSolve,nVertLevels)        )
 if(.not.allocated(det_mf_p) ) allocate(det_mf_p(nCellsSolve,nVertLevels)        )

!allocation of arrays for 2D and 3D restarts:
 if(.not.allocated(phy_f2d_p)) allocate(phy_f2d_p(nCellsSolve,numP2d)            )
 if(.not.allocated(phy_f3d_p)) allocate(phy_f3d_p(nCellsSolve,nVertLevels,numP3d))

!allocation of arrays for parameterization of long- and short-wave radiation:
 if(.not.allocated(htrlw)    ) allocate(htrlw(nCellsSolve,nVertLevels)           )
 if(.not.allocated(htrsw)    ) allocate(htrsw(nCellsSolve,nVertLevels)           )

!allocation of arrays for parameterization of gravity wave drag over orography:
 if(.not.allocated(hprime)   ) allocate(hprime(nCellsSolve,nGwd)                 )

!allocation of arrays for parameterization of land surface scheme:
 if(.not.allocated(slc)      ) allocate(slc(nCellsSolve,nSoilLevels)             )
 if(.not.allocated(smc)      ) allocate(smc(nCellsSolve,nSoilLevels)             )
 if(.not.allocated(stc)      ) allocate(stc(nCellsSolve,nSoilLevels)             )

 end subroutine allocate_gbphys

!=================================================================================================================
 subroutine deallocate_gbphys
!=================================================================================================================

 if(allocated(dt3dt)    ) deallocate(dt3dt    )
 if(allocated(dq3dt)    ) deallocate(dq3dt    )
 if(allocated(du3dt)    ) deallocate(du3dt    )
 if(allocated(dv3dt)    ) deallocate(dv3dt    )

!deallocation of arrays for miscellaneous:
 if(allocated(fscav)    ) deallocate(fscav    )
 if(allocated(slmsk)    ) deallocate(slmsk    )
 if(allocated(nlons)    ) deallocate(nlons    )
 if(allocated(sinlat)   ) deallocate(sinlat   )
 if(allocated(coslat)   ) deallocate(coslat   )
 if(allocated(rann)     ) deallocate(rann     )
 if(allocated(rnp_p)    ) deallocate(rnp_p    )

!deallocation of arrays for ozone parameterization:
 if(allocated(poz)      ) deallocate(poz      )
 if(allocated(prdout)   ) deallocate(prdout   )

!deallocation of arrays for parameterization of vertical diffusion (gocart):
 if(allocated(dkh_p)    ) deallocate(dkh_p    )

!deallocation of arrays for parameterized convection:
 if(allocated(dpshc)    ) deallocate(dpshc    )
 if(allocated(upd_mf_p) ) deallocate(upd_mf_p )
 if(allocated(dwn_mf_p) ) deallocate(dwn_mf_p )
 if(allocated(det_mf_p) ) deallocate(det_mf_p )

!deallocation of arrays for 2D and 3D restarts:
 if(allocated(phy_f2d_p)) deallocate(phy_f2d_p)
 if(allocated(phy_f3d_p)) deallocate(phy_f3d_p)

!deallocation of arrays for parameterization of long- and short-wave radiation:
 if(allocated(htrlw)    ) deallocate(htrlw    )
 if(allocated(htrsw)    ) deallocate(htrsw    )

!deallocation of arrays for parameterization of gravity wave drag over orography:
 if(allocated(hprime)   ) deallocate(hprime   )

!deallocation of arrays for parameterization of land surface scheme:
 if(allocated(slc)      ) deallocate(slc      )
 if(allocated(smc)      ) deallocate(smc      )
 if(allocated(stc)      ) deallocate(stc      )  

 end subroutine deallocate_gbphys

!=================================================================================================================
 subroutine gbphys_from_MPAS(mesh,atm_input,sfc_input,diag_physics,tend_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: atm_input

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nSoilLevels,numP2d,numP3d,nGwd
 integer,pointer:: nOznLevels,nOznPlCoeff
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer    :: latCell
 real(kind=RKIND),dimension(:),pointer    :: slope,xice
 real(kind=RKIND),dimension(:),pointer    :: zorl,dugwd,dvgwd,acv,acvb,acvt,pl_pres
 real(kind=RKIND),dimension(:),pointer    :: dlwsfci,ulwsfci,dswsfci,uswsfci,dqsfci,dtsfci,epi,gfluxi
 real(kind=RKIND),dimension(:,:),pointer  :: smcrel,smois,tslb
 real(kind=RKIND),dimension(:,:),pointer  :: dkh,rnp,phy_f2d,orogwd
 real(kind=RKIND),dimension(:,:),pointer  :: rthratenlw,rthratensw
 real(kind=RKIND),dimension(:,:,:),pointer:: phy_f3d,o3clim

!local variables:
 integer:: iCell,k,n

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)
 call mpas_pool_get_dimension(mesh,'numP2d'     ,numP2d     )
 call mpas_pool_get_dimension(mesh,'numP3d'     ,numP3d     )
 call mpas_pool_get_dimension(mesh,'nGwd'       ,nGwd       )

 call mpas_pool_get_array(mesh,'latCell',latCell)

 call mpas_pool_get_array(atm_input   ,'pl_pres' ,pl_pres )
 call mpas_pool_get_array(sfc_input   ,'landmask',landmask)
 call mpas_pool_get_array(sfc_input,   'orogwd'  ,orogwd  )
 call mpas_pool_get_array(sfc_input   ,'smcrel'  ,smcrel  )
 call mpas_pool_get_array(sfc_input   ,'smois'   ,smois   )
 call mpas_pool_get_array(sfc_input   ,'tslb'    ,tslb    )
 call mpas_pool_get_array(sfc_input   ,'slope'   ,slope   )
 call mpas_pool_get_array(sfc_input   ,'xice'    ,xice    )
 call mpas_pool_get_array(diag_physics,'acv'     ,acv     )
 call mpas_pool_get_array(diag_physics,'acvb'    ,acvb    )
 call mpas_pool_get_array(diag_physics,'acvt'    ,acvt    )
 call mpas_pool_get_array(diag_physics,'zorl'    ,zorl    )
 call mpas_pool_get_array(diag_physics,'dkh'     ,dkh     )
 call mpas_pool_get_array(diag_physics,'rnp'     ,rnp     )
 call mpas_pool_get_array(diag_physics,'dugwd'   ,dugwd   )
 call mpas_pool_get_array(diag_physics,'dvgwd'   ,dvgwd   )  
 call mpas_pool_get_array(diag_physics,'phy_f2d' ,phy_f2d )
 call mpas_pool_get_array(diag_physics,'phy_f3d' ,phy_f3d )
 call mpas_pool_get_array(diag_physics,'o3clim'  ,o3clim  )
 call mpas_pool_get_array(diag_physics,'dlwsfci' ,dlwsfci )
 call mpas_pool_get_array(diag_physics,'ulwsfci' ,ulwsfci )
 call mpas_pool_get_array(diag_physics,'dswsfci' ,dswsfci )
 call mpas_pool_get_array(diag_physics,'uswsfci' ,uswsfci )
 call mpas_pool_get_array(diag_physics,'dqsfci'  ,dqsfci  )
 call mpas_pool_get_array(diag_physics,'dtsfci'  ,dtsfci  )
 call mpas_pool_get_array(diag_physics,'epi'     ,epi     )
 call mpas_pool_get_array(diag_physics,'gfluxi'  ,gfluxi  )

 call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
 call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

 do iCell = 1, nCellsSolve
    slmsk(iCell) = landmask(iCell)
    if(xice(iCell) .gt. 0._RKIND) then
       slmsk(iCell) = 2._RKIND
       zorl(iCell) = z0_sib(13)
    endif

!temporary:
    slope(iCell) = min(slope(iCell),9._RKIND)

    nlons(iCell)  = nCellsSolve
    sinlat(iCell) = sin(latCell(iCell))
    coslat(iCell) = cos(latCell(iCell))
 enddo

!initialization of local physics tendencies:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       dt3dt(iCell,k,:) = 0._RKIND
       dq3dt(iCell,k,:) = 0._RKIND
       du3dt(iCell,k,:) = 0._RKIND
       dv3dt(iCell,k,:) = 0._RKIND
    enddo
 enddo

 do iCell = 1, nCellsSolve
    !initialization of 2D arrays for restarts:
    do n = 1, numP2d
       phy_f2d_p(iCell,n) = phy_f2d(n,iCell)
    enddo
    do k = 1, nVertLevels
       !initialization of 3D arrays for restarts:
       do n = 1, numP3d
          phy_f3d_p(iCell,k,n) = phy_f3d(n,k,iCell)
       enddo
       !initialization of in cloud precipitation rate:
       rnp_p(iCell,k) = rnp(k,iCell)
       !initialization of arrays for parameterization of vertical diffusion (gocart):
       dkh_p(iCell,k) = dkh(k,iCell)
    enddo
 enddo

!initialization of the maximum pressure depth for shallow convection, others...:
 write(0,*) '--- gen_coord_hybrid =', gen_coord_hybrid
 write(0,*) '--- thermodyn_id     =', thermodyn_id
 do iCell = 1, nCellsSolve
    phil(iCell,nVertLevels) = 0._RKIND !forces the calculation of geopotential in gbphys.
    dpshc(iCell) = 0.3 * pgr_hyd(iCell)

    if(gen_coord_hybrid .and. thermodyn_id == 3) then
       prslk(iCell,1) = 0._RKIND !forces the calculation of geopotential in gbphys.
       prsik(iCell,1) = 0._RKIND !forces the calculation of geopotential in gbphys.
    endif
 enddo

!initialization for parameterization of gravity wave drag over orography:
 do iCell = 1, nCellsSolve
    do n = 1, nGwd
       hprime(iCell,n)  = orogwd(n,iCell)
    enddo
    dugwd(iCell) = 0._RKIND
    dvgwd(iCell) = 0._RKIND
 enddo

!initialization for parameterized convection: CHECK IF INITIALIZATION IS NEEDED FOR ACV,ACVB,ACVT ??
 do iCell = 1, nCellsSolve
!   acv(iCell)  = 0._RKIND
!   acvb(iCell) = 0._RKIND
!   acvt(iCell) = 0._RKIND
    do k = 1, nVertLevels
       upd_mf_p(iCell,k) = 0._RKIND
       dwn_mf_p(iCell,k) = 0._RKIND
       det_mf_p(iCell,k) = 0._RKIND
    enddo
 enddo

!initialization for parameterization of long-wave and short-wave radiation:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       htrlw(iCell,k) = rthratenlw(k,iCell)
       htrsw(iCell,k) = rthratensw(k,iCell)
    enddo
 enddo

!initialization of local liquid soil moisture, soil moisture, and soil temperature for land
!surface scheme:
 do iCell = 1, nCellsSolve
    do k = 1, nSoilLevels
       slc(iCell,k) = smcrel(k,iCell)
       smc(iCell,k) = smois(k,iCell)
       stc(iCell,k) = tslb(k,iCell)
    enddo
 enddo

!initialization for ozone:
 do k = 1, nOznLevels
    poz(k) = pl_pres(k)
 enddo
 do iCell = 1, nCellsSolve
    do k = 1, nOznLevels
       do n = 1, nOznPlCoeff
          prdout(iCell,k,n) = o3clim(k,n,iCell)
       enddo
    enddo
 enddo

!initialization of temperature, momentum, water vapor, water/ice condensates, and ozone mixing
!ratios updated in gbphys:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       gt0(iCell,k)   = 0._RKIND
       gu0(iCell,k)   = 0._RKIND
       gv0(iCell,k)   = 0._RKIND
       gq0(iCell,k,:) = 0._RKIND       
    enddo
 enddo

!outputs to subroutine to gbphys:
 do iCell = 1, nCellsSolve
    dlwsfci(iCell) = 0._RKIND
    ulwsfci(iCell) = 0._RKIND
    dswsfci(iCell) = 0._RKIND
    uswsfci(iCell) = 0._RKIND
    dqsfci(iCell)  = 0._RKIND
    dtsfci(iCell)  = 0._RKIND
    epi(iCell)     = 0._RKIND
    gfluxi(iCell)  = 0._RKIND
 enddo

 end subroutine gbphys_from_MPAS

!=================================================================================================================
 subroutine gbphys_to_MPAS(dminfo,mesh,diag_physics,tend_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nSoilLevels,numP2d,numP3d

 real(kind=RKIND),dimension(:,:),pointer:: thradlw,thradsw,thdiffusion,thcudeep,thcushallow,thmicrop
 real(kind=RKIND),dimension(:,:),pointer:: qvdiffusion,qvcudeep,qvcushallow,qvmicrop
 real(kind=RKIND),dimension(:,:),pointer:: qo3diffusion,qo3product
 real(kind=RKIND),dimension(:,:),pointer:: udiffusion,ugwdo,ucudeep,ugwdc
 real(kind=RKIND),dimension(:,:),pointer:: vdiffusion,vgwdo,vcudeep,vgwdc
 real(kind=RKIND),dimension(:,:),pointer:: thphysics,qvphysics,qo3physics,qcphysics,uphysics,vphysics

 real(kind=RKIND),dimension(:,:),pointer  :: smcrel,smois,tslb
 real(kind=RKIND),dimension(:,:),pointer  :: rnp,dkh,upd_mf,dwn_mf,det_mf
 real(kind=RKIND),dimension(:,:),pointer  :: phy_f2d
 real(kind=RKIND),dimension(:,:,:),pointer:: phy_f3d

!local variables:
 integer:: iCell,k,n
 real(kind=RKIND):: dti

 real(kind=RKIND):: tend_min,tend_max
 real(kind=RKIND):: global_tend_min,global_tend_max

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine gbphys_to_MPAS:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'numP2d'     ,numP2d     )
 call mpas_pool_get_dimension(mesh,'numP3d'     ,numP3d     )

 call mpas_pool_get_array(tend_physics,'thradlw'     ,thradlw     )
 call mpas_pool_get_array(tend_physics,'thradsw'     ,thradsw     )
 call mpas_pool_get_array(tend_physics,'thdiffusion' ,thdiffusion )
 call mpas_pool_get_array(tend_physics,'thcudeep'    ,thcudeep    )
 call mpas_pool_get_array(tend_physics,'thcushallow' ,thcushallow )
 call mpas_pool_get_array(tend_physics,'thmicrop'    ,thmicrop    )
 
 call mpas_pool_get_array(tend_physics,'qvdiffusion' ,qvdiffusion )
 call mpas_pool_get_array(tend_physics,'qvcudeep'    ,qvcudeep    )
 call mpas_pool_get_array(tend_physics,'qvcushallow' ,qvcushallow )
 call mpas_pool_get_array(tend_physics,'qvmicrop'    ,qvmicrop    )
 call mpas_pool_get_array(tend_physics,'qo3diffusion',qo3diffusion)
 call mpas_pool_get_array(tend_physics,'qo3product'  ,qo3product  )

 call mpas_pool_get_array(tend_physics,'udiffusion'  ,udiffusion  )
 call mpas_pool_get_array(tend_physics,'ugwdo'       ,ugwdo       )
 call mpas_pool_get_array(tend_physics,'ucudeep'     ,ucudeep     )
 call mpas_pool_get_array(tend_physics,'ugwdc'       ,ugwdc       )
 call mpas_pool_get_array(tend_physics,'vdiffusion'  ,vdiffusion  )
 call mpas_pool_get_array(tend_physics,'vgwdo'       ,vgwdo       )
 call mpas_pool_get_array(tend_physics,'vcudeep'     ,vcudeep     )
 call mpas_pool_get_array(tend_physics,'vgwdc'       ,vgwdc       )

 call mpas_pool_get_array(tend_physics,'thphysics '  ,thphysics   )
 call mpas_pool_get_array(tend_physics,'qvphysics '  ,qvphysics   )
 call mpas_pool_get_array(tend_physics,'qo3physics'  ,qo3physics  )
 call mpas_pool_get_array(tend_physics,'qcphysics '  ,qcphysics   )
 call mpas_pool_get_array(tend_physics,'uphysics'    ,uphysics    )
 call mpas_pool_get_array(tend_physics,'vphysics'    ,vphysics    )

 call mpas_pool_get_array(sfc_input   ,'smcrel' ,smcrel )
 call mpas_pool_get_array(sfc_input   ,'smois'  ,smois  )
 call mpas_pool_get_array(sfc_input   ,'tslb'   ,tslb   )
 call mpas_pool_get_array(diag_physics,'rnp'    ,rnp    )
 call mpas_pool_get_array(diag_physics,'dkh'    ,dkh    )
 call mpas_pool_get_array(diag_physics,'upd_mf' ,upd_mf )
 call mpas_pool_get_array(diag_physics,'dwn_mf' ,dwn_mf )
 call mpas_pool_get_array(diag_physics,'det_mf' ,det_mf )
 call mpas_pool_get_array(diag_physics,'phy_f2d',phy_f2d)
 call mpas_pool_get_array(diag_physics,'phy_f3d',phy_f3d)

!fill MPAS physics tendencies:
 dti = 1._RKIND/dt_dyn
!potential temperature:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       thradlw(k,iCell)     = dt3dt(iCell,k,1)*dti/prslk(iCell,k)
       thradsw(k,iCell)     = dt3dt(iCell,k,2)*dti/prslk(iCell,k)
       thdiffusion(k,iCell) = dt3dt(iCell,k,3)*dti/prslk(iCell,k)
       thcudeep(k,iCell)    = dt3dt(iCell,k,4)*dti/prslk(iCell,k)
       thcushallow(k,iCell) = dt3dt(iCell,k,5)*dti/prslk(iCell,k)
       thmicrop(k,iCell)    = dt3dt(iCell,k,6)*dti/prslk(iCell,k)
    enddo
 enddo

!water vapor and ozone mixing ratios:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       qvdiffusion(k,iCell)  = dq3dt(iCell,k,1)*dti
       qvcudeep(k,iCell)     = dq3dt(iCell,k,2)*dti
       qvcushallow(k,iCell)  = dq3dt(iCell,k,3)*dti
       qvmicrop(k,iCell)     = dq3dt(iCell,k,4)*dti
       qo3diffusion(k,iCell) = dq3dt(iCell,k,5)*dti
       qo3product(k,iCell)   = dq3dt(iCell,k,6)*dti
    enddo
 enddo
 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_min = min(tend_min,qvdiffusion(k,iCell))
       tend_max = max(tend_max,qvdiffusion(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvdiffusion',global_tend_min,global_tend_max

 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_min = min(tend_min,qvcudeep(k,iCell))
       tend_max = max(tend_max,qvcudeep(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvcudeep',global_tend_min,global_tend_max

 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_min = min(tend_min,qvcushallow(k,iCell))
       tend_max = max(tend_max,qvcushallow(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvcushallow',global_tend_min,global_tend_max

 tend_min =  999._RKIND
 tend_max = -999._RKIND
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       tend_min = min(tend_min,qvmicrop(k,iCell))
       tend_max = max(tend_max,qvmicrop(k,iCell))
    enddo
 enddo
 call mpas_dmpar_min_real(dminfo,tend_min,global_tend_min)
 call mpas_dmpar_min_real(dminfo,tend_max,global_tend_max)
 write(0,*)
 write(0,*) 'global min, max tend_qvmicrop',global_tend_min,global_tend_max

!zonal and meridional winds:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       udiffusion(k,iCell) = du3dt(iCell,k,1)*dti
       ugwdo(k,iCell)      = du3dt(iCell,k,2)*dti 
       ucudeep(k,iCell)    = du3dt(iCell,k,3)*dti
       ugwdc(k,iCell)      = du3dt(iCell,k,4)*dti 

       vdiffusion(k,iCell) = dv3dt(iCell,k,1)*dti
       vgwdo(k,iCell)      = dv3dt(iCell,k,2)*dti
       vcudeep(k,iCell)    = dv3dt(iCell,k,3)*dti 
       vgwdc(k,iCell)      = dv3dt(iCell,k,4)*dti
    enddo
 enddo

!total tendencies using updated temperature,updated water vapor mixing ratio,updated total cloud
!condensate,updated ozone mixing ratio, and updated u- and v-momentum:
 write(0,*) "max tend thradlw     = ",maxval(thradlw(:,1:nCellsSolve)    )
 write(0,*) "min tend thradlw     = ",minval(thradlw(:,1:nCellsSolve)    )
 write(0,*) "max tend thradsw     = ",maxval(thradsw(:,1:nCellsSolve)    )
 write(0,*) "min tend thradsw     = ",minval(thradsw(:,1:nCellsSolve)    )
 write(0,*) "max tend thdiffusion = ",maxval(thdiffusion(:,1:nCellsSolve))
 write(0,*) "min tend thdiffusion = ",minval(thdiffusion(:,1:nCellsSolve))
 write(0,*) "max tend thcudeep    = ",maxval(thcudeep(:,1:nCellsSolve)   )
 write(0,*) "min tend thcudeep    = ",minval(thcudeep(:,1:nCellsSolve)   )
 write(0,*) "max tend thcushallow = ",maxval(thcushallow(:,1:nCellsSolve))
 write(0,*) "min tend thcushallow = ",minval(thcushallow(:,1:nCellsSolve))
 write(0,*) "max tend thmicrop    = ",maxval(thmicrop(:,1:nCellsSolve)   )
 write(0,*) "min tend thmicrop    = ",minval(thmicrop(:,1:nCellsSolve)   )


 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       thphysics(k,iCell)  = (gt0(iCell,k)-tgrs(iCell,k))*dti/prslk(iCell,k)
       qvphysics(k,iCell)  = (gq0(iCell,k,1)-qgrs(iCell,k,1))*dti 
       qo3physics(k,iCell) = (gq0(iCell,k,2)-qgrs(iCell,k,2))*dti
       qcphysics(k,iCell)  = (gq0(iCell,k,3)-qgrs(iCell,k,3))*dti
       uphysics(k,iCell)   = (gu0(iCell,k)-ugrs(iCell,k))*dti
       vphysics(k,iCell)   = (gv0(iCell,k)-vgrs(iCell,k))*dti
    enddo
 enddo
 write(0,*) "max tend thphysics   = ",maxval(thphysics(:,1:nCellsSolve))
 write(0,*) "min tend thphysics   = ",minval(thphysics(:,1:nCellsSolve))
 write(0,*)

 write(0,*)
 write(0,*) '--- numP3d =', numP3d
 do iCell = 1, nCellsSolve
    !update 2D arrays for restarts:
    do n = 1, numP2d
       phy_f2d(n,iCell) = phy_f2d_p(iCell,n)
    enddo
    do k = 1, nVertLevels
       !update 3D arrays for restarts:
       do n = 1, numP3d
          phy_f3d(n,k,iCell) = phy_f3d_p(iCell,k,n)
       enddo
       !update in-cloud precipitation rate:
       rnp(k,iCell) = rnp_p(iCell,k)
       !update vertical diffusion coefficients (gocart):
       dkh(k,iCell) = dkh_p(iCell,k)
    enddo
 enddo

!fill MPAS arrays from parameterized deep convection:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       upd_mf(k,iCell) = upd_mf_p(iCell,k)
       dwn_mf(k,iCell) = dwn_mf_p(iCell,k)
       det_mf(k,iCell) = det_mf_p(iCell,k)
    enddo
 enddo

!update MPAS arrays from land surface scheme:
 do iCell = 1, nCellsSolve
    do k = 1, nSoilLevels
       smcrel(k,iCell) = slc(iCell,k)
       smois(k,iCell)  = smc(iCell,k)
       tslb (k,iCell)  = stc(iCell,k)
    enddo
 enddo

 write(0,*) '--- end subroutine gbphys_to_MPAS'

 end subroutine gbphys_to_MPAS

!=================================================================================================================
 subroutine mpas_atmphys_gbphys(itimestep,xtime_s,dminfo,mesh,atm_input,sfc_input, &
                                diag_physics,tend_physics)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(mpas_pool_type),intent(inout):: atm_input
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nSoilLevels,numP2d,numP3d
 integer,dimension(:),pointer:: isltyp,ivgtyp

 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: slope,shdmin,shdmax,snoalb,skintemp,snow,snowc,snowh, &
                                         tmn,ter,xice
 real(kind=RKIND),dimension(:),pointer:: coszen,sfcdsw,sfcnsw,sfcdlw,tsflw,sfcemis,sfcalb,     &
                                         vegfra,uustar,hice,tisfc,tprcp,cv,cvb,cvt,srflag,     &
                                         zorl,canopy,ffmm,ffhh,f10m,srunoff,evbsa,evcwa,       &
                                         snohfa,transa,sbsnoa,snowca,soilm,tmpmin,tmpmax,      &
                                         dusfc,dvsfc,dtsfc,dqsfc,totprcp,gflux,dlwsfc,ulwsfc,  &
                                         suntim,runoff,ep,cldwrk,dugwd,dvgwd,psmean,cnvprcp,   &
                                         spfhmin,spfhmax,acv,acvb,acvt
 real(kind=RKIND),dimension(:),pointer:: t2m,q2m,u10m,v10m,hpbl,pwat,t1,q1,u1,v1,chh,cmm,      &
                                         dlwsfci,ulwsfci,dswsfci,uswsfci,dtsfci,dqsfci,gfluxi, &
                                         epi,smcwlt2,smcref2,gsoil,gtmp2m,gustar,gpblh,gu10m,  &
                                         gv10m,gzorl,goro

!local variables:
 logical:: flipv

 integer:: iCell,k,n
 integer:: me
 integer:: lat,lonf,latg,kdt,ko3

 real(kind=RKIND),parameter:: cons_24 = 24._RKIND
 real(kind=RKIND),parameter:: cons_99 = 99._RKIND

 real(kind=RKIND):: xtime_h
 real(kind=RKIND):: dtp,dtf,solhr
 real(kind=RKIND):: clstp


!-----------------------------------------------------------------------------------------------------------------
 write(0,100) itimestep

!initialization of local pointers:
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'numP2d'     ,numP2d     )
 call mpas_pool_get_dimension(mesh,'numP3d'     ,numP3d     )

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)

 call mpas_pool_get_array(sfc_input   ,'isltyp'  ,isltyp  )
 call mpas_pool_get_array(sfc_input   ,'ivgtyp'  ,ivgtyp  ) 
 call mpas_pool_get_array(sfc_input   ,'slope'   ,slope   )
 call mpas_pool_get_array(sfc_input   ,'shdmin'  ,shdmin  )
 call mpas_pool_get_array(sfc_input   ,'shdmax'  ,shdmax  )
 call mpas_pool_get_array(sfc_input   ,'snoalb'  ,snoalb  )
 call mpas_pool_get_array(sfc_input   ,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input   ,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input   ,'snowc'   ,snowc   )
 call mpas_pool_get_array(sfc_input   ,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input   ,'ter'     ,ter     )
 call mpas_pool_get_array(sfc_input   ,'tmn'     ,tmn     )
 call mpas_pool_get_array(sfc_input   ,'xice'    ,xice    )

 call mpas_pool_get_array(diag_physics,'vegfra',vegfra    )
 call mpas_pool_get_array(diag_physics,'uustar',uustar    )
 call mpas_pool_get_array(diag_physics,'coszen',coszen    )
 call mpas_pool_get_array(diag_physics,'sfcdsw',sfcdsw    )
 call mpas_pool_get_array(diag_physics,'sfcnsw',sfcnsw    )
 call mpas_pool_get_array(diag_physics,'sfcdlw',sfcdlw    )
 call mpas_pool_get_array(diag_physics,'tsflw' ,tsflw     )
 call mpas_pool_get_array(diag_physics,'sfcemis',sfcemis  )
 call mpas_pool_get_array(diag_physics,'sfcalb' ,sfcalb   )
 call mpas_pool_get_array(diag_physics,'hice'   ,hice     )
 call mpas_pool_get_array(diag_physics,'tisfc'  ,tisfc    )
 call mpas_pool_get_array(diag_physics,'tprcp'  ,tprcp    )
 call mpas_pool_get_array(diag_physics,'cv'     ,cv       )
 call mpas_pool_get_array(diag_physics,'cvb'    ,cvb      )
 call mpas_pool_get_array(diag_physics,'cvt'    ,cvt      )
 call mpas_pool_get_array(diag_physics,'srflag' ,srflag   )
 call mpas_pool_get_array(diag_physics,'zorl'   ,zorl     )
 call mpas_pool_get_array(diag_physics,'canopy' ,canopy   )
 call mpas_pool_get_array(diag_physics,'ffmm'   ,ffmm     )
 call mpas_pool_get_array(diag_physics,'ffhh'   ,ffhh     )
 call mpas_pool_get_array(diag_physics,'f10m'   ,f10m     )
 call mpas_pool_get_array(diag_physics,'srunoff',srunoff  )
 call mpas_pool_get_array(diag_physics,'evbsa'  ,evbsa    )
 call mpas_pool_get_array(diag_physics,'evcwa'  ,evcwa    )
 call mpas_pool_get_array(diag_physics,'snohfa' ,snohfa   )
 call mpas_pool_get_array(diag_physics,'transa' ,transa   )
 call mpas_pool_get_array(diag_physics,'sbsnoa' ,sbsnoa   )
 call mpas_pool_get_array(diag_physics,'snowca' ,snowca   )
 call mpas_pool_get_array(diag_physics,'soilm'  ,soilm    )
 call mpas_pool_get_array(diag_physics,'tmpmin' ,tmpmin   )
 call mpas_pool_get_array(diag_physics,'tmpmax' ,tmpmax   )
 call mpas_pool_get_array(diag_physics,'dusfc'  ,dusfc    )
 call mpas_pool_get_array(diag_physics,'dvsfc'  ,dvsfc    )
 call mpas_pool_get_array(diag_physics,'dtsfc'  ,dtsfc    )
 call mpas_pool_get_array(diag_physics,'dqsfc'  ,dqsfc    )
 call mpas_pool_get_array(diag_physics,'totprcp',totprcp  )
 call mpas_pool_get_array(diag_physics,'gflux'  ,gflux    )
 call mpas_pool_get_array(diag_physics,'dlwsfc' ,dlwsfc   )
 call mpas_pool_get_array(diag_physics,'ulwsfc' ,ulwsfc   )
 call mpas_pool_get_array(diag_physics,'suntim' ,suntim   )
 call mpas_pool_get_array(diag_physics,'runoff' ,runoff   )
 call mpas_pool_get_array(diag_physics,'ep'     ,ep       )
 call mpas_pool_get_array(diag_physics,'cldwrk' ,cldwrk   )
 call mpas_pool_get_array(diag_physics,'dugwd'  ,dugwd    )
 call mpas_pool_get_array(diag_physics,'dvgwd'  ,dvgwd    )
 call mpas_pool_get_array(diag_physics,'psmean' ,psmean   )
 call mpas_pool_get_array(diag_physics,'cnvprcp',cnvprcp  )
 call mpas_pool_get_array(diag_physics,'spfhmin',spfhmin  )
 call mpas_pool_get_array(diag_physics,'spfhmax',spfhmax  )
 call mpas_pool_get_array(diag_physics,'acv'    ,acv      )
 call mpas_pool_get_array(diag_physics,'acvb'   ,acvb     )
 call mpas_pool_get_array(diag_physics,'acvt'   ,acvt     )

 call mpas_pool_get_array(diag_physics,'t2m'    ,t2m      )
 call mpas_pool_get_array(diag_physics,'q2m'    ,q2m      )
 call mpas_pool_get_array(diag_physics,'u10m'   ,u10m     )
 call mpas_pool_get_array(diag_physics,'v10m'   ,v10m     )
 call mpas_pool_get_array(diag_physics,'hpbl'   ,hpbl     )
 call mpas_pool_get_array(diag_physics,'pwat'   ,pwat     )
 call mpas_pool_get_array(diag_physics,'t1'     ,t1       )
 call mpas_pool_get_array(diag_physics,'q1'     ,q1       )
 call mpas_pool_get_array(diag_physics,'u1'     ,u1       )
 call mpas_pool_get_array(diag_physics,'v1'     ,v1       )
 call mpas_pool_get_array(diag_physics,'chh'    ,chh      )
 call mpas_pool_get_array(diag_physics,'cmm'    ,cmm      )
 call mpas_pool_get_array(diag_physics,'dlwsfci',dlwsfci  )
 call mpas_pool_get_array(diag_physics,'ulwsfci',ulwsfci  )
 call mpas_pool_get_array(diag_physics,'dswsfci',dswsfci  )
 call mpas_pool_get_array(diag_physics,'uswsfci',uswsfci  )
 call mpas_pool_get_array(diag_physics,'dtsfci' ,dtsfci   )
 call mpas_pool_get_array(diag_physics,'dqsfci' ,dqsfci   )
 call mpas_pool_get_array(diag_physics,'gfluxi' ,gfluxi   )
 call mpas_pool_get_array(diag_physics,'epi'    ,epi      )
 call mpas_pool_get_array(diag_physics,'smcwlt2',smcwlt2  )
 call mpas_pool_get_array(diag_physics,'smcref2',smcref2  )
 call mpas_pool_get_array(diag_physics,'gsoil'  ,gsoil    )
 call mpas_pool_get_array(diag_physics,'gtmp2m' ,gtmp2m   )
 call mpas_pool_get_array(diag_physics,'gustar' ,gustar   )
 call mpas_pool_get_array(diag_physics,'gpblh'  ,gpblh    )
 call mpas_pool_get_array(diag_physics,'gu10m'  ,gu10m    )
 call mpas_pool_get_array(diag_physics,'gv10m'  ,gv10m    )
 call mpas_pool_get_array(diag_physics,'gzorl'  ,gzorl    )
 call mpas_pool_get_array(diag_physics,'goro'   ,goro     )

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

!initialization of local variables:
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
    lonCell(1:nCellsSolve)      , latCell(1:nCellsSolve)      , slope(1:nCellsSolve)   ,        &
    shdmin(1:nCellsSolve)       , shdmax(1:nCellsSolve)       , snoalb(1:nCellsSolve)  ,        &
    tmn(1:nCellsSolve)          , slmsk                       , vegfra(1:nCellsSolve)  ,        & 
    DBLE(ivgtyp(1:nCellsSolve)) , DBLE(isltyp(1:nCellsSolve)) , uustar(1:nCellsSolve)  ,        &
    ter(1:nCellsSolve)          , coszen(1:nCellsSolve)       , sfcdsw(1:nCellsSolve)  ,        &
    sfcnsw(1:nCellsSolve)       , sfcdlw(1:nCellsSolve)       , tsflw(1:nCellsSolve)   ,        &
    sfcemis(1:nCellsSolve)      , sfcalb(1:nCellsSolve)       ,                                 &
    htrsw, htrlw, ras, pre_rad, ldiag3d, lggfs3d, lssav, lssav_cc, bkgd_vdif_m, bkgd_vdif_h,    &
    bkgd_vdif_s, psautco, prautco, evpco, wminco, flipv, old_monin, cnvgwd, shal_cnv, sashal,   &
    newsas, cal_pre, mom4ice, mstrat, trans_trac, nst_fcst, moist_adj, fscav, thermodyn_id,     & 
    sfcpress_id, gen_coord_hybrid,                                                              &
    !... inputs/outputs:
    hice(1:nCellsSolve)         , xice(1:nCellsSolve)         , tisfc(1:nCellsSolve)   ,        &
    skintemp(1:nCellsSolve)     , tprcp(1:nCellsSolve)        , cv(1:nCellsSolve)      ,        &
    cvb(1:nCellsSolve)          , cvt(1:nCellsSolve)          , srflag(1:nCellsSolve)  ,        &
    snowh(1:nCellsSolve)        , snow(1:nCellsSolve)         , snowc(1:nCellsSolve)   ,        &
    zorl(1:nCellsSolve)         , canopy(1:nCellsSolve)       , ffmm(1:nCellsSolve)    ,        &
    ffhh(1:nCellsSolve)         , f10m(1:nCellsSolve)         , srunoff(1:nCellsSolve) ,        &
    evbsa(1:nCellsSolve)        , evcwa(1:nCellsSolve)        , snohfa(1:nCellsSolve)  ,        &
    transa(1:nCellsSolve)       , sbsnoa(1:nCellsSolve)       , snowca(1:nCellsSolve)  ,        &
    soilm(1:nCellsSolve)        , tmpmin(1:nCellsSolve)       , tmpmax(1:nCellsSolve)  ,        &
    dusfc(1:nCellsSolve)        , dvsfc(1:nCellsSolve)        , dtsfc(1:nCellsSolve)   ,        &
    dqsfc(1:nCellsSolve)        , totprcp(1:nCellsSolve)      , gflux(1:nCellsSolve)   ,        &
    dlwsfc(1:nCellsSolve)       , ulwsfc(1:nCellsSolve)       , suntim(1:nCellsSolve)  ,        &
    runoff(1:nCellsSolve)       , ep(1:nCellsSolve)           , cldwrk(1:nCellsSolve)  ,        &
    dugwd(1:nCellsSolve)        , dvgwd(1:nCellsSolve)        , psmean(1:nCellsSolve)  ,        &    
    cnvprcp(1:nCellsSolve)      , spfhmin(1:nCellsSolve)      , spfhmax(1:nCellsSolve) ,        &
    dt3dt , dq3dt , du3dt , dv3dt ,                                                             &
    acv(1:nCellsSolve)          , acvb(1:nCellsSolve)         , acvt(1:nCellsSolve)    ,        &
    slc, smc, stc, upd_mf_p, dwn_mf_p, det_mf_p, dkh_p, rnp_p, phy_f3d_p, phy_f2d_p,            &
    dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys,                     &
    dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys, dumphys,                     &
    dumphys, dumphys, dumphys, dumphys,                                                         &
    !... outputs:
    gt0, gq0, gu0, gv0,                                                                         &
    t2m(1:nCellsSolve)          , q2m(1:nCellsSolve)          , u10m(1:nCellsSolve)    ,        &
    v10m(1:nCellsSolve)         , dumphys                     , dumphys                ,        &
    hpbl(1:nCellsSolve)         , pwat(1:nCellsSolve)         , t1(1:nCellsSolve)      ,        &
    q1(1:nCellsSolve)           , u1(1:nCellsSolve)           , v1(1:nCellsSolve)      ,        &
    chh(1:nCellsSolve)          , cmm(1:nCellsSolve)          , dlwsfci(1:nCellsSolve) ,        &
    ulwsfci(1:nCellsSolve)      , dswsfci(1:nCellsSolve)      , uswsfci(1:nCellsSolve) ,        &
    dtsfci(1:nCellsSolve)       , dqsfci(1:nCellsSolve)       , gfluxi(1:nCellsSolve)  ,        &
    epi(1:nCellsSolve)          , smcwlt2(1:nCellsSolve)      , smcref2(1:nCellsSolve) ,        &
    gsoil(1:nCellsSolve)        , gtmp2m(1:nCellsSolve)       , gustar(1:nCellsSolve)  ,        &
    gpblh(1:nCellsSolve)        , gu10m(1:nCellsSolve)        , gv10m(1:nCellsSolve)   ,        &
    gzorl(1:nCellsSolve)        , goro(1:nCellsSolve)         , dumphys                ,        &
    dumphys                     , dumphys                     , dumphys                ,        &
    dumphys                     , dumphys                     , dumphys                ,        &
    dumphys                     , dumphys                     , dumphys                ,        &
    dumphys                     , dumphys                                                       &
            )
 write(0,*) '--- end subroutine gbphys:'

!... copy all GFS arrays to MPAS arrays:
 call gbphys_to_MPAS(dminfo,mesh,diag_physics,tend_physics,sfc_input)

 write(0,*) '--- end subroutine mpas_atmphys_gbphys:'

 end subroutine mpas_atmphys_gbphys

!=================================================================================================================
 end module mpas_atmphys_gfs_driver_gbphys
!=================================================================================================================
