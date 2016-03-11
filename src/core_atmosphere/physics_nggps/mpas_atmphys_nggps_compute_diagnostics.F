! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_compute_diagnostics
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_timekeeping

 use funcphys
 use physcons, only: eps => con_eps,epsm1 => con_epsm1
 use mpas_atmphys_nggps_constants
 use mpas_atmphys_nggps_vars, only: l_diags

 implicit none
 private
 public:: dynamics_nggps_compute_diagnostics, &
          physics_nggps_compute_diagnostics,  &
          physics_nggps_fill_diagnostics,     &
          physics_nggps_zero_diagnostics,     &
          physics_nggps_compute_relhum

!modifications to sourcecode:
!----------------------------
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.
! * added dynamics and physics diagnostics for nggps diagnostics. Outputs of those diagnostics to I/O are based
!   on the output_interval of the stream "nggps_tavg" in streams.atmosphere.
!   Laura D. Fowler (laura@ucar.edu) / 2016-02-07.

 contains
 

!=================================================================================================================
 subroutine physics_nggps_fill_diagnostics(mesh,nggps_sfc_properties,nggps_tbd_ddt,nggps_diagnostics_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_sfc_properties
 type(mpas_pool_type),intent(in):: nggps_tbd_ddt

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_diagnostics_physics

!local pointers:
 integer,pointer:: nCellsSolve

 integer,dimension(:),pointer:: vtypesfc,stypesfc
 real(kind=RKIND),dimension(:),pointer:: icesfc,slmsksfc,spfh2m,tmp2m,tmpsfc,weasdsfc,snodsfc,zorlsfc,vfracsfc, &
                                         f10msfc
 real(kind=RKIND),dimension(:),pointer:: soilt1,soilt2,soilt3,soilt4,soilw1,soilw2,soilw3,soilw4

 real(kind=PHYSKIND),dimension(:),pointer  :: fice,slmsk,q2m,t2m,tsfc,weasd,snwdph,zorl,vfrac,f10m,vtype,stype
 real(kind=PHYSKIND),dimension(:,:),pointer:: slc,stc

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_nggps_fill_diagnostics:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(nggps_sfc_properties,'fice'  ,fice  )
 call mpas_pool_get_array(nggps_sfc_properties,'slmsk' ,slmsk )
 call mpas_pool_get_array(nggps_sfc_properties,'q2m'   ,q2m   )
 call mpas_pool_get_array(nggps_sfc_properties,'t2m'   ,t2m   )
 call mpas_pool_get_array(nggps_sfc_properties,'tsfc'  ,tsfc  )
 call mpas_pool_get_array(nggps_sfc_properties,'weasd' ,weasd )
 call mpas_pool_get_array(nggps_sfc_properties,'snwdph',snwdph)
 call mpas_pool_get_array(nggps_sfc_properties,'zorl'  ,zorl  )
 call mpas_pool_get_array(nggps_sfc_properties,'vfrac' ,vfrac )
 call mpas_pool_get_array(nggps_sfc_properties,'f10m'  ,f10m  )
 call mpas_pool_get_array(nggps_sfc_properties,'vtype' ,vtype )
 call mpas_pool_get_array(nggps_sfc_properties,'stype' ,stype )

 call mpas_pool_get_array(nggps_tbd_ddt,'slc',slc)
 call mpas_pool_get_array(nggps_tbd_ddt,'stc',stc)

 call mpas_pool_get_array(nggps_diagnostics_physics,'icesfc'  ,icesfc  )
 call mpas_pool_get_array(nggps_diagnostics_physics,'slmsksfc',slmsksfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'spfh2m'  ,spfh2m  )
 call mpas_pool_get_array(nggps_diagnostics_physics,'tmp2m'   ,tmp2m   )
 call mpas_pool_get_array(nggps_diagnostics_physics,'tmpsfc'  ,tmpsfc  )
 call mpas_pool_get_array(nggps_diagnostics_physics,'weasdsfc',weasdsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'snodsfc' ,snodsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'zorlsfc' ,zorlsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'vfracsfc',vfracsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'f10msfc' ,f10msfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'vtypesfc',vtypesfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'stypesfc',stypesfc)

 call mpas_pool_get_array(nggps_diagnostics_physics,'soilt1',soilt1)
 call mpas_pool_get_array(nggps_diagnostics_physics,'soilt2',soilt2)
 call mpas_pool_get_array(nggps_diagnostics_physics,'soilt3',soilt3)
 call mpas_pool_get_array(nggps_diagnostics_physics,'soilt4',soilt4)
 call mpas_pool_get_array(nggps_diagnostics_physics,'soilw1',soilw1)
 call mpas_pool_get_array(nggps_diagnostics_physics,'soilw2',soilw2)
 call mpas_pool_get_array(nggps_diagnostics_physics,'soilw3',soilw3)
 call mpas_pool_get_array(nggps_diagnostics_physics,'soilw4',soilw4)

 do iCell = 1, nCellsSolve
    icesfc(iCell) = 0._RKIND
 enddo

 do iCell = 1, nCellsSolve
    if(fice(iCell) .gt. 0._RKIND) icesfc(iCell) = 1._RKIND
    slmsksfc(iCell) = slmsk(iCell)
    spfh2m(iCell)   = q2m(iCell)
    tmp2m(iCell)    = t2m(iCell)
    tmpsfc(iCell)   = tsfc(iCell)
    weasdsfc(iCell) = weasd(iCell)
    snodsfc(iCell)  = snwdph(iCell) * 0.001_RKIND
    zorlsfc(iCell)  = zorl(iCell) * 0.01_RKIND
    vfracsfc(iCell) = vfrac(iCell)
    f10msfc(iCell)  = f10m(iCell)
    vtypesfc(iCell) = int(vtype(iCell))
    stypesfc(iCell) = int(stype(iCell))
 enddo

 do iCell = 1, nCellsSolve
    soilt1(iCell) = stc(1,iCell) 
    soilt2(iCell) = stc(2,iCell)
    soilt3(iCell) = stc(3,iCell)
    soilt4(iCell) = stc(4,iCell)
    soilw1(iCell) = slc(1,iCell) 
    soilw2(iCell) = slc(2,iCell)
    soilw3(iCell) = slc(3,iCell)
    soilw4(iCell) = slc(4,iCell)
 enddo

 write(0,*) '--- end subroutine physics_nggps_fill_diagnostics:'

 end subroutine physics_nggps_fill_diagnostics

!=================================================================================================================
 subroutine dynamics_nggps_compute_diagnostics(mesh,state,time_lev,diag,nggps_diagnostics_dynamics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_diagnostics_dynamics

!local pointers:
 integer,pointer:: nCells,nCellsSolve,nVertLevels
 integer,pointer:: index_qv,index_qc,index_qo3

 real(kind=RKIND),dimension(:,:),pointer  :: zgrid,zz
 real(kind=RKIND),dimension(:,:),pointer  :: rho_zz
 real(kind=RKIND),dimension(:,:),pointer  :: exner,pressure_b,pressure_p
 real(kind=RKIND),dimension(:,:),pointer  :: uZonal,uMeridional
 real(kind=PHYSKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND),dimension(:),pointer  :: hgtsfc,pressfc
 real(kind=RKIND),dimension(:,:),pointer:: theta_m
 real(kind=PHYSKIND),dimension(:,:),pointer:: qv,qc,qo3
 real(kind=RKIND),dimension(:,:),pointer:: ugrd,vgrd,dpres,pres
 real(kind=RKIND),dimension(:,:),pointer:: temperature,spfh,o3mr,clwmr

!local variables:
 integer:: iCell,k
 real(kind=RKIND):: rho1,rho2,tem1,tem2
 real(kind=RKIND):: w1,w2,z0,z1,z2

 real(kind=RKIND),dimension(:,:),allocatable:: fzm,fzp,prsi

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine dynamics_nggps_compute_dynamics:'

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'zgrid',zgrid)
 call mpas_pool_get_array(mesh,'zz'   ,zz   )

 call mpas_pool_get_dimension(state,'index_qv' ,index_qv )
 call mpas_pool_get_dimension(state,'index_qc' ,index_qc )
 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)

 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv  => scalars(index_qv,:,:)
 qc  => scalars(index_qc,:,:)
 qo3 => scalars(index_qo3,:,:)

 call mpas_pool_get_array(diag,'exner'                 ,exner      )
 call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p )
 call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b )
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,uZonal     )
 call mpas_pool_get_array(diag,'uReconstructMeridional',uMeridional)

 call mpas_pool_get_array(nggps_diagnostics_dynamics,'temperature',temperature)
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'pres'   ,pres   )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'dpres'  ,dpres  )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'ugrd'   ,ugrd   )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'vgrd'   ,vgrd   )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'spfh'   ,spfh   )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'o3mr'   ,o3mr   )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'clwmr'  ,clwmr  )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'hgtsfc' ,hgtsfc )
 call mpas_pool_get_array(nggps_diagnostics_dynamics,'pressfc',pressfc)

!...
 do k = 1,nVertLevels
 do iCell = 1, nCellsSolve
    ugrd(k,iCell) = uZonal(k,iCell)
    vgrd(k,iCell) = uMeridional(k,iCell)

    pres(k,iCell) = pressure_p(k,iCell) + pressure_b(k,iCell)

    temperature(k,iCell) = theta_m(k,iCell)/(1._RKIND+R_v/R_d*qv(k,iCell))
    temperature(k,iCell) = temperature(k,iCell)*exner(k,iCell)

    spfh(k,iCell)  = qv(k,iCell)  / (1._RKIND+qv(k,iCell))
    o3mr(k,iCell)  = qo3(k,iCell) / (1._RKIND+qv(k,iCell))      
    clwmr(k,iCell) = qc(k,iCell)  / (1._RKIND+qv(k,iCell))
 enddo
 enddo

 301 format(i9,i4,10(1x,e15.8))
!do iCell = 1, 1
!do k = 1, nVertLevels
!   write(0,301) iCell,k,ugrd(k,iCell),vgrd(k,iCell),pres(k,iCell)/100.,temperature(k,iCell), &
!                spfh(k,iCell),o3mr(k,iCell),clwmr(k,iCell)
!enddo
!write(0,*)
!enddo

!... pressure at layer interfaces:
 if(.not.allocated(fzm) ) allocate(fzm(1:nVertLevels,nCells)  )
 if(.not.allocated(fzp) ) allocate(fzp(1:nVertLevels,nCells)  )
 if(.not.allocated(prsi)) allocate(prsi(1:nVertLevels+1,nCells))

 do k = 2, nVertLevels
 do iCell = 1, nCellsSolve
    tem1 = 1./(zgrid(k+1,iCell)-zgrid(k-1,iCell))
    fzm(k,iCell)  = (zgrid(k,iCell)-zgrid(k-1,iCell)) * tem1
    fzp(k,iCell)  = (zgrid(k+1,iCell)-zgrid(k,iCell)) * tem1
    prsi(k,iCell) = fzm(k,iCell)*pres(k,iCell) + fzp(k,iCell)*pres(k-1,iCell)
 enddo
 enddo
!... extrapolation of temperature and pressure to the top-of-the-model:
 k = nVertLevels+1
 do iCell = 1, nCellsSolve
    z0 = zgrid(k,iCell)
    z1 = 0.5_RKIND*(zgrid(k,iCell)+zgrid(k-1,iCell))
    z2 = 0.5_RKIND*(zgrid(k-1,iCell)+zgrid(k-2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    prsi(k,iCell) = exp(w1*log(pres(k-1,iCell)) + w2*log(pres(k-2,iCell)))
 enddo
!... extrapolation of temperature and pressure to the surface; calculation of
!the surface pressure:
 k = 1
 do iCell = 1, nCellsSolve
    z0 = zgrid(k,iCell)
    z1 = 0.5*(zgrid(k,iCell)+zgrid(k+1,iCell))
    z2 = 0.5*(zgrid(k+1,iCell)+zgrid(k+2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    tem1 = zgrid(2,iCell)-zgrid(1,iCell)
    tem2 = zgrid(3,iCell)-zgrid(2,iCell)
    rho1 = rho_zz(1,iCell)*zz(1,iCell)*(1._RKIND + qv(1,iCell))
    rho2 = rho_zz(2,iCell)*zz(2,iCell)*(1._RKIND + qv(2,iCell))
    prsi(k,iCell) = 0.5*gravity*(zgrid(2,iCell)-zgrid(1,iCell)) &
                  * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    prsi(k,iCell) = prsi(k,iCell) + pres(k,iCell)
 enddo

!... layer thickness:
 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    dpres(k,iCell) = prsi(k,iCell)-prsi(k+1,iCell)
 enddo
 enddo

 if(allocated(fzm) ) deallocate(fzm)
 if(allocated(fzp) ) deallocate(fzp)
 if(allocated(prsi)) deallocate(prsi)

!... surface geopotential:
 do iCell = 1, nCellsSolve
    hgtsfc(iCell) = gravity*zgrid(1,iCell)
 enddo

!... surface pressure using hydrostatic assumption down to the surface:
 do iCell = 1, nCellsSolve
    tem1 = zgrid(2,iCell) - zgrid(1,iCell)
    tem2 = zgrid(3,iCell) - zgrid(2,iCell)

    rho1 = rho_zz(1,iCell)*zz(1,iCell)*(1._RKIND + qv(1,iCell))
    rho2 = rho_zz(2,iCell)*zz(2,iCell)*(1._RKIND + qv(2,iCell))

    pressfc(iCell) = 0.5*gravity*(zgrid(2,iCell)-zgrid(1,iCell)) &
                   * (rho1 + 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    pressfc(iCell) = pressfc(iCell) + pressure_p(1,iCell) + pressure_b(1,iCell)
 enddo

 end subroutine dynamics_nggps_compute_diagnostics

!=================================================================================================================
 subroutine physics_nggps_compute_diagnostics(configs,mesh,nggps_interval,nggps_diagnostics, &
                                  nggps_diagnostics_radiation,nggps_diagnostics_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: nggps_diagnostics

 character(len=StrKIND),intent(in):: nggps_interval

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_diagnostics_radiation
 type(mpas_pool_type),intent(inout):: nggps_diagnostics_physics

!local pointers:
 character(len=StrKIND),pointer:: config_radtlw_interval, &
                                  config_radtsw_interval

 integer,pointer:: nFxr,nCellsSolve
 real(kind=RKIND),dimension(:),pointer:: albdosfc,dlwrfsfc,ulwrfsfc,dswrfsfc,uswrfsfc,ulwrftoa,dswrftoa,uswrftoa
 real(kind=RKIND),dimension(:),pointer:: tcdcclm,tcdchcl,tcdcmcl,tcdclcl 

 real(kind=RKIND),dimension(:),pointer:: cnvprcp,totprcp,gflux,dqsfc,dtsfc,dugwd,dvgwd,dusfc,dvsfc
 real(kind=RKIND),dimension(:),pointer:: cpratsfc,pratesfc,gfluxsfc,lhtflsfc,shtflsfc,ugwdsfc,vgwdsfc, &
                                         uflxsfc,vflxsfc

 real(kind=RKIND),dimension(:,:),pointer:: fluxr

!local variables:
 type(MPAS_TimeInterval_Type):: intv

 integer:: iCell,iFx
 real(kind=RKIND):: dt_radlw,dt_radsw,dt_diags
 real(kind=RKIND):: ratio_lw,ratio_sw,ratio_cl

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_nggps_compute_diagnostics:'

 call mpas_pool_get_config(configs,'config_radtlw_interval',config_radtlw_interval)
 call mpas_pool_get_config(configs,'config_radtsw_interval',config_radtsw_interval)

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nFxr'       ,nFxr       )

!... radiation diagnostics:
 call mpas_pool_get_array(nggps_diagnostics,'fluxr',fluxr)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'albdosfc',albdosfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'dlwrfsfc',dlwrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'ulwrfsfc',ulwrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'dswrfsfc',dswrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'uswrfsfc',uswrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'ulwrftoa',ulwrftoa)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'dswrftoa',dswrftoa)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'uswrftoa',uswrftoa)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdcclm' ,tcdcclm )
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdchcl' ,tcdchcl )
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdcmcl' ,tcdcmcl )
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdclcl' ,tcdclcl )

!... physics diagnostics:
 call mpas_pool_get_array(nggps_diagnostics,'cnvprcp',cnvprcp)
 call mpas_pool_get_array(nggps_diagnostics,'totprcp',totprcp)
 call mpas_pool_get_array(nggps_diagnostics,'gflux'  ,gflux  )
 call mpas_pool_get_array(nggps_diagnostics,'dqsfc'  ,dqsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dtsfc'  ,dtsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dugwd'  ,dugwd  )
 call mpas_pool_get_array(nggps_diagnostics,'dvgwd'  ,dvgwd  )
 call mpas_pool_get_array(nggps_diagnostics,'dusfc'  ,dusfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dvsfc'  ,dvsfc  )

 call mpas_pool_get_array(nggps_diagnostics_physics,'cpratsfc',cpratsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'pratesfc',pratesfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'gfluxsfc',gfluxsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'lhtflsfc',lhtflsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'shtflsfc',shtflsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'ugwdsfc' ,ugwdsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'vgwdsfc' ,vgwdsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'uflxsfc' ,uflxsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'vflxsfc' ,vflxsfc )

 call mpas_set_timeInterval(intv,timeString=config_radtlw_interval)
 call mpas_get_timeInterval(intv,dt=dt_radlw)
 call mpas_set_timeInterval(intv,timeString=config_radtsw_interval)
 call mpas_get_timeInterval(intv,dt=dt_radsw)
 call mpas_set_timeInterval(intv,timeString=nggps_interval)
 call mpas_get_timeInterval(intv,dt=dt_diags)

 ratio_lw = dt_radlw / dt_diags
 ratio_sw = dt_radsw / dt_diags
 ratio_cl = min(dt_radlw,dt_radsw) / dt_diags

 write(0,*) '--- config_radtlw_interval = ', dt_radlw
 write(0,*) '--- config_radtsw_interval = ', dt_radsw
 write(0,*) '--- nggps_interval         = ', dt_diags
!write(0,*)
!write(0,*) '--- ratio_lw = ', ratio_lw
!write(0,*) '--- ratio_sw = ', ratio_sw
!write(0,*) '--- ratio_cl = ', ratio_cl
!write(0,*)
!do iCell = 1,20
!   write(0,301) iCell,fluxr(19,iCell),fluxr(20,iCell),fluxr(4,iCell),fluxr(3,iCell), &
!      fluxr(23,iCell),fluxr(2,iCell),fluxr(1,iCell)
!enddo
!write(0,*)
!do iCell = 1,20
!   write(0,301) iCell,fluxr(17,iCell),fluxr(5,iCell),fluxr(7,iCell),fluxr(6,iCell)
!enddo
 301 format(i6,11(1x,e15.8))

!--- average the long-wave and short-wave fluxes, and the total, high, middle, and low cloud cover over the
!    output time interval:
 do iCell = 1, nCellsSolve
    dlwrfsfc(iCell) = fluxr(19,iCell) / dt_diags
    ulwrfsfc(iCell) = fluxr(20,iCell) / dt_diags
    ulwrftoa(iCell) = fluxr(1,iCell)  / dt_diags

    dswrfsfc(iCell) = fluxr(4,iCell)  / dt_diags
    uswrfsfc(iCell) = fluxr(3,iCell)  / dt_diags
    dswrftoa(iCell) = fluxr(23,iCell) / dt_diags
    uswrftoa(iCell) = fluxr(2,iCell)  / dt_diags

    if(dswrfsfc(iCell) .gt. 0._RKIND) albdosfc(iCell) = fluxr(3,iCell) / fluxr(4,iCell)
    albdosfc(iCell) = albdosfc(iCell) * 100._RKIND

    tcdcclm(iCell)  = fluxr(17,iCell) * 100._RKIND / dt_diags
    tcdchcl(iCell)  = fluxr(5,iCell)  * 100._RKIND / dt_diags
    tcdcmcl(iCell)  = fluxr(6,iCell)  * 100._RKIND / dt_diags
    tcdclcl(iCell)  = fluxr(7,iCell)  * 100._RKIND / dt_diags
 enddo
!write(0,*)
!do iCell = 1, 20
!   write(0,301) iCell,dlwrfsfc(iCell),ulwrfsfc(iCell),dswrfsfc(iCell),uswrfsfc(iCell), &
!      dswrftoa(iCell),uswrftoa(iCell),ulwrftoa(iCell)
!enddo
!write(0,*)
!do iCell = 1, 20
!   write(0,301) iCell,tcdcclm(iCell),tcdchcl(iCell),tcdclcl(iCell),tcdcmcl(iCell)
!enddo

!... average the physics fluxes over the output time interval:
 do iCell = 1, nCellsSolve
    cpratsfc(iCell) = cnvprcp(iCell) / dt_diags
    pratesfc(iCell) = totprcp(iCell) / dt_diags
    gfluxsfc(iCell) = gflux(iCell)   / dt_diags
    lhtflsfc(iCell) = dqsfc(iCell)   / dt_diags
    shtflsfc(iCell) = dtsfc(iCell)   / dt_diags
    ugwdsfc(iCell)  = dugwd(iCell)   / dt_diags
    vgwdsfc(iCell)  = dvgwd(iCell)   / dt_diags
    uflxsfc(iCell)  = dusfc(iCell)   / dt_diags
    vflxsfc(iCell)  = dvsfc(iCell)   / dt_diags
 enddo 
 
!... other surface diagnostics:
 



 write(0,*) '--- exit subroutine physics_nggps_compute_diagnostics'
 write(0,*)

 end subroutine physics_nggps_compute_diagnostics

!=================================================================================================================
 subroutine physics_nggps_zero_diagnostics(mesh,nggps_diagnostics,nggps_diagnostics_physics, &
                                           nggps_diagnostics_radiation)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_diagnostics
 type(mpas_pool_type),intent(inout):: nggps_diagnostics_physics
 type(mpas_pool_type),intent(inout):: nggps_diagnostics_radiation

!local pointers:
 integer,pointer:: nFxr,nCellsSolve,nVertLevels

 real(kind=RKIND),dimension(:),pointer:: psmean,suntim,dlwsfc,ulwsfc,gflux,evbsa,evcwa,transa,  &
                                         sbsnoa,snowca,snohfa,ep,dusfc,dvsfc,dtsfc,dqsfc,dugwd, &  
                                         dvgwd,cldwrk,cnvprcp,totprcp,runoff,srunoff,tmpmax,    &
                                         tmpmin,spfhmax,spfhmin

 real(kind=RKIND),dimension(:),pointer:: albdosfc,dlwrfsfc,ulwrfsfc,dswrfsfc,uswrfsfc,ulwrftoa,dswrftoa,uswrftoa
 real(kind=RKIND),dimension(:),pointer:: tcdcclm,tcdchcl,tcdcmcl,tcdclcl 

 real(kind=RKIND),dimension(:),pointer:: cpratsfc,pratesfc,gfluxsfc,lhtflsfc,shtflsfc,ugwdsfc,vgwdsfc, &
                                         uflxsfc,vflxsfc

 real(kind=RKIND),dimension(:,:),pointer:: fluxr

!local variables:
 integer:: iCell,iFx,k

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine physics_nggps_update_diagnostics:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nFxr'       ,nFxr       )

!... radiation diagnostics:
 call mpas_pool_get_array(nggps_diagnostics,'fluxr',fluxr)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'albdosfc',albdosfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'dlwrfsfc',dlwrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'ulwrfsfc',ulwrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'dswrfsfc',dswrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'uswrfsfc',uswrfsfc)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'ulwrftoa',ulwrftoa)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'dswrftoa',dswrftoa)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'uswrftoa',uswrftoa)
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdcclm' ,tcdcclm )
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdchcl' ,tcdchcl )
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdcmcl' ,tcdcmcl )
 call mpas_pool_get_array(nggps_diagnostics_radiation,'tcdclcl' ,tcdclcl )

!... physics diagnostics:
 call mpas_pool_get_array(nggps_diagnostics_physics,'cpratsfc',cpratsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'pratesfc',pratesfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'gfluxsfc',gfluxsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'lhtflsfc',lhtflsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'shtflsfc',shtflsfc)
 call mpas_pool_get_array(nggps_diagnostics_physics,'ugwdsfc' ,ugwdsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'vgwdsfc' ,vgwdsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'uflxsfc' ,uflxsfc )
 call mpas_pool_get_array(nggps_diagnostics_physics,'vflxsfc' ,vflxsfc )

 call mpas_pool_get_array(nggps_diagnostics,'cnvprcp',cnvprcp)
 call mpas_pool_get_array(nggps_diagnostics,'totprcp',totprcp)
 call mpas_pool_get_array(nggps_diagnostics,'gflux'  ,gflux  )
 call mpas_pool_get_array(nggps_diagnostics,'dqsfc'  ,dqsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dtsfc'  ,dtsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dugwd'  ,dugwd  )
 call mpas_pool_get_array(nggps_diagnostics,'dvgwd'  ,dvgwd  )
 call mpas_pool_get_array(nggps_diagnostics,'dusfc'  ,dusfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dvsfc'  ,dvsfc  )

 call mpas_pool_get_array(nggps_diagnostics,'psmean' ,psmean )
 call mpas_pool_get_array(nggps_diagnostics,'suntim' ,suntim )
 call mpas_pool_get_array(nggps_diagnostics,'dlwsfc' ,dlwsfc )
 call mpas_pool_get_array(nggps_diagnostics,'ulwsfc' ,ulwsfc ) 
 call mpas_pool_get_array(nggps_diagnostics,'evbsa'  ,evbsa  )
 call mpas_pool_get_array(nggps_diagnostics,'evcwa'  ,evcwa  )
 call mpas_pool_get_array(nggps_diagnostics,'transa' ,transa )
 call mpas_pool_get_array(nggps_diagnostics,'sbsnoa' ,sbsnoa )
 call mpas_pool_get_array(nggps_diagnostics,'snowca' ,snowca )
 call mpas_pool_get_array(nggps_diagnostics,'snohfa' ,snohfa )
 call mpas_pool_get_array(nggps_diagnostics,'ep'     ,ep     )
 call mpas_pool_get_array(nggps_diagnostics,'cldwrk' ,cldwrk )
 call mpas_pool_get_array(nggps_diagnostics,'runoff' ,runoff )
 call mpas_pool_get_array(nggps_diagnostics,'srunoff',srunoff)
 call mpas_pool_get_array(nggps_diagnostics,'tmpmax' ,tmpmax )
 call mpas_pool_get_array(nggps_diagnostics,'tmpmin' ,tmpmin )
 call mpas_pool_get_array(nggps_diagnostics,'spfhmax',spfhmax)
 call mpas_pool_get_array(nggps_diagnostics,'spfhmin',spfhmin)

!... zero out radiation diagnostics time-averaged in subroutine physics_nggps_compute_diagnostics:
 do iCell = 1, nCellsSolve
    do iFx = 1, nFxr
       fluxr(iFx,iCell) = 0._RKIND
    enddo
    albdosfc(iCell) = 0._RKIND
    dlwrfsfc(iCell) = 0._RKIND
    ulwrfsfc(iCell) = 0._RKIND
    dswrfsfc(iCell) = 0._RKIND
    uswrfsfc(iCell) = 0._RKIND
    ulwrftoa(iCell) = 0._RKIND
    dswrftoa(iCell) = 0._RKIND
    uswrftoa(iCell) = 0._RKIND
    tcdcclm(iCell)  = 0._RKIND
    tcdchcl(iCell)  = 0._RKIND
    tcdcmcl(iCell)  = 0._RKIND
    tcdclcl(iCell)  = 0._RKIND
 enddo

!... zero out physics diagnostics time-averaged in subroutine physics_nggps_compute_diagnostics:
 do iCell = 1, nCellsSolve
    cnvprcp(iCell) = 0._RKIND
    totprcp(iCell) = 0._RKIND
    gflux(iCell)   = 0._RKIND
    dqsfc(iCell)   = 0._RKIND
    dtsfc(iCell)   = 0._RKIND
    dugwd(iCell)   = 0._RKIND
    dvgwd(iCell)   = 0._RKIND
    dusfc(iCell)   = 0._RKIND
    dvsfc(iCell)   = 0._RKIND

    cpratsfc(iCell) = 0._RKIND
    pratesfc(iCell) = 0._RKIND
    gfluxsfc(iCell) = 0._RKIND
    lhtflsfc(iCell) = 0._RKIND
    shtflsfc(iCell) = 0._RKIND
    ugwdsfc(iCell)  = 0._RKIND
    vgwdsfc(iCell)  = 0._RKIND
    uflxsfc(iCell)  = 0._RKIND
    vflxsfc(iCell)  = 0._RKIND
 enddo 
 
!... zero out all other physics diagnostics:
 do iCell = 1, nCellsSolve
    psmean(iCell)  = 0._RKIND
    suntim(iCell)  = 0._RKIND
    dlwsfc(iCell)  = 0._RKIND
    ulwsfc(iCell)  = 0._RKIND
    evbsa(iCell)   = 0._RKIND
    evcwa(iCell)   = 0._RKIND
    transa(iCell)  = 0._RKIND
    sbsnoa(iCell)  = 0._RKIND
    snowca(iCell)  = 0._RKIND
    snohfa(iCell)  = 0._RKIND
    ep(iCell)      = 0._RKIND
    cldwrk(iCell)  = 0._RKIND
    runoff(iCell)  = 0._RKIND
    srunoff(iCell) = 0._RKIND
    tmpmax(iCell)  = 0._RKIND
    tmpmin(iCell)  = 0._RKIND
    spfhmax(iCell) = 0._RKIND
    spfhmin(iCell) = 0._RKIND
 enddo

 write(0,*) '--- end subroutine physics_nggps_update_diagnostics'
 write(0,*)

 end subroutine physics_nggps_zero_diagnostics

!=================================================================================================================
 subroutine physics_nggps_compute_relhum(mesh,state,diag,time_lev)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state

 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: index_qv,index_qc,index_qo3

 real(kind=RKIND),dimension(:,:),pointer  :: relhum,theta_m
 real(kind=RKIND),dimension(:,:),pointer  :: exner,pressure_b,pressure_p
 real(kind=PHYSKIND),dimension(:,:),pointer  :: qv
 real(kind=PHYSKIND),dimension(:,:,:),pointer:: scalars

!local variables:
 integer:: iCell,k
 real(kind=PHYSKIND):: es,pres,qs,rs,temp

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_dimension(state,'index_qv',index_qv )

 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv  => scalars(index_qv,:,:)

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'relhum'       ,relhum    )

 if(l_diags) then

    do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       temp = theta_m(k,iCell)/(1._RKIND+R_v/R_d*qv(k,iCell))
       temp = temp*exner(k,iCell)
       pres = pressure_p(k,iCell) + pressure_b(k,iCell)
       es   = min(0.99_RKIND*pres,fpvs(temp))
       qs   = eps*es / (pres + epsm1*es)
       rs   = qs / (1._RKIND - qs)

       relhum(k,iCell) = qv(k,iCell)/rs * 100._RKIND
    enddo
    enddo
 endif
 205 format(2i9,10(1x,e15.8))

 end subroutine physics_nggps_compute_relhum

!=================================================================================================================
 end module mpas_atmphys_nggps_compute_diagnostics
!=================================================================================================================
