! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_gfs_initialize_real
 use mpas_dmpar
 use mpas_grid_types
 use mpas_kind_types
 use mpas_atmphys_date_time
 use mpas_atmphys_gfs_date_time
 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_gfs_initialize_real


!>\brief MPAS initialization of surface properties for real case initialization using GFS physics.
!>       Subroutine physics_gfs_initialize_real follows the same steps as physics_initialize_real,
!>       except that static fields needed in the GFS physics are used.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-01-27.
!>
!>\details
!>
!> subroutines in mpas_atmphys_initialize_real:
!> --------------------------------------------
!> physics_gfs_initialize_real: main subroutine (called from subroutine init_atm_setup_test_case in
!                               ./src/core_init_atmosphere/mpas_init_atm_test_cases.F).
!> init_soil_layers           : main subroutine for initialization of soil !properties. 
!> init_soil_layers_depth     : initialize height and depth of soil layers needed in NOAH scheme.
!> init_soil_layers_properties: initialize soil temperature, soil moisture, etc.
!> adjust_input_soiltemp      : adjust the deep soil temperature to sea-level values.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Converted snowh (actual snow depth over land/sea ice) from m to mm.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-04-19.
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine physics_gfs_initialize_real(configs,dims,mesh,fg,dminfo)
!=================================================================================================================
!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: mesh
 type(dm_info),intent(in)  :: dminfo

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg 

!local pointer:
 character(len=StrKIND),pointer:: config_start_time

 integer,pointer:: config_nsoillevels
 integer,pointer:: nCellsSolve

 real(kind=RKIND),dimension(:),pointer:: alnsf,alnwf,alvsf,alvwf,snow,snowc,snowh,vegfra
 real(kind=RKIND),dimension(:,:),pointer:: alnsf4s,alnwf4s,alvsf4s,alvwf4s,greenfrac

!local variables:
 character(len=StrKIND):: initial_date
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter physics_gfs_initialize_real:'

 call mpas_pool_get_config(configs,'config_start_time' ,config_start_time )
 call mpas_pool_get_config(configs,'config_nsoillevels',config_nsoillevels)

 call mpas_pool_get_dimension(dims,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(mesh,'alnsf4s'  ,alnsf4s  )
 call mpas_pool_get_array(mesh,'alnwf4s'  ,alnwf4s  )
 call mpas_pool_get_array(mesh,'alvsf4s'  ,alvsf4s  )
 call mpas_pool_get_array(mesh,'alvwf4s'  ,alvwf4s  )
 call mpas_pool_get_array(mesh,'greenfrac',greenfrac)

 call mpas_pool_get_array(fg,'alnsf' ,alnsf )
 call mpas_pool_get_array(fg,'alnwf' ,alnwf )
 call mpas_pool_get_array(fg,'alvsf' ,alvsf )
 call mpas_pool_get_array(fg,'alvwf' ,alvwf )
 call mpas_pool_get_array(fg,'snow'  ,snow  )
 call mpas_pool_get_array(fg,'snowc' ,snowc )
 call mpas_pool_get_array(fg,'snowh ',snowh )
 call mpas_pool_get_array(fg,'vegfra',vegfra)

 initial_date = trim(config_start_time)

!initialization of the green-ness (vegetation) fraction: interpolation of the monthly mean values
!to the initial date. get the min/max for each cell for the monthly green-ness fraction:
 call monthly_interp_to_date(nCellsSolve,initial_date,greenfrac,vegfra)

!initialization of the seasonal mean near-infrared and visible spectral surface albedos with strong
!and weak dependence on the solar zenith angle
 call season_interp_to_date(nCellsSolve,initial_date,alnsf4s,alnsf)
 call season_interp_to_date(nCellsSolve,initial_date,alnwf4s,alnwf)
 call season_interp_to_date(nCellsSolve,initial_date,alvsf4s,alvsf)
 call season_interp_to_date(nCellsSolve,initial_date,alvwf4s,alvwf)

!initialization of the flag indicating the presence of snow (0 or 1) and of the snow depth (m) as
!functions of the input snow water content (kg/m2). we use a 5:1 ratio from liquid water equivalent
!to snow depth:
 do iCell = 1, nCellsSolve
    if(snow(iCell) .ge. 10._RKIND) then
       snowc(iCell) = 1._RKIND
    else
       snowc(iCell) = 0._RKIND
    endif
    snowh(iCell) = snow(iCell) * 5.0_RKIND / 1000._RKIND

    !In the GFS physics (see subroutine sfc_drv.f), snowh needs to be in mm:
    snowh(iCell) = snowh(iCell) * 1000._RKIND
 enddo

!initialization of soil layers properties:
 call init_soil_layers(config_nsoillevels,dims,mesh,fg,dminfo)

 write(0,*) '--- end physics_gfs_initialize_real:'

 end subroutine physics_gfs_initialize_real

!=================================================================================================================
 subroutine init_soil_layers(config_nsoillevels,dims,mesh,fg,dminfo)
!=================================================================================================================

!input arguments:
 integer,intent(in):: config_nsoillevels
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: mesh
 type(dm_info),intent(in)  :: dminfo

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!-----------------------------------------------------------------------------------------------------------------

!adjust the annual mean deep soil temperature:
 call adjust_input_soiltemps(dims,mesh,fg)

!initialize the depth of the soil layers:
 call init_soil_layers_depth(config_nsoillevels,dims,mesh,fg)

!initialize the temperature, moisture, and liquid water of the individual soil layers:
 call init_soil_layers_properties(config_nsoillevels,dims,mesh,fg,dminfo)

 end subroutine init_soil_layers

!=================================================================================================================
 subroutine adjust_input_soiltemps(dims,mesh,fg)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: dims 
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local pointers:
 integer,pointer:: nCellsSolve,nFGSoilLevels
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer  :: soilz,ter
 real(kind=RKIND),dimension(:),pointer  :: skintemp,soiltemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: st_fg

!local variables:
 integer:: iCell,ifgSoil

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(dims,'nCellsSolve'  ,nCellsSolve  )
 call mpas_pool_get_dimension(dims,'nFGSoilLevels',nFGSoilLevels)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'soiltemp',soiltemp)
 call mpas_pool_get_array(mesh,'ter'     ,ter     )

 call mpas_pool_get_array(fg,'skintemp',skintemp)
 call mpas_pool_get_array(fg,'tmn'     ,tmn     )
 call mpas_pool_get_array(fg,'st_fg'   ,st_fg   )
 call mpas_pool_get_array(fg,'soilz'   ,soilz   )

 do iCell = 1, nCellsSolve
    if(landmask(iCell) .eq. 1) then

       !adjust the annual deep mean soil temperature and skin temperatures over land:
       tmn(iCell) = soiltemp(iCell) - 0.0065_RKIND * ter(iCell)
       skintemp(iCell) = skintemp(iCell) - 0.0065_RKIND * (ter(iCell)-soilz(iCell))

       !adjust the soil layer temperatures:
       do ifgSoil = 1, nFGSoilLevels
          st_fg(ifgSoil,iCell) = st_fg(ifgSoil,iCell) - 0.0065_RKIND * (ter(iCell)-soilz(iCell))
       enddo

    elseif(landmask(iCell) .eq. 0) then

       tmn(iCell) = skintemp(iCell)

    endif
 enddo

 end subroutine adjust_input_soiltemps

!=================================================================================================================
 subroutine init_soil_layers_depth(config_nsoillevels,dims,mesh,fg)
!=================================================================================================================

!input arguments:
 integer,intent(in):: config_nsoillevels
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg 

!local pointers:
 integer,pointer:: nCellsSolve,nSoilLevels,nFGSoilLevels

 real(kind=RKIND),dimension(:,:),pointer:: dzs_fg,zs_fg
 real(kind=RKIND),dimension(:,:),pointer:: dzs,zs

!local variables:
 integer:: iCell,iSoil

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_depth:'

 call mpas_pool_get_dimension(dims,'nCellsSolve'  ,nCellsSolve  )
 call mpas_pool_get_dimension(dims,'nSoilLevels'  ,nSoilLevels  )
 call mpas_pool_get_dimension(dims,'nFGSoilLevels',nFGSoilLevels)

 call mpas_pool_get_array(fg,'dzs_fg',dzs_fg)
 call mpas_pool_get_array(fg,'zs_fg' ,zs_fg )
 call mpas_pool_get_array(fg,'dzs'   ,dzs   )
 call mpas_pool_get_array(fg,'zs'    ,zs    )

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 do iCell = 1, nCellsSolve
    iSoil = 1
    zs_fg(iSoil,iCell) = 0.5_RKIND*dzs_fg(iSoil,iCell)
    do iSoil = 2, nFGSoilLevels
       zs_fg(iSoil,iCell) = zs_fg(iSoil-1,iCell)            &
                          + 0.5_RKIND*dzs_fg(iSoil-1,iCell) &
                          + 0.5_RKIND*dzs_fg(iSoil,iCell)
    enddo
 enddo

 do iCell = 1, nCellsSolve
    dzs(1,iCell) = 0.10_RKIND
    dzs(2,iCell) = 0.30_RKIND
    dzs(3,iCell) = 0.60_RKIND
    dzs(4,iCell) = 1.00_RKIND

    iSoil = 1
    zs(iSoil,iCell)  = 0.5_RKIND * dzs(iSoil,iCell)
    do iSoil = 2, nSoilLevels
       zs(iSoil,iCell) = zs(iSoil-1,iCell)            &
                       + 0.5_RKIND*dzs(iSoil-1,iCell) &
                       + 0.5_RKIND*dzs(iSoil,iCell)
    enddo

 enddo

 end subroutine init_soil_layers_depth

!=================================================================================================================
 subroutine init_soil_layers_properties(config_nsoillevels,dims,mesh,fg,dminfo)
!=================================================================================================================

!input arguments:
 integer,intent(in):: config_nsoillevels
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: mesh
 type(dm_info),intent(in)  :: dminfo

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local pointers:
 integer,pointer:: nCellsSolve
 integer,pointer:: nSoilLevels,nFGSoilLevels
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer  :: skintemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: dzs,zs,tslb,smois,sh2o,smcrel
 real(kind=RKIND),dimension(:,:),pointer:: sm_fg,st_fg,zs_fg

!local variables:
 integer:: iCell,ifgSoil,iSoil
 integer:: num_sm,num_st

 real(kind=RKIND),dimension(:,:),allocatable:: zhave,sm_input,st_input

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_properties:'

 call mpas_pool_get_dimension(dims,'nCellsSolve'  ,nCellsSolve  )
 call mpas_pool_get_dimension(dims,'nSoilLevels'  ,nSoilLevels  )
 call mpas_pool_get_dimension(dims,'nFGSoilLevels',nFGSoilLevels)

 call mpas_pool_get_array(mesh,'landmask',landmask)

 call mpas_pool_get_array(fg,'skintemp',skintemp)
 call mpas_pool_get_array(fg,'tmn'     ,tmn     )
 call mpas_pool_get_array(fg,'dzs'     ,dzs     )
 call mpas_pool_get_array(fg,'zs'      ,zs      )
 call mpas_pool_get_array(fg,'tslb'    ,tslb    )
 call mpas_pool_get_array(fg,'smois'   ,smois   )
 call mpas_pool_get_array(fg,'sh2o'    ,sh2o    )
 call mpas_pool_get_array(fg ,'smcrel' ,smcrel  )
 call mpas_pool_get_array(fg,'sm_fg'   ,sm_fg   )
 call mpas_pool_get_array(fg,'st_fg'   ,st_fg   )
 call mpas_pool_get_array(fg,'zs_fg'   ,zs_fg   )

 write(0,*) 'nSoilLevels   =',nSoilLevels
 write(0,*) 'nFGSoilLevels =',nFGSoilLevels

!check that interpolation of the meteorological data to the MPAS grid did not create negative
!values for the first-guess soil temperatures and soil moistures.
 num_sm = 0
 num_st = 0
 do iCell = 1, nCellsSolve
    do ifgSoil = 1, nFGSoilLevels
       if(st_fg(ifgSoil,iCell) .le. 0._RKIND) num_st = num_st + 1
       if(sm_fg(ifgSoil,iCell) .lt. 0._RKIND) num_sm = num_sm + 1
    enddo
 enddo
 if(num_st .gt. 0) then
    write(0,*) 'Error in interpolation of st_fg to MPAS grid: num_st =', num_st
    call mpas_dmpar_abort(dminfo)
 elseif(num_sm .gt. 0) then
    write(0,*) 'Error in interpolation of sm_fg to MPAS grid: num_sm =', num_sm
    call mpas_dmpar_abort(dminfo)
 endif

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')


 if(.not.allocated(zhave)   ) allocate(zhave(nFGSoilLevels+2,nCellsSolve)   )
 if(.not.allocated(st_input)) allocate(st_input(nFGSoilLevels+2,nCellsSolve))
 if(.not.allocated(sm_input)) allocate(sm_input(nFGSoilLevels+2,nCellsSolve))

 do iCell = 1, nCellsSolve

    ifgSoil = 1
    zhave(ifgSoil,iCell)    = 0._RKIND
    st_input(ifgSoil,iCell) = skintemp(iCell)
    sm_input(ifgSoil,iCell) = sm_fg(ifgSoil+1,iCell)

    do ifgSoil = 1, nFGSoilLevels
       zhave(ifgSoil+1,iCell) = zs_fg(ifgSoil,iCell) / 100._RKIND
       st_input(ifgSoil+1,iCell) = st_fg(ifgSoil,iCell)
       sm_input(ifgSoil+1,iCell) = sm_fg(ifgSoil,iCell)
    enddo

    zhave(nFGSoilLevels+2,iCell) = 300._RKIND/100._RKIND
    st_input(nFGSoilLevels+2,iCell) = tmn(iCell)
    sm_input(nFGSoilLevels+2,iCell) = sm_input(nFGSoilLevels,iCell)

    if(iCell .eq. 1) then
       do ifgSoil = 1,nFGSoilLevels+2
          write(0,101) ifgSoil,zhave(ifgSoil,iCell)
       enddo
    endif

 enddo

!... interpolate the soil temperature, soil moisture, and soil liquid temperature to the four
!    layers used in the NOAH land surface scheme:

 do iCell = 1, nCellsSolve

    if(landmask(iCell) .eq. 1) then

       noah: do iSoil = 1 , nSoilLevels
          input: do ifgSoil = 1 , nFGSoilLevels+2-1
             if(iCell .eq. 1) write(0,102) iSoil,ifgSoil,zs(iSoil,iCell), &
                zhave(ifgSoil,iCell),zhave(ifgSoil+1,iCell)

             if(zs(iSoil,iCell).ge.zhave(ifgSoil,iCell) .and. &
                zs(iSoil,iCell).le.zhave(ifgSoil+1,iCell)) then

                tslb(iSoil,iCell) = &
                      (st_input(ifgSoil,iCell) * (zhave(ifgSoil+1,iCell)-zs(iSoil,iCell))    &
                    +  st_input(ifgSoil+1,iCell) * (zs(iSoil,iCell)-zhave(ifgSoil,iCell)))   &
                         / (zhave(ifgSoil+1,iCell)-zhave(ifgSoil,iCell))
                if(iCell .eq. 1) write(6,102) iSoil,ifgSoil,zs(iSoil,iCell), &
                   zhave(ifgSoil,iCell),zhave(ifgSoil+1,iCell)
                         
                smois(iSoil,iCell) = &
                       (sm_input(ifgSoil,iCell) * (zhave(ifgSoil+1,iCell)-zs(iSoil,iCell))   &
                    +  sm_input(ifgSoil+1,iCell) * (zs(iSoil,iCell)-zhave(ifgSoil,iCell)))   &
                    / (zhave(ifgSoil+1,iCell)-zhave(ifgSoil,iCell))

                sh2o(iSoil,iCell)   = 0._RKIND
                smcrel(iSoil,iCell) = 0._RKIND

                exit input
             endif
          enddo input
          if(iCell.eq. 1) write(0,*)
       enddo noah

    elseif(landmask(iCell) .eq. 0) then

       !fill the soil temperatures with the skin temperatures over oceans:
       do iSoil = 1, nSoilLevels
          tslb(iSoil,iCell)    = skintemp(iCell)
          smois(iSoil,iCell)   = 1._RKIND
          sh2o(iSoil,iCell)    = 1._RKIND
          smcrel(iSoil,iCell)  = 0._RKIND
       enddo

    endif

 enddo

!... final checks:

 do iCell = 1, nCellsSolve

    if(landmask(iCell).eq.1 .and. tslb(1,iCell).gt.170._RKIND .and. tslb(1,iCell).lt.400._RKIND &
       .and. smois(1,iCell).lt.0.005_RKIND) then
       do iSoil = 1, nSoilLevels
          smois(iSoil,iCell) = 0.005_RKIND
       enddo
    endif

 enddo

!formats:
 101 format(i4,4(1x,e15.8))
 102 format(2i5,5(1x,e15.8))

 if(allocated(zhave)   ) deallocate(zhave )
 if(allocated(st_input)) deallocate(st_input)
 if(allocated(sm_input)) deallocate(sm_input)

 end subroutine init_soil_layers_properties

!=================================================================================================================
 end module mpas_atmphys_gfs_initialize_real
!=================================================================================================================

