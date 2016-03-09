! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_initialize_real
 use mpas_kind_types
 use mpas_dmpar
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_date_time
 use mpas_atmphys_nggps_date_time
 use mpas_atmphys_utilities

 implicit none
 private
 public:: physics_nggps_initialize_real


 contains


!=================================================================================================================
 subroutine physics_nggps_initialize_real(mesh,fg,dminfo,dims,configs)
!=================================================================================================================
!input arguments:
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(dm_info),intent(in):: dminfo

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg 

!local variables:
 character(len=StrKIND):: initial_date
 integer:: iCell

 character(len=StrKIND),pointer:: config_start_time
 integer,pointer:: nCellsSolve,nCells

 integer,dimension(:),pointer:: landmask
 real(kind=PHYSKIND),dimension(:),pointer:: alnsf,alnwf,alvsf,alvwf,snow,snowc,vegfra
 real(kind=PHYSKIND),dimension(:,:),pointer:: alnsf4s,alnwf4s,alvsf4s,alvwf4s,greenfrac

 real(kind=PHYSKIND),dimension(:),pointer:: facsf,facwf
 real(kind=PHYSKIND),dimension(:),pointer:: canopy,sltype,vgtype
 real(kind=PHYSKIND),dimension(:),pointer:: shdmax,shdmin,slope,snoalb

 real(kind=PHYSKIND):: sfc_max,sfc_min
 real(kind=PHYSKIND),dimension(:),allocatable:: sfc

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter physics_nggps_initialize_real:'

 call mpas_pool_get_config(configs,'config_start_time',config_start_time)
 call mpas_pool_get_dimension(dims,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(dims,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)

 call mpas_pool_get_array(fg,'snow'  ,snow  )
 call mpas_pool_get_array(fg,'snowc' ,snowc )

!NGGPS PROJECT: when reading intermediate surface fields, we do not need to update alnsf, alnwf, alvsf, alvwf,
!vegfra using the static fields alnsf4s, alnwf4s, alvsf4s, alvwf4s, and greenfrac:
!Laura D. Fowler (2015-10-19).

!call mpas_pool_get_array(mesh,'alnsf4s'  ,alnsf4s  )
!call mpas_pool_get_array(mesh,'alnwf4s'  ,alnwf4s  )
!call mpas_pool_get_array(mesh,'alvsf4s'  ,alvsf4s  )
!call mpas_pool_get_array(mesh,'alvwf4s'  ,alvwf4s  )
!call mpas_pool_get_array(mesh,'greenfrac',greenfrac)

 call mpas_pool_get_array(fg,'alnsf' ,alnsf )
 call mpas_pool_get_array(fg,'alnwf' ,alnwf )
 call mpas_pool_get_array(fg,'alvsf' ,alvsf )
 call mpas_pool_get_array(fg,'alvwf' ,alvwf )
 call mpas_pool_get_array(fg,'vegfra',vegfra)

 do iCell = 1, nCellsSolve
    if(landmask(iCell) == 0) then
       if(alnsf(iCell) > 0.06_PHYSKIND) then
          write(0,*) '--- alnsf greater than 0.06 over oceans'
          write(0,*) iCell,alnsf(iCell)
          stop
       elseif(alnwf(iCell) > 0.06_PHYSKIND) then
          write(0,*) '--- alnwf greater than 0.06 over oceans'
          write(0,*) iCell,alnwf(iCell)
          stop
       elseif(alvsf(iCell) > 0.06_PHYSKIND) then
          write(0,*) '--- alvsf greater than 0.06 over oceans'
          write(0,*) iCell,alvsf(iCell)
          stop
       elseif(alvwf(iCell) > 0.06_PHYSKIND) then
          write(0,*) '--- alvwf greater than 0.06 over oceans'
          write(0,*) iCell,alvwf(iCell)
          stop
       elseif(vegfra(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- vegfra greater than 0. over oceans'
          write(0,*) iCell,vegfra(iCell)
          stop
       endif
    endif
 enddo

!initial_date = trim(config_start_time)

!initialization of the green-ness (vegetation) fraction: interpolation of the monthly mean values
!to the initial date. get the min/max for each cell for the monthly green-ness fraction:
!call monthly_interp_to_date(nCellsSolve,initial_date,greenfrac,vegfra)

!initialization of the seasonal mean near-infrared and visible spectral surface albedos with strong
!and weak dependence on the solar zenith angle
!call season_interp_to_date(nCellsSolve,initial_date,alnsf4s,alnsf)
!call season_interp_to_date(nCellsSolve,initial_date,alnwf4s,alnwf)
!call season_interp_to_date(nCellsSolve,initial_date,alvsf4s,alvsf)
!call season_interp_to_date(nCellsSolve,initial_date,alvwf4s,alvwf)

!additional checks after initialization of surface fields using intermediate surface data:
 call mpas_pool_get_array(mesh,'sltype',sltype)
 call mpas_pool_get_array(mesh,'vgtype',vgtype)
 call mpas_pool_get_array(mesh,'facsf' ,facsf )
 call mpas_pool_get_array(mesh,'facwf' ,facwf )
 call mpas_pool_get_array(mesh,'shdmax',shdmax)
 call mpas_pool_get_array(mesh,'shdmin',shdmin)
 call mpas_pool_get_array(mesh,'slope' ,slope )
 call mpas_pool_get_array(mesh,'snoalb',snoalb)
 call mpas_pool_get_array(fg  ,'canopy',canopy)

 do iCell = 1, nCellsSolve
    if(landmask(iCell) == 0) then
       if(vgtype(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- vgtype greater than 0. over oceans'
          write(0,*) iCell,vgtype(iCell)
          stop
       elseif(sltype(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- sltype greater than 0. over oceans'
          write(0,*) iCell,sltype(iCell)
          stop
       elseif(canopy(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- canopy greater than 0. over oceans'
          write(0,*) iCell,canopy(iCell)
          stop
       elseif(facsf(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- facsf greater than 0. over oceans'
          write(0,*) iCell,facsf(iCell)
          stop
       elseif(facwf(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- facwf greater than 0. over oceans'
          write(0,*) iCell,facwf(iCell)
          stop
       elseif(shdmax(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- shdmax greater than 0. over oceans'
          write(0,*) iCell,shdmax(iCell)
          stop
       elseif(shdmin(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- shdmin greater than 0. over oceans'
          write(0,*) iCell,shdmin(iCell)
          stop
       elseif(snoalb(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- snoalb greater than 0. over oceans'
          write(0,*) iCell,snoalb(iCell)
          stop
       elseif(slope(iCell) > 0.0_PHYSKIND) then
          write(0,*) '--- slope greater than 0. over oceans'
          write(0,*) iCell,slope(iCell)
          stop
        endif
    endif
 enddo

 allocate(sfc(1:nCells+1))
 sfc_max = 0._PHYSKIND
 sfc_min = 0._PHYSKIND
 where(landmask == 0) sfc = alnsf
 write(0,*) minval(sfc)
 write(0,*) maxval(sfc) 
 deallocate(sfc)

!END NGGPS PROJECT.

!initializes the "snow flag" as a function of the accumulated water equivalent snow depth (snow):
 do iCell = 1, nCellsSolve
    snowc(iCell) = 0._PHYSKIND
    if(snow(iCell) > 0._PHYSKIND) snowc(iCell) = 1._PHYSKIND
 enddo

!initialization of soil layers properties:
 call init_soil_layers(mesh,fg,dminfo,dims,configs)

 write(0,*) '--- end physics_nggps_initialize_real:'

 end subroutine physics_nggps_initialize_real

!=================================================================================================================
 subroutine init_soil_layers(mesh,fg,dminfo,dims,configs)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!-----------------------------------------------------------------------------------------------------------------

!adjust the annual mean deep soil temperature:
 call adjust_input_soiltemps(mesh,fg,dims)

!initialize the depth of the soil layers:
 call init_soil_layers_depth(mesh,fg,dims,configs)

!initialize the temperature, moisture, and liquid water of the individual soil layers:
!call init_soil_layers_properties(mesh,fg,dminfo,dims,configs)
 call init_soil_layers_properties_nggps(mesh,fg,dminfo,dims,configs)

 end subroutine init_soil_layers

!=================================================================================================================
 subroutine adjust_input_soiltemps(mesh,fg,dims)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: dims

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local variables:
 integer:: iCell,ifgSoil
 integer,pointer:: nCellsSolve,nFGSoilLevels
 integer,dimension(:),pointer:: landmask

 real(kind=PHYSKIND),dimension(:),pointer  :: soilz,ter
 real(kind=PHYSKIND),dimension(:),pointer  :: skintemp,soiltemp,tmn
 real(kind=PHYSKIND),dimension(:,:),pointer:: st_fg

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine adjust_input_soiltemp:'

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
       tmn(iCell) = soiltemp(iCell) - 0.0065_PHYSKIND * ter(iCell)
       skintemp(iCell) = skintemp(iCell) - 0.0065_PHYSKIND * (ter(iCell)-soilz(iCell))

       !adjust the soil layer temperatures:
       do ifgSoil = 1, nFGSoilLevels
          st_fg(ifgSoil,iCell) = st_fg(ifgSoil,iCell) - 0.0065_PHYSKIND * (ter(iCell)-soilz(iCell))
       enddo

    elseif(landmask(iCell) .eq. 0) then

       tmn(iCell) = skintemp(iCell)

       !--- when initialize MPAS with sigma-level and surface data, soil temperatures are set to the surface
       !temperature over oceans; in subroutine init_atm_case_gfs, soil temperatures are initially set to 285.
       do ifgSoil = 1, nFGSoilLevels
          st_fg(ifgSoil,iCell) = skintemp(iCell)
       enddo

    endif
 enddo

 write(0,*) '--- enter subroutine adjust_input_soiltemp'

 end subroutine adjust_input_soiltemps

!=================================================================================================================
 subroutine init_soil_layers_depth(mesh,fg,dims,configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in) :: mesh
 type(mpas_pool_type),intent(in) :: dims
 type (mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local variables:
 integer:: iCell,iSoil
 integer,pointer:: nCellsSolve,nSoilLevels,nFGSoilLevels
 integer,pointer:: config_nsoillevels

 real(kind=PHYSKIND),dimension(:,:),pointer:: dzs_fg,zs_fg
 real(kind=PHYSKIND),dimension(:,:),pointer:: dzs,zs

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_depth:'

 call mpas_pool_get_dimension(dims,'nCellsSolve'  ,nCellsSolve  )
 call mpas_pool_get_dimension(dims,'nSoilLevels'  ,nSoilLevels  )
 call mpas_pool_get_dimension(dims,'nFGSoilLevels',nFGSoilLevels)

 call mpas_pool_get_array(fg,'zs_fg' , zs_fg )
 call mpas_pool_get_array(fg,'dzs_fg', dzs_fg)
 call mpas_pool_get_array(fg,'zs'    , zs    )
 call mpas_pool_get_array(fg,'dzs'   , dzs   )

 call mpas_pool_get_config(configs,'config_nsoillevels',config_nsoillevels)

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 do iCell = 1, nCellsSolve

    iSoil = 1
    zs_fg(iSoil,iCell) = 0.5_PHYSKIND * dzs_fg(iSoil,iCell)
    do iSoil = 2, nFGSoilLevels
       zs_fg(iSoil,iCell) = zs_fg(iSoil-1,iCell)        &
                          + 0.5_PHYSKIND * dzs_fg(iSoil-1,iCell) &
                          + 0.5_PHYSKIND * dzs_fg(iSoil,iCell)
    enddo

 enddo

 do iCell = 1, nCellsSolve
    dzs(1,iCell) = 0.10_PHYSKIND
    dzs(2,iCell) = 0.30_PHYSKIND
    dzs(3,iCell) = 0.60_PHYSKIND
    dzs(4,iCell) = 1.00_PHYSKIND

    iSoil = 1
    zs(iSoil,iCell) = 0.5_PHYSKIND * dzs(iSoil,iCell)
    do iSoil = 2, nSoilLevels
       zs(iSoil,iCell) = zs(iSoil-1,iCell)              &
                       + 0.5_PHYSKIND * dzs(iSoil-1,iCell) &
                       + 0.5_PHYSKIND * dzs(iSoil,iCell)
    enddo

 enddo

 end subroutine init_soil_layers_depth

!=================================================================================================================
 subroutine init_soil_layers_properties_nggps(mesh,fg,dminfo,dims,configs)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in) :: mesh
 type(mpas_pool_type),intent(in) :: dims
 type (mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local variables:
 integer:: iCell,ifgSoil,iSoil
 integer:: num_sl,num_sm,num_st
 integer,pointer:: nCellsSolve,nFGSoilLevels,nSoilLevels
 integer,pointer:: config_nsoillevels

 real(kind=PHYSKIND),dimension(:,:),pointer:: sh2o,smois,tslb
 real(kind=PHYSKIND),dimension(:,:),pointer:: sl_fg,sm_fg,st_fg

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_properties_nggps:'

 call mpas_pool_get_dimension(dims,'nCellsSolve'  ,nCellsSolve  )
 call mpas_pool_get_dimension(dims,'nSoilLevels'  , nSoilLevels )
 call mpas_pool_get_dimension(dims,'nFGSoilLevels',nFGSoilLevels)

 call mpas_pool_get_array(fg,'sl_fg',sl_fg)
 call mpas_pool_get_array(fg,'sm_fg',sm_fg)
 call mpas_pool_get_array(fg,'st_fg',st_fg)

 call mpas_pool_get_array(fg,'sh2o' ,sh2o )
 call mpas_pool_get_array(fg,'smois',smois)
 call mpas_pool_get_array(fg,'tslb' ,tslb )

 call mpas_pool_get_config(configs, 'config_nsoillevels', config_nsoillevels)

 write(0,*) 'nSoilLevels   =',nSoilLevels
 write(0,*) 'nFGSoilLevels =',nFGSoilLevels

 if(nSoilLevels /= nFGSoilLevels) then
    write(0,*) 'Error in number of soil layers. Set nSoilLevels to nFGSoilLevels.'
    call mpas_dmpar_abort(dminfo)
 endif

!check that interpolation of the meteorological data to the MPAS grid did not create negative
!values for the first-guess soil temperatures and soil moistures.
 num_sl = 0
 num_sm = 0
 num_st = 0
 do iCell = 1, nCellsSolve
    do ifgSoil = 1, nFGSoilLevels
       if(sl_fg(ifgSoil,iCell) .lt. 0._PHYSKIND) num_sl = num_sl + 1
       if(sm_fg(ifgSoil,iCell) .lt. 0._PHYSKIND) num_sm = num_sm + 1
       if(st_fg(ifgSoil,iCell) .le. 0._PHYSKIND) num_st = num_st + 1
    enddo
 enddo
 if(num_sl .gt. 0) then
    write(0,*) 'Error in interpolation of sl_fg to MPAS grid: num_sl =', num_sl
    call mpas_dmpar_abort(dminfo)
 elseif(num_sm .gt. 0) then
    write(0,*) 'Error in interpolation of sm_fg to MPAS grid: num_sm =', num_sm
    call mpas_dmpar_abort(dminfo)
 elseif(num_st .gt. 0) then
    write(0,*) 'Error in interpolation of st_fg to MPAS grid: num_st =', num_st
    call mpas_dmpar_abort(dminfo)
 endif

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 do iCell = 1, nCellsSolve
    do iSoil = 1, nSoilLevels
       sh2o(iSoil,iCell)  = sl_fg(iSoil,iCell)
       smois(iSoil,iCell) = sm_fg(iSoil,iCell)
       tslb(iSoil,iCell)  = st_fg(iSoil,iCell)
    enddo
 enddo

 write(0,*) '--- end subroutine init_soil_layers_properties_nggps'
 
!formats:
 101 format(i4,4(1x,e15.8))
 102 format(2i5,5(1x,e15.8))

 end subroutine init_soil_layers_properties_nggps

!=================================================================================================================
 subroutine init_soil_layers_properties(mesh,fg,dminfo,dims,configs)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in) :: mesh
 type(mpas_pool_type),intent(in) :: dims
 type (mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local variables:
 integer:: iCell,ifgSoil,iSoil
 integer:: num_sm,num_st
 integer,pointer:: nCellsSolve,nFGSoilLevels,nSoilLevels
 integer,pointer:: config_nsoillevels
 integer,dimension(:),pointer:: landmask

 real(kind=PHYSKIND),dimension(:,:),allocatable:: zhave,sm_input,st_input

 real(kind=PHYSKIND),dimension(:),pointer  :: skintemp,tmn
 real(kind=PHYSKIND),dimension(:,:),pointer:: dzs,zs,tslb,smois,sh2o,smcrel
 real(kind=PHYSKIND),dimension(:,:),pointer:: sm_fg,st_fg,zs_fg

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_properties:'

 call mpas_pool_get_dimension(dims, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(dims, 'nSoilLevels', nSoilLevels)
 call mpas_pool_get_dimension(dims, 'nFGSoilLevels', nFGSoilLevels)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(fg,'zs_fg',zs_fg)
 call mpas_pool_get_array(fg,'st_fg',st_fg)
 call mpas_pool_get_array(fg,'sm_fg',sm_fg)
 call mpas_pool_get_array(fg,'zs'      ,zs      )
 call mpas_pool_get_array(fg,'dzs'     ,dzs     )
 call mpas_pool_get_array(fg,'sh2o'    ,sh2o    )
 call mpas_pool_get_array(fg,'smcrel'  ,smcrel  )
 call mpas_pool_get_array(fg,'smois'   ,smois   )
 call mpas_pool_get_array(fg,'tslb'    ,tslb    )
 call mpas_pool_get_array(fg,'skintemp',skintemp)
 call mpas_pool_get_array(fg,'tmn'     ,tmn     )

 call mpas_pool_get_config(configs, 'config_nsoillevels', config_nsoillevels)

 write(0,*) 'nSoilLevels   =',nSoilLevels
 write(0,*) 'nFGSoilLevels =',nFGSoilLevels

!check that interpolation of the meteorological data to the MPAS grid did not create negative
!values for the first-guess soil temperatures and soil moistures.
 num_sm = 0
 num_st = 0
 do iCell = 1, nCellsSolve
    do ifgSoil = 1, nFGSoilLevels
       if(st_fg(ifgSoil,iCell) .le. 0._PHYSKIND) num_st = num_st + 1
       if(sm_fg(ifgSoil,iCell) .lt. 0._PHYSKIND) num_sm = num_sm + 1
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
    zhave(ifgSoil,iCell)    = 0._PHYSKIND
    st_input(ifgSoil,iCell) = skintemp(iCell)
    sm_input(ifgSoil,iCell) = sm_fg(ifgSoil+1,iCell)

    do ifgSoil = 1, nFGSoilLevels
       zhave(ifgSoil+1,iCell) = zs_fg(ifgSoil,iCell) / 100._PHYSKIND
       st_input(ifgSoil+1,iCell) = st_fg(ifgSoil,iCell)
       sm_input(ifgSoil+1,iCell) = sm_fg(ifgSoil,iCell)
    enddo

    zhave(nFGSoilLevels+2,iCell) = 300._PHYSKIND/100._PHYSKIND
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

                sh2o(iSoil,iCell)   = 0._PHYSKIND
                smcrel(iSoil,iCell) = 0._PHYSKIND

                exit input
             endif
          enddo input
          if(iCell.eq. 1) write(0,*)
       enddo noah

    elseif(landmask(iCell) .eq. 0) then

       !fill the soil temperatures with the skin temperatures over oceans:
       do iSoil = 1, nSoilLevels
          tslb(iSoil,iCell)    = skintemp(iCell)
          smois(iSoil,iCell)   = 1._PHYSKIND
          sh2o(iSoil,iCell)    = 1._PHYSKIND
          smcrel(iSoil,iCell)  = 0._PHYSKIND
       enddo

    endif

 enddo

!... final checks:

 do iCell = 1, nCellsSolve

    if(landmask(iCell).eq.1 .and. tslb(1,iCell).gt.170._PHYSKIND .and. tslb(1,iCell).lt.400._PHYSKIND &
       .and. smois(1,iCell).lt.0.005_PHYSKIND) then
       do iSoil = 1, nSoilLevels
          smois(iSoil,iCell) = 0.005_PHYSKIND
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
 end module mpas_atmphys_nggps_initialize_real
!=================================================================================================================

