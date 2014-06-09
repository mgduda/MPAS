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
 subroutine physics_gfs_initialize_real(mesh,fg,dminfo)
!=================================================================================================================
!input arguments:
 type(mesh_type),intent(in):: mesh
 type(dm_info),intent(in)  :: dminfo

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 character(len=StrKIND):: initial_date
 integer:: iCell,nCellsSolve

 real(kind=RKIND),dimension(:),pointer:: alnsf,alnwf,alvsf,alvwf,snow,snowc,snowh,vegfra
 real(kind=RKIND),dimension(:,:),pointer:: alnsf4s,alnwf4s,alvsf4s,alvwf4s,greenfrac

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter physics_gfs_initialize_real:'

 nCellsSolve = mesh % nCellsSolve

 alnsf4s   => mesh % alnsf4s   % array
 alnwf4s   => mesh % alnwf4s   % array
 alvsf4s   => mesh % alvsf4s   % array
 alvwf4s   => mesh % alvwf4s   % array
 greenfrac => mesh % greenfrac % array

 alnsf  => fg % alnsf  % array
 alnwf  => fg % alnwf  % array
 alvsf  => fg % alvsf  % array
 alvwf  => fg % alvwf  % array
 snow   => fg % snow   % array
 snowc  => fg % snowc  % array
 snowh  => fg % snowh  % array
 vegfra => fg % vegfra % array

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
 call init_soil_layers(mesh,fg,dminfo)

 write(0,*) '--- end physics_gfs_initialize_real:'

 end subroutine physics_gfs_initialize_real

!=================================================================================================================
 subroutine init_soil_layers(mesh,fg,dminfo)
!=================================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(dm_info),intent(in)  :: dminfo

!inout arguments:
 type(fg_type),intent(inout):: fg

!-----------------------------------------------------------------------------------------------------------------

!adjust the annual mean deep soil temperature:
 call adjust_input_soiltemps(mesh,fg)

!initialize the depth of the soil layers:
 call init_soil_layers_depth(mesh,fg)

!initialize the temperature, moisture, and liquid water of the individual soil layers:
 call init_soil_layers_properties(mesh,fg,dminfo)

 end subroutine init_soil_layers

!=================================================================================================================
 subroutine adjust_input_soiltemps(mesh,fg)
!=================================================================================================================

!input arguments:
 type(mesh_type),intent(in) :: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg

!local variables:
 integer:: iCell,ifgSoil
 integer:: nCellsSolve,nFGSoilLevels
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer  :: soilz,ter
 real(kind=RKIND),dimension(:),pointer  :: skintemp,soiltemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: st_fg

!-----------------------------------------------------------------------------------------------------------------

 nCellsSolve   = mesh % nCellsSolve
 nFGSoilLevels = mesh % nFGSoilLevels

 landmask => mesh % landmask % array
 soiltemp => mesh % soiltemp % array
 ter      => mesh % ter      % array

 skintemp => fg % skintemp % array
 tmn      => fg % tmn      % array
 st_fg    => fg % st_fg    % array
 soilz    => fg % soilz    % array

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
 subroutine init_soil_layers_depth(mesh,fg)
!=================================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(fg_type),intent(inout):: fg 

!local variables:
 integer:: iCell,iSoil
 integer:: nCellsSolve,nSoilLevels,nFGSoilLevels

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_depth:'

 nCellsSolve   = mesh % nCellsSolve
 nSoilLevels   = mesh % nSoilLevels
 nFGSoilLevels = mesh % nFGSoilLevels

 if(config_nsoillevels .ne. 4) &
    call physics_error_fatal('NOAH lsm uses 4 soil layers. Correct config_nsoillevels.')

 do iCell = 1, nCellsSolve

    iSoil = 1
    fg % zs_fg % array(iSoil,iCell) = 0.5_RKIND * fg % dzs_fg % array(iSoil,iCell)
    do iSoil = 2, nFGSoilLevels
       fg % zs_fg % array(iSoil,iCell) = fg % zs_fg % array(iSoil-1,iCell)        &
                                 + 0.5_RKIND * fg % dzs_fg % array(iSoil-1,iCell) &
                                 + 0.5_RKIND * fg % dzs_fg % array(iSoil,iCell)
    enddo

 enddo

 do iCell = 1, nCellsSolve
    fg % dzs % array(1,iCell) = 0.10_RKIND
    fg % dzs % array(2,iCell) = 0.30_RKIND
    fg % dzs % array(3,iCell) = 0.60_RKIND
    fg % dzs % array(4,iCell) = 1.00_RKIND

    iSoil = 1
    fg % zs % array(iSoil,iCell)  = 0.5_RKIND * fg % dzs % array(iSoil,iCell)
    do iSoil = 2, nSoilLevels
       fg % zs % array(iSoil,iCell) = fg % zs % array(iSoil-1,iCell)              &
                                    + 0.5_RKIND * fg % dzs % array(iSoil-1,iCell) &
                                    + 0.5_RKIND * fg % dzs % array(iSoil,iCell)
    enddo

 enddo

 end subroutine init_soil_layers_depth

!=================================================================================================================
 subroutine init_soil_layers_properties(mesh,fg,dminfo)
!=================================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(dm_info),intent(in)  :: dminfo

!inout arguments:
 type(fg_type),intent(inout):: fg

!local variables:
 integer:: iCell,ifgSoil,iSoil
 integer:: nCellsSolve,nFGSoilLevels,nSoilLevels
 integer:: num_sm,num_st
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:,:),allocatable:: zhave,sm_input,st_input

 real(kind=RKIND),dimension(:),pointer  :: skintemp,tmn
 real(kind=RKIND),dimension(:,:),pointer:: dzs,zs,tslb,smois,sh2o,smcrel
 real(kind=RKIND),dimension(:,:),pointer:: sm_fg,st_fg,zs_fg

!-----------------------------------------------------------------------------------------------------------------

!write(0,*)
 write(0,*) '--- enter subroutine init_soil_layers_properties:'

 nCellsSolve   = mesh % nCellsSolve
 nSoilLevels   = mesh % nSoilLevels
 nFGSoilLevels = mesh % nFGSoilLevels
 write(0,*) 'nSoilLevels   =',nSoilLevels
 write(0,*) 'nFGSoilLevels =',nFGSoilLevels

 landmask => mesh % landmask % array

 zs_fg => fg % zs_fg % array
 st_fg => fg % st_fg % array
 sm_fg => fg % sm_fg % array

 zs     => fg % zs  % array
 dzs    => fg % dzs % array
 sh2o   => fg % sh2o    % array
 smcrel => fg % smcrel % array
 smois  => fg % smois   % array
 tslb   => fg % tslb    % array
 skintemp => fg % skintemp % array
 tmn      => fg % tmn      % array

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

