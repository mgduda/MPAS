! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_gfs_vars
 use mpas_kind_types
 
 use date_def
 use gfsmisc_def
 use layout_grid_tracers
 use mpi_def, only: liope
 use namelist_def
 use namelist_soilveg
 use resol_def
 use ozne_def

 use module_radlw_cntr_para
 use module_radsw_cntr_para

 implicit none
 public
 save

!>\brief contains all local variables and arrays used in the GFS physics parameterizations.
!>\author Laura D. Fowler (send comments to laura@mmm.ucar.edu).
!>\date 2013-11-21.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Removed non-needed variables.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-20.


!==================================================================================================

!... parameters that should be defined using an input namelist file:
 integer, parameter:: &
    icld    = 1       !cloud computation method flag.
                      !icld   = 0: model use diagnostic cloud scheme.
                      !icld   = 1: model use prognostic cloud scheme.
 integer, parameter:: &
    np3d    = 4       !cloud parameterization scheme in parameterization of shortwave radiation.
                      !np3    = 3: model use Ferrier cloud microphysics scheme.
                      !np3    = 4: model use Zhao/Carr/Sundqvist cloud microphysics scheme.

!... fixed parameters:
 integer, parameter:: &
    iflip   = 1       !control flag for direction of vertical index.
                      !iflip = 0: index is from top-of-the-atmosphere to surface.
                      !iflip = 1: index is from surface to top-of-the-atmosphere.
 integer, parameter:: &
    irad1st = 0       !control the call to subroutines needed to initialize the radiation schemes.
                      !by setting irad1st to zero, we force the GFS initialization to always skip
                      !that step. 
 integer, parameter:: &
    isubclw = 0       !use standard grid-average clouds, no sub-column clouds approximation (not
                      !available in that parameterization of longwave radiation.
 integer, parameter:: &
    isubcsw = 0       !use standard grid-average clouds, no sub-column clouds approximation (not
                      !available in that parameterization of shortwave radiation.

!... parameters set for gbphys:
!integer, parameter:: &
!   jcap      = 1       !nb pf spectral wave truncation used only by sascnv and shalcnv.

 real(kind=RKIND), parameter:: &
    xkzm_h    = 1.,    &!background vertical diffusion for heat, q
    xkzm_m    = 1.,    &!background vertical diffusion for momentum
    xkzm_s    = 1.      !background vertical diffusion for momentum diffusion

 logical, parameter:: &
    lssav_cc   = .false. !logical flag for save data for ocean coupling.

!==================================================================================================
!INTERFACE BETWEEN MPAS AND THE GFS PHYSICS:
!==================================================================================================

 logical:: l_radtlw                   !controls call to longwave radiation parameterization.
 logical:: l_radtsw                   !controls call to shortwave radiation parameterization.

 integer,public:: num_months          !number of months.
 integer,public:: num_soils           !number of soil layers                                    [-]
 integer,public:: n_microp            !number of cloud microphysics sub timesteps

 real(kind=RKIND),public:: dt_dyn     !time-step for dynamics
 real(kind=RKIND),public:: dt_microp  !time-step for cloud microphysics parameterization.
 real(kind=RKIND),public:: dt_radtlw  !time-step for longwave radiation parameterization  [mns]
 real(kind=RKIND),public:: dt_radtsw  !time-step for shortwave radiation parameterization [mns]

!... dummy array:
 real(kind=RKIND),dimension(:),allocatable:: &
    dumphys            !dummy array

!... surface pressure:
 real(kind=RKIND),dimension(:),allocatable:: &
    pgr,              &!surface pressure                                                       [Pa]
    pgr_hyd            !hydrostatic surface pressure                                           [Pa]

!... weights:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    fzm,              &!weight for interpolation to w points                                    [-]
    fzp                !weight for interpolation to w points                                    [-]

!... zonal and meridional wind, vertical velocity:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    ugrs,             &!u-velocity interpolated to theta points                               [m/s]
    vgrs,             &!v-velocity interpolated to theta points                               [m/s]
    vvel               !vertical velocity                                                     [m/s]
    
!... random number generation:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    rann               !random number                                                           [-]  

!... layer mean temperature, layer mean pressure, layer mean exner function, and layer mean
!... geopotential height:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    prsl,             &!pressure                                                               [Pa]
    prslk,            &!Exner function                                                          [-]
    phil,             &!geopotential height of layer                                            [m]
    tgrs,             &!temperature                                                             [K]
    rhos,             &!air density                                                         [kg/m3]
    dzgrid             !thickness of layer                                                      [m]
 real(kind=RKIND),dimension(:,:),allocatable:: &
    prsl_hyd           !hydrostatic pressure                                                   [Pa]

!... temperature, pressure, exner function, and geopotential height at layer interfaces:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    prsi,             &!pressure                                                               [Pa]
    prsik,            &!Exner function                                                          [-]
    phii,             &!geopotential height of layer                                            [m]
    tgrsi              !temperature                                                             [K]
 real(kind=RKIND),dimension(:,:),allocatable:: &
    prsi_hyd           !hydrostatic pressure                                                   [Pa]

!... layer-mean pressure and pressure at layer interfaces for radiation:
 real(kind=RKIND),dimension(:),allocatable:: &
    pgr_rad,           &!surface pressure                                                      [cb]
    sig_rad             !reference sigma levels                                                 [-]
 real(kind=RKIND),dimension(:,:),allocatable:: &
    prsl_rad,         &!mean-layer pressure                                                    [cb]
    prsi_rad           !pressure at layer interface                                            [cb]

!... tracers:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    qvrs               !water vapor       mixing ratio                                      [kg/kg]
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    qgrs,             &!tracer mixing ratios (including water vapor)                        [kg/kg]
    qtrcs              !tracer mixing ratios (excluding water vapor)                        [kg/kg]

!... updated arrays after call to physics:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    gt0,              &!updated temperature                                                     [K]
    gu0,              &!updated zonal wind                                                    [m/s]
    gv0                !updated meridional wind                                               [m/s]

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    gq0                !updated tracer array                                                    [-]

!...
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    dt3dt,            &!change in temperature due to physics processes                          [K]
    dq3dt,            &!change in tracers due to physics processes                          [kg/kg]
    du3dt,            &!change in u momentum due to physics processes                         [m/s]
    dv3dt              !change in v momentum due to physics processes                         [m/s]

!==================================================================================================
!OZONE INPUT DATA: This array should be moved to Registry.
!==================================================================================================

 real(kind=RKIND),dimension(:,:,:,:),allocatable:: &
    O3data             !input ozone data                                                        [?]

!==================================================================================================
!GFS VARIABLES AND ARRAYS:
!==================================================================================================

!... initial values for the surface roughness length as a function of the vegetation type:
 real(kind=RKIND),dimension(13):: z0_sib
 data z0_sib /2.653,0.826,0.563,1.089,0.854,0.856,0.035,0.238,0.065,0.076,0.011,0.035,0.011/

 real(kind=RKIND),dimension(:),allocatable:: &
    slmsk              !land/sea/sea-ice mask (sea=0;land=1;seaice=2)

 real(kind=RKIND):: &
    solhr,            &!
    clstp

!... miscellaneous:
 integer,dimension(:),allocatable:: &
    nlons              !nb of total grid-points
 real(kind=RKIND),dimension(:),allocatable:: &
    sinlat,           &!sin of latitude 
    coslat             !cos of latitude

!... temporary:
 real(kind=RKIND),dimension(:),allocatable:: &
    poz
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    prdout

!... arrays related to land surface schemes:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    slc,              &!liquid soil moisture
    smc,              &!total soil moisture
    stc                !soil temperature                                                        [K]

!... arrays related to gravity wave drag over orography:
 real(kind=RKIND),dimension(:),allocatable:: &
    hprime_rad         !input orography statistics for radiation
 real(kind=RKIND),dimension(:,:),allocatable:: &
    hprime             !input orography statistics. 

!... arrays related to shallow convection:
 real(kind=RKIND),dimension(:),allocatable:: &
    dpshc              !maximum pressure depth for shallow convection                         [hPa]

!... arrays related to deep convection:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    upd_mf,           &!convective updraft mass flux
    dwn_mf,           &!convective downdraft mass flux
    det_mf             !detrainment mass flux

!... arrays related to vertical diffusion:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    dkh                !vertical diffusion coefficient

!... arrays related to precipitation:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    rnp                !n cloud precipitation rate

!... array for cloud fraction:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    cldfrac            !cloud fraction                                                          [-]

!... arrays saved for restarts:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    phy_f2d
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    phy_f3d

!==================================================================================================
!... variables and arrays related to parameterization of longwave and shortwave radiation:
!==================================================================================================

 integer:: &
    idt                !

 integer:: &
    ipsd0              !initial value for random seed number generator.

 real(kind=RKIND),dimension(:),allocatable:: &
    sheleg

 real(kind=RKIND),dimension(:,:),allocatable:: &
    htrsw,            &!
    htrlw,            &!
    fluxr,            &!
    cldcov             !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    htrswb,           &!
    htrlwb             !

 real(kind=RKIND),dimension(:,:),allocatable:: &
    frain,            &!
    fcice,            &!
    rrime              !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    ozone

!==================================================================================================
 end module mpas_atmphys_gfs_vars
!==================================================================================================
