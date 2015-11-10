! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_vars
 use mpas_kind_types

 use nuopc_physics
 use machine
 use funcphys
 use mersenne_twister
 use module_ras
 use physcons
 use physparam
 use ozne_def
 use cs_conv,only: nctp

 implicit none
 public
 save

!=================================================================================================================
!... logicals defined in MPAS:
 logical:: l_radtlw                   !controls call to longwave radiation parameterization.
 logical:: l_radtsw                   !controls call to shortwave radiation parameterization.

 real(kind=RKIND),public:: dt_dyn     !time-step for dynamics
 real(kind=RKIND),public:: dt_radtlw  !time-step for longwave radiation parameterization  [mns]
 real(kind=RKIND),public:: dt_radtsw  !time-step for shortwave radiation parameterization [mns]

!... initial values for the surface roughness length as a function of the vegetation type:
 real(kind=RKIND),dimension(13):: z0_sib
 data z0_sib /2.653,0.826,0.563,1.089,0.854,0.856,0.035,0.238,0.065,0.076,0.011,0.035,0.011/

!=================================================================================================================
!... TYPES USED IN THE NUOPC DRIVER AND EXTENDED TYPES USED IN MPAS NGGPS PHYSICS:
!=================================================================================================================

 type(model_parameters)  :: mdl_parameters
 type(dynamic_parameters):: dyn_parameters


!=================================================================================================================
!... ADDITIONAL INPUTS NEEDED FOR INITIALIZATION OF THE NGGPS DRIVER
!... (SHOULD BE MOVED TO REGISTRY.PHYSICS_NGGPS.NAMELISTS.XML):
!=================================================================================================================

 integer,dimension(8):: idat                            !initial date but different size and ordering (radupdate).

!... long- and short-wave radiation only:
 integer,parameter:: &
 ialb     = 0              !control flag for the surface albedo.
                           !ialb  = 0: climatological surface albedo, based on surface types.
                           !ialb  = 1: read from an input file.
 integer,parameter:: &
 iaer     = 111            !control flag for aerosol effects.
                           !iaer   = 0  : turn all aerosol effects of (sw,lw,volc)
                           !iaer   = 111: volc and clim trops aerosol for both sw and lw.
                           !iaer   = .  : see subroutine aerinit for description of all other options.
 integer,parameter:: &
 ico2     = 2              !control flag for input CO2 data.
                           !ico2   = 0  : use prescribed global mean CO2.
                           !ico2   = 1  : use observed CO2 global annual mean value.
                           !ico2   = 2  : use observed CO2 monthly data with 2D variation.
 integer,parameter:: &
 ictm     = 1              !external data time/date control flag.
                           !ictm   = 0  : use data at initial time.
                           !ictm   = *  : see all *init subroutines for description of all other options.
 integer,parameter:: &
 iems     = 1              !control flag for the surface emissivity.
                           !iems   = 0  : the surface emissivity is set to 1.
                           !iems   = 1  : climatological surface emissivity, based on surface type.
 integer,parameter:: &
 iovr_lw  = 1              !cloud overlapping control flag in parameterization of longwave radiation.
                           !iovr_lw = 0: random cloud overlap.
                           !iovr_lw = 1: maximum/random cloud overlap.
 integer,parameter:: &
 iovr_sw  = 1              !cloud overlapping control flag.
                           !iovr_sw = 0: random cloud overlap.
                           !iovr_sw = 1: maximum/random cloud overlap.
 integer,parameter:: &
 isol     = 2              !control flag for solar constants:
                           !isol   = 0  : use a fixed solar constant.
                           !isol   = 1  : use 11-year cycle constant table.
                           !isol   = 1  :
 integer,parameter:: &
 isubc_lw = 2              !lw clouds without sub-grid approximation.
                           !isubc_lw = 1: sub-grid clouds with prescribed seeds.
                           !isubc_lw = 2: sub-grid clouds with randomly generated seeds.
 integer,parameter:: &
 isubc_sw = 2              !sw clouds without sub-grid approximation.
                           !isubc_sw = 1: sub-grid clouds with prescribed seeds.
                           !isubc_sw = 2: sub-grid clouds with randomly generated seeds.


!=================================================================================================================
!... LOCAL VARIABLES AND ARRAYS DEFINED IN THE NGGPS TYPE "DYNAMIC_PARAMETERS":
!=================================================================================================================

!... FOR RADIATION ONLY:
 integer,dimension(8):: jdate                           !current forecast time (year,month,day,t-zone,hour,minute,
                                                        !second,milli-second).
!=================================================================================================================
!... LOCAL VARIABLES AND ARRAYS DEFINED IN THE NGGPS TYPE "STATE_FIELDS_IN":
!=================================================================================================================

!... FOR PHYSICS AND RADIATION:
 real(kind=kind_phys),dimension(:,:),allocatable:: &
    prsi,             &!pressure at layer interface                                                          [Pa].
    prsl,             &!pressure                                                                             [Pa].
    prslk,            &!Exner function                                                                        [-].
    tgrs               !temperature                                                                           [K].

!... FOR RADIATION ONLY:
 real(kind=kind_phys),dimension(:,:),allocatable:: &
    qgrs_rad           !water vapor mixing ratio for radiation                                            [kg/kg].

 real(kind=kind_phys),dimension(:,:,:),allocatable:: &
    tracer

 real(kind=kind_phys),dimension(:,:),allocatable:: &
    vvl                !vertical velocity                                                                   [m/s].

!... FOR PHYSICS ONLY:
 real(kind=kind_phys),dimension(:),allocatable:: &
    pgr,              &!surface pressure                                                                     [Pa].
    adjtrc             !

 real(kind=kind_phys),dimension(:,:),allocatable:: &
    ugrs,             &!u-velocity interpolated to theta points                                             [m/s].
    vgrs,             &!v-velocity interpolated to theta points                                             [m/s].
    prsik,            &!
    phii,             &!
    phil               !

 real(kind=kind_phys),dimension(:,:,:),allocatable:: &
    qgrs               !

!... ADDITIONAL MPAS ARRAYS NEEDED FOR THE INITIALIZATION OF "STATE_FIELDS_IN" ARRAYS:
 real(kind=RKIND),dimension(:),allocatable:: &
    pgr_hyd            !

 real(kind=kind_phys),dimension(:,:),allocatable:: &
    rhos,             &!
    dzgrid,           &!
    fzm,              &!
    fzp,              &!
    prsl_hyd,         &!    
    tgrsi,            &!
    prsi_hyd           !


!=================================================================================================================
!... LOCAL VARIABLES AND ARRAYS DEFINED IN THE NGGPS TYPE "STATE_FIELDS_OUT":
!=================================================================================================================

!... FOR PHYSICS ONLY:
 real(kind=kind_phys),dimension(:,:),allocatable:: &
    gt0,             &!
    gu0,             &!
    gv0               !

 real(kind=kind_phys),dimension(:,:,:),allocatable:: &
    gq0


!=================================================================================================================
!... LOCAL VARIABLES AND ARRAYS DEFINED IN THE NGGPS TYPE "DIAGNOSTICS":
!=================================================================================================================

!... input/output - used in radiation:
 real(kind=kind_phys),dimension(:,:),allocatable:: &
    fluxr,          &!
    dswcmp,         &!
    uswcmp           !

!... input/output - used in physics:
 real(kind=kind_phys),dimension(:,:),allocatable:: &
    dqdt_v           !
 real(kind=kind_phys),dimension(:,:,:),allocatable:: &
    dt3dt,          &!
    dq3dt,          &!
    du3dt,          &!
    dv3dt            !

!=================================================================================================================
!... LOCAL VARIABLES AND ARRAYS DEFINED IN THE NGGPS TYPE "CLOUD_PROPERTIES":
!=================================================================================================================

 real(kind=kind_phys),dimension(2),parameter:: flgmin  = (/0.180,0.220/)
 real(kind=kind_phys),parameter:: sup = 1.1                      !supersaturation in pdf cloud when t is very low

!=================================================================================================================
!... LOCAL VARIABLES AND ARRAYS DEFINED IN THE NGGPS TYPE "TBD_DDT":
!=================================================================================================================

!... INPUTS:
 real(kind=kind_phys):: xkzm_m
 real(kind=kind_phys):: xkzm_h
 real(kind=kind_phys):: xkzm_s

 real(kind=kind_phys):: evpco
 real(kind=kind_phys),dimension(2):: psautco
 real(kind=kind_phys),dimension(2):: prautco
 real(kind=kind_phys),dimension(2):: wminco

!real(kind=kind_phys),dimension(:),allocatable:: &
!   dpshc,           &!
!   poz               !
!real(kind=kind_phys),dimension(:,:),allocatable:: &
!   rann              !
!real(kind=kind_phys),dimension(:,:,:),allocatable:: &
!   prdout            !

!real(kind=kind_phys),dimension(:,:),allocatable:: &
!   phy_f2d           !
!real(kind=kind_phys),dimension(:,:,:),allocatable:: &
!   phy_f3d           !

!=================================================================================================================
 end module mpas_atmphys_nggps_vars
!=================================================================================================================
