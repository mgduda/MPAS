! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_model_parameters

 use machine
 use ozne_def

 implicit none
 public
 save
 

!=================================================================================================================
!... THESE PARAMETERS SHOULD BE MOVED TO REGISTRY.PHYSICS_NGGPS.NAMELISTS.XML:
!=================================================================================================================

 logical,parameter:: ras              = .false. !flag for ras convection parameterization.
 logical,parameter:: newsas           = .true.  !flag for calling new sas convection scheme.
 logical,parameter:: cal_pre          = .true.  !flag for controlling precipitation type algorithm.
 logical,parameter:: ccnorm           = .false. !
 logical,parameter:: cnvgwd           = .true.  !flag for convective wave drag.
 logical,parameter:: crick_proof      = .false. !
 logical,parameter:: dspheat          = .false. !flag for TKE dissipative heating. 
 logical,parameter:: gen_coord_hybrid = .true.  !flag for Henry's generalized coordinate.
 logical,parameter:: hybedmf          = .false. !flag for hybrid edmf PBL scheme.
 logical,parameter:: ldiag3d          = .false. !flag for computing diagnostics.
 logical,parameter:: lgocart          = .false. !flag for 3d diagnostic fields for gocart 1
 logical,parameter:: lsidea           = .false. !not sure what it does yet.
 logical,parameter:: lssav_cpl        = .false. !flag for save data for A/O/I coupling.
 logical,parameter:: moist_adj        = .false. !flag for moist convective adjustment scheme.
 logical,parameter:: mom4ice          = .false. !flag controls mom4 sea ice option.
 logical,parameter:: norad_precip     = .false. !
 logical,parameter:: pdfcld           = .false. !flag for stratiform clouds (does not seem to be used).
 logical,parameter:: pre_rad          = .false. !flag for testing purposes.
 logical,parameter:: redrag           = .true.  !flag for reduced drag coefficient over oceans.
 logical,parameter:: shal_cnv         = .false. !flag for calling shallow convection scheme.
 logical,parameter:: shcnvcw          = .false. !flag for saving shallow convective cloud diagnostics.
 logical,parameter:: trans_trac       = .true.  !flag for convective transport of tracers.
 logical,parameter:: zhao_mic         = .true.  !flag for cloud microphysics.

!... new logicals added in nuopc_physics.F90 version 2:
 logical,parameter:: shoc_cld         = .false. !flag for SHOC parameterization of cloud fraction.
 logical,parameter:: do_shoc          = .false. !flag for SHOC parameterization of cloud fraction.
 logical,parameter:: shocaftcnv       = .false. !flag for SHOC parameterization of cloud fraction.

 logical,parameter:: cscnv            = .false. !flag for Chikira-Sugiyama convection scheme.

 integer,parameter:: ntrac            = 3       !number of tracers (water vapor+ozone+cloud condensates).
 integer,parameter:: ncld             = 1       !number of cloud condensates.
 integer,parameter:: ntoz             = 2       !ozone index location in the tracer array.
 integer,parameter:: ntcw             = 3       !cloud condensate index location in the tracer array.
 integer,parameter:: ntke             = 0       ! ?
 integer,parameter:: thermodyn_id     = 1       !    
 integer,parameter:: sfcpress_id      = 2       !

 integer,parameter:: iflip = 1                  !control flag for direction of vertical index.
                                                !iflip = 0: index is from top-of-the-atmosphere to surface.
                                                !iflip = 1: index is from surface to top-of-the-atmosphere.
 integer,parameter:: lsm = 1                    !lsm = 0: use OSU land surface scheme.
                                                !lsm = 1: use NOAH land surface scheme.

!=================================================================================================================
!... OTHER PARAMETERS THAT ARE DEFINED DURING AT INITIALIZATION TIME-STEP:
!=================================================================================================================

 logical:: &
    flipv,                &!
    old_monin,            &!flag for old surface monin scheme.
    mstrat,               &!flag for moorthi approach for stratus cloud scheme.
    sashal,               &!flag for calling new shallow convection scheme.
    sas_shal               !

 integer,parameter:: nst_fcst = 0

 integer,parameter:: latr   = 1  !number of latitude points.
 integer,parameter:: jcap   = 1  !number of spectral wave truncation for sascnv and shalcnv.
 integer,parameter:: nlunit = 35 !file unit.

 integer,dimension(2),parameter:: ncw = (/20,120/) !for Ferrier microphysics.

 integer:: &
    levr,                 &!
    levs,                 &!
    lonr,                 &!
    lsoil,                &!
    me,                   &!
    npdf3d,               &!
    nmtvr,                &!
    nrcm,                 &!second dimension for the random number array rann.
    ntot2d,               &!
    ntot3d,               &!
    num_p3d,              &!
    num_p2d                !
 integer,dimension(4):: &
    idate                  !

 real(kind=kind_phys),parameter:: prslrd0 = 200.
 real(kind=kind_phys),parameter,dimension(2):: cgwf    = (/0.5,0.05/)
 real(kind=kind_phys),parameter,dimension(2):: cdmbgwd = (/2.0,0.25/)       !for Zhao microphysics.
 real(kind=kind_phys),parameter,dimension(3):: crtrh   = (/0.90,0.90,0.90/) !for Zhao microphysics.

 real(kind=kind_phys):: dxmaxin, dxminin, dxinvin
 real(kind=kind_phys),dimension(2):: ccwf
 real(kind=kind_phys),dimension(2):: ctei_rm
 real(kind=kind_phys),dimension(2):: dlqf

 real(kind=kind_phys),dimension(:),allocatable:: si

 logical,parameter:: nocnv = .false.

 end module mpas_atmphys_nggps_model_parameters
!=================================================================================================================
