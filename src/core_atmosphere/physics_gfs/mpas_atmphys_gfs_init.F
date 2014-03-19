! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_init
 use mpas_configure
 use mpas_grid_types
 use mpas_timekeeping

 use mpas_atmphys_gfs_grrad_init
 use mpas_atmphys_gfs_o3climatology
 use mpas_atmphys_gfs_vars

!physics_gfs:
 use funcphys
 use mersenne_twister, only: random_number,random_setseed
 use tracer_const

 implicit none
 private
 public:: physics_init


!>\brief MPAS main initialization subroutine for the GFS physics suite.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-01-27.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Added the calculation of fhlwr and fhswr (frequencies of long- and short-wave radiation 
!>      calls per hour to subroutine GFS_initialize in mpas_atmphys_gfs_init.F, to compute the
!>      variables nslwr and nsswr (number of long- and short-wave radiation calls per hour). 
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-19. 


 contains

!==================================================================================================
 subroutine physics_init(dminfo,clock,config_do_restart,mesh,diag,state_1,state_2, &
                         diag_physics,atm_input,sfc_input)
!==================================================================================================

!input arguments:
 logical,intent(in):: config_do_restart
 type (dm_info), intent(in):: dminfo
 type(mesh_type),intent(in):: mesh
 type(diag_type),intent(in):: diag
 type(MPAS_Clock_type),intent(in):: clock

!inout arguments:
 type(state_type),intent(inout):: state_1,state_2
 type(diag_physics_type),intent(inout):: diag_physics
 type(atm_input_type),intent(inout):: atm_input
 type(sfc_input_type),intent(inout):: sfc_input

!local variables:
 type(MPAS_Time_Type):: initTime

 integer:: iCell,ierr
 integer:: hour,day,month,year

 integer:: seed0
 real(kind=kind_evod) wrk(1)

!--------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine physics_init:'

!initialization of east-north directions to convert u-tendencies from cell centers to cell
!edges:
 call init_dirs_forphys(mesh)

!retrieve the initial starting date and initialize the array idate: 
 call mpas_set_time(curr_time=initTime,dateTimeString=mesh%initial_date%scalar,ierr=ierr)
 call mpas_get_time(curr_time=initTime,YYYY=year,MM=month,DD=day,H=hour)
 idate(1) = hour
 idate(2) = month
 idate(3) = day
 idate(4) = year
 write(0,*) '... year  =', idate(4)
 write(0,*) '... month =', idate(2)
 write(0,*) '... day   =', idate(3)
 write(0,*) '... hour  =', idate(1)
 write(0,*) '--- end mpas_set_time:'

! currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
!call mpas_get_time(curr_time=currTime,DoY=julday,ierr=ierr)

!initialization of GFS non-radiation physics schemes:
 call GFS_initialize(dminfo,mesh)

!abort MPAS if the new SAS convection scheme is not used. If RAS or old SAS are used, we need
!to compute a random number seed (see additional sourcecode in gloopb.f):
 if(.not.newsas) then
    if(random_clds) then
       write(0,*) 'Error in choice of convection scheme, set newsas to true for now.'
       call mpas_dmpar_abort(dminfo)
    endif
 endif

!initialization of surface roughness. for now, we initialize and update the surface roughness
!length using the monthly-mean climatological surface roughness length zorl12m (see call to
!subroutine physics_gfs_update_surface):
 if(.not.config_do_restart) then
    do iCell = 1, mesh % nCellsSolve
       if(sfc_input%ivgtyp%array(iCell) >= 1 .and. sfc_input%ivgtyp%array(iCell) <= 13) then
          diag_physics%zorl%array(iCell) = z0_sib(sfc_input%ivgtyp%array(iCell))
       elseif(sfc_input%ivgtyp%array(iCell) == 0) then
          diag_physics%zorl%array(iCell) = 0.01
       endif
    enddo
 endif

!initialization of surface temperature over sea-ice cells:
 if(.not. config_do_restart) then
    do iCell = 1, mesh % nCellsSolve
       diag_physics % tisfc % array(iCell) = sfc_input % skintemp % array(iCell)
    enddo
 endif

!initialization of all physics function tables:
 call gfuncphys

!initialization of ozone climatology:
 call init_o3climatology(dminfo,mesh,atm_input)

!initialization of GFS radiation schemes:
 call radiation_init_forMPAS(dminfo,clock,mesh)

 write(0,*) '--- end subroutine physics_init.'
 write(0,*)

 end subroutine physics_init

!==================================================================================================
 subroutine GFS_initialize(dminfo,mesh)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh

!local variables:
!character(len=StrKIND):: gfs_namelist

 logical:: lgoc3d
 integer:: iret,me
 integer:: nmtvr_l

 real(kind=RKIND):: deltim
 real(kind=RKIND):: tfiltc,tol

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine GFS_initialize: dt_dyn=',dt_dyn
 write(0,*) '--- my_proc_id = ',dminfo % my_proc_id
 write(0,*) '--- IO_NODE    = ',IO_NODE

!if(dminfo % my_proc_id == IO_NODE) me = 0
 me = dminfo % my_proc_id

!... initialization of some parameters. some parameters are needed for the GFS physics; some
!other parameters are also defined, but not needed in MPAS:

 call compns(deltim,iret,ntrac,jcapg,jcap,levs,levr,lonf,lonr,latg,latr,ntoz,ntcw,ncld,lsoil,&
             nmtvr_l,num_p3d,num_p2d)

!... initialize logicals and arrays that are set through the namelist gfs_physics when running
!the GFS physics. Those are default values used when running the GFS operational model, using
!gfs_namelist.t574 from FangLin. Those options will be moved to namelist.input.

 deltim = config_dt        !model time-step.
 levs  = mesh % nVertLevels!number of model layers.
 lsoil = mesh % nSoilLevels!number of soil layers in the land-surface scheme.

 lssav            = .true. !logical that controls all outputs in subroutine gbphys.
 ldiag3d          = .true. !if set to true, output physics tendencies in subroutine gbphys.
 ccnorm           = .false.!if set to true, cloud condensates are normalized by cloud fraction.
 run_enthalpy     = .false.!not used.
 hybrid           = .true.
 out_virttemp     = .true.
 zflxtvd          = .true.
 gen_coord_hybrid = .false.

 mstrat           = .false.
 old_monin        = .false.
 mom4ice          = .false.
 cnvgwd           = .false.

 ntrac    = 3 !number of tracers (water vapor + ozone + cloud condensates).
 ncld     = 1 !number of cloud condensates.
 ntoz     = 2 !ozone location in the tracer array (should have ntoz = 2 instead).
 ntcw     = 3 !cloud condensate location in the tracer array.
 if(.not.allocated(fscav)) allocate(fscav(ntrac-ncld-1))

!... option for parameterization of deep and shallow convection:
 newsas           = .true. !if true, call subroutine sascnvn in subroutine gbphys.
 ras              = .false.!if true, call subroutine rascnv in subroutine gbphys.
 sashal           = .true. !
!... option for parameterization of cloud microphysics:
 zhao_mic         = .true. !if true, call the Zhao/Carr/Sundqvist microphysics scheme.
                           !if false, call the Ferrier cloud microphysics scheme.
!... option for the land surface scheme:
 lsm = 1                   !lsm = 0: use OSU land surface scheme.
                           !lsm = 1: use NOAH land surface scheme.
                            
!... initialization for radiation calculations:
 fhswr   = dt_radtsw / 3600._RKIND !frequency of shortwave radiation calls per hour.
 fhlwr   = dt_radtlw / 3600._RKIND !frequency of longwave radiation calls per hour.

 nfxr    = 33      !second dimension of input/output array fluxr.

 ialb    = 0               !control flag for the surface albedo.
                           !ialb  = 0: climatological surface albedo, based on surface types.
                           !ialb  = 1: read from an input file.
 iaer    = 111             !control flag for aerosol effects.
                           !iaer   = 0  : turn all aerosol effects of (sw,lw,volc)
                           !iaer   = 111: volc and clim trops aerosol for both sw and lw.
                           !iaer   = .  : see subroutine aerinit for description of all other options.
 ico2    = 0               !control flag for input CO2 data.
                           !ico2   = 0  : use prescribed global mean CO2.
                           !ico2   = 1  : use observed CO2 global annual mean value.
                           !ico2   = 2  : use observed CO2 monthly data with 2D variation.
 ictm    = 1               !external data time/date control flag.
                           !ictm   = 0  : use data at initial time.
                           !ictm   = .  : see all *init subroutines for description of all other options.
 iems    = 0               !control flag for the surface emissivity.
                           !iems   = 0  : the surface emissivity is set to 1.
                           !iems   = 1  : climatological surface emissivity, based on surface type.
 iovr_lw  = 1              !cloud overlapping control flag in parameterization of longwave radiation.
                           !iovr_lw = 0: random cloud overlap.
                           !iovr_lw = 1: maximum/random cloud overlap.
 iovr_sw  = 1              !cloud overlapping control flag.
                           !iovr_sw = 0: random cloud overlap.
                           !iovr_sw = 1: maximum/random cloud overlap.
 isol    = 0               !control flag for solar constants:
                           !isol   = 0  : use a fixed solar constant.
                           !isol   = 1  : use 11-year cycle constant table.
 isubc_lw = 2              !lw clouds without sub-grid approximation.
                           !isubc_lw = 1: sub-grid clouds with prescribed seeds.
                           !isubc_lw = 2: sub-grid clouds with randomly generated seeds.
 isubc_sw = 2              !sw clouds without sub-grid approximation.
                           !isubc_sw = 1: sub-grid clouds with prescribed seeds.
                           !isubc_sw = 2: sub-grid clouds with randomly generated seeds.
!... other dimensions needed in the long- and short-wave radiation codes:
 lonf    = mesh%nCellsSolve
 lonr    = mesh%nCellsSolve
 if(ngptc.gt.lonr) ngptc = lonr
 latg   = 1
 latr   = 1
 jcap   = 1

!logical that controls outputs:

!range of droplet number concentration for Ferrier microphysics:
 ncw(1)   =  50.
 ncw(2)   = 150.
!range of minimum large ice fraction for Ferrier microphysics:
 flgmin   = 0.220 
!critical relative humidity at the surface,PBL top,and top-of-the-atmosphere:
 crtrh(1) = 0.85
 crtrh(2) = 0.85
 crtrh(3) = 0.85
!critical cloud-top entrainment instability:
 ctei_rm  = 0.50

!number of statistical fields for the parameterization of gravity wave drag:
 nmtvr    = mesh % nMtvr
!random number:
 nrcm     = 2 

!originally defined in nam_mrf:
 lgoc3d   = .false.
 tfiltc   = 0.85
 tol      = 0.01

!... this part is identical to the bottom of subroutine ./physics_gfs/compns, i.e. after
!reading gfs_namelist, except for parts that were commented out:
 filta = tfiltc
 phigs = phigs_d * (atan(1.0)/45.0)
 if(me == 0 .and. semilag) then
    write(0,*)' phigs_d=',phigs_d,' phigs=',phigs
 endif

 if(ndsl) then
    if (jcapg < 62 .or. jcapg > jcap) jcapg = jcap
 endif

 lggfs3d = lgoc3d

 if (levs > 100 .and. ref_temp < 400.0) ref_temp = 1500.0

 if(me == 0) then
    write(0,*) ' The time filter coefficient tfiltc=',tfiltc
    if(adiab) then
       write(0,*)' This is an adiabatic run'
    else
       if(lsm == 1) then
          write(0,*)' NOAH Land Surface Model used'
       elseif (lsm == 0) then
          write(0,*)' OSU Land Surface Model used'
       else
          write(0,*)' Unsupported LSM type - job aborted - lsm=',lsm
          call mpi_quit(2222)
       endif

       if(ras) then
          write(0,*) ' RAS Convection scheme used with ccwf=',ccwf
       else
          if(newsas) then
             write(0,*) ' New modified SAS Convection scheme used'
          else
             write(0,*) ' OPR SAS Convection scheme used'
          endif
       endif
       if(.not. old_monin) write(0,*)' New PBL scheme MONINQ used'
       if(.not. shal_cnv) then
          write(0,*) ' No shallow convection used'
       else
          if(sashal) write(0,*) ' New Massflux based shallow convection used'
       endif
       if(cnvgwd) write(0,*) ' Convective GWD parameterization used'
       if(crick_proof) write(0,*) ' CRICK-Proof cloud water used in radiation '
       if (ccnorm) write(0,*) ' Cloud condensate normalized by cloud cover for radiation'
    endif
 endif

 if(levr == 0) then
    levr = levs
 endif

 if(me ==  0) then
    if(.not. adiab) then
       write(0,*) ' Radiative heating calculated at',levr, ' layers'
       if(iovr_sw == 0) then
          write(0,*) ' random cloud overlap for Shortwave IOVR_SW=', iovr_sw
       else
          write(0,*) ' max-random cloud overlap for Shortwave IOVR_SW=',iovr_sw
       endif
       if(iovr_lw == 0) then
          write(0,*) ' random cloud overlap for Longwave IOVR_LW=',iovr_lw
       else
          write(0,*) ' max-random cloud overlap for Longwave IOVR_LW=',iovr_lw
       endif
       if(isubc_sw == 0) then
          write(0,*) ' no sub-grid cloud for Shortwave ISUBC_SW=',isubc_sw
       else
          write(0,*) ' sub-grid cloud for Shortwave ISUBC_SW=',isubc_sw
       endif
       if(isubc_lw == 0) then
          write(0,*) ' no sub-grid cloud for Longwave ISUBC_LW=',isubc_lw
       else
          write(0,*) ' sub-grid cloud for Longwave ISUBC_LW=',isubc_lw
       endif
    endif
 endif

 if(zhao_mic) then  ! default setup for Zhao Microphysics
    num_p3d = 4
    num_p2d = 3
    if(me == 0 .and. .not. adiab) &
       write(0,*) ' Using Zhao Microphysics : nump3d=',num_p3d, &
                  ' num_p2d=',num_p2d,' crtrh=',crtrh
 else               ! Brad Ferrier's Microphysics
    num_p3d = 3
    num_p2d = 1
    if(me == 0 .and. .not. adiab) &
       write(0,*) ' Using Ferrier Microphysics : nump3d=',num_p3d, &
                  ' num_p2d=',num_p2d,' crtrh=',crtrh,' ncw=',ncw,' flgmin=',flgmin
 endif

!
!sela - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tol=0.01
!  Check rule 1.
!     if(deltim.le.0) then
!       iret=1
!       return
!     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsout and check rule 2.
!     if(nsout.gt.0) fhout=nsout*deltim/3600.
!     nsout=nint(fhout*3600./deltim)
!     if(nsout.le.0.or.abs(nsout-fhout*3600./deltim).gt.tol) then
!       iret=2
!       return
!     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsout_hf and check rule 21.
!     if(nsout_hf.gt.0) fhout=nsout_hf*deltim/3600.
!     nsout_hf=nint(fhout_hf*3600./deltim)
!     if(nsout_hf <= 0.or.abs(nsout_hf-fhout_hf*3600./deltim)>tol) then
!       iret=9
!       return
!     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsswr and check rule 3.
      nsswr=nint(fhswr*3600./deltim)
      if(nsswr.le.0.or.abs(nsswr-fhswr*3600./deltim).gt.tol) then
        iret=3
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nslwr and check rule 4.
      nslwr=nint(fhlwr*3600./deltim)
      if(nslwr.le.0.or.abs(nslwr-fhlwr*3600./deltim).gt.tol.or. &
         mod(nslwr,nsswr).ne.0) then
        iret=4
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nszer and check rule 5.
!     nszer=nint(fhzer*3600./deltim)
!     if(nszer.le.0.or.abs(nszer-fhzer*3600./deltim).gt.tol.or.
!    &   mod(nszer,nsout).ne.0) then
!       iret=5
!       return
!     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsres and check rule 6.
!     nsres=nint(fhres*3600./deltim)
!     if(nsres.le.0.or.abs(nsres-fhres*3600./deltim).gt.tol.or.
!    &   mod(nsres,nslwr).ne.0.or.mod(nsres,nszer).ne.0) then
!       iret=6
!       return
!     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsdfi and check rule 7.
!     if(fhdfi.eq.0.) then
!       nsdfi=0
!     else
!       nsdfi=nint(fhdfi*3600./deltim)
!       if(nsdfi.le.0.or.abs(nsdfi-fhdfi*3600./deltim).gt.tol.or.
!    &     mod(nsdfi,nslwr).ne.0.or.nsdfi.gt.nsres) then
!         iret=7
!         return
!       endif
!     endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nscyc and check rule 8.
!     if(fhcyc.eq.0.) then
!       nscyc=0
!     else
!       nscyc=nint(fhcyc*3600./deltim)
!       if(nscyc.le.0.or.abs(nscyc-fhcyc*3600./deltim).gt.tol.or.
!    &     mod(nscyc,nslwr).ne.0) then
!         iret=8
!         return
!       endif
!     endif

      IF (NGPTC.GT.LONR) THEN
         NGPTC=LONR
         WRITE(*,*) "NGPTC IS TOO BIG, RESET NGPTC TO LONR",NGPTC
      ENDIF
      IF (ME.EQ.0)   WRITE(*,*) "NGPTC IS SET TO NGPTC :",NGPTC
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  All checks are successful.
!     iret=0

!     return
!
! 999 print *,' error reading  namelist - execution terminated by user'
!      if (me.eq.0) then
!       write(6,nam_mrf)
!       print*,' YHALO =', YHALO
!     endif
!     call mpi_quit(999)
!

!... read soil and vegetation data. Logical LPARAM was initially defined in the GFS namelist:
 lparam = .false.
 call set_soilveg(me)

!... check the number of tracers initialized at the top of GFS_initialize:
 call set_tracer_const(me,ntrac)

!... initialization of the 2nd dimension for the random number rann needed in the "RAS"
!    parameterization of deep convection:
 if(ras) then
    if(fix_ncld_hr) then
       nrcm = min(nrcmax, levs-1) * (deltim/1200) + 0.10001
    else
       nrcm = min(nrcmax, levs-1)
    endif
 else
    nrcm = 1
 endif

 write(0,*) '--- end subroutine GFS_initialize'
 write(0,*)

 end subroutine GFS_initialize

!==================================================================================================
 subroutine init_dirs_forphys(mesh)
!==================================================================================================

!inout arguments:
!----------------
 type(mesh_type),intent(in):: mesh

!local variables:
 integer:: iCell
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!--------------------------------------------------------------------------------------------------

 latCell => mesh % latCell % array
 lonCell => mesh % lonCell % array
 east    => mesh % east  % array
 north   => mesh % north % array

!Compute unit vectors in east and north directions for each cell:
 do iCell = 1, mesh % nCells

    east(1,iCell) = -sin(lonCell(iCell))
    east(2,iCell) =  cos(lonCell(iCell))
    east(3,iCell) =  0.0
    call r3_normalize(east(1,iCell), east(2,iCell), east(3,iCell))

    north(1,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
    north(2,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
    north(3,iCell) =  cos(latCell(iCell))
    call r3_normalize(north(1,iCell), north(2,iCell), north(3,iCell))

 end do

 end subroutine init_dirs_forphys

!==================================================================================================
 subroutine r3_normalize(ax, ay, az)
!==================================================================================================
!normalizes the vector (ax, ay, az)

 real (kind=RKIND), intent(inout) :: ax, ay, az
 real (kind=RKIND) :: mi

!--------------------------------------------------------------------------------------------------

 mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

 end subroutine r3_normalize

!==================================================================================================
 end module mpas_atmphys_init
!==================================================================================================


