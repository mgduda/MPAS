!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_init
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_timekeeping

 use mpas_atmphys_nggps_model_parameters
 use mpas_atmphys_nggps_o3climatology
 use mpas_atmphys_nggps_vars
 use nuopc_physics

 implicit none
 private
 public:: physics_init


!modifications to sourcecode:
!----------------------------
! * in the subroutine physics_init, added the initialization of the variables fhour, nlonsCell, coslatCell, and
!   sinlatCell. removed write statements after call to the subroutine nuopc_phys_init.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-28.
! * in the subroutine physics_init, added deallocation of array si after call to subroutine nuopc_phys_init.
!   Laura D. Fowler (laura@ucar.edu) / 2015-12-01.
! * in the subroutine physics_init, added the calculation of the total number of cells (nCellsTotal).
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-10.
! * in the subroutine physics_init, changed the definition of oro.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-21.
! * in the subroutine physics_init, commented out the initialization of zorl since it is initialized when we
!   build the initial conditions for the NGGPS project. we may have to correct this if zorl is not available
!   in the intermediate files. also, multiply zorl by 0.01 for conversion to cm.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-25.
! * in the subroutine physics_init, modified the initialization of the "tank" of random numbers rannum_tank so
!   that random numbers can be updated in the physics driver.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-29.

 contains


!=================================================================================================================
 subroutine physics_init(dminfo,clock,configs,mesh,nggps_atm_input,nggps_dynamic_parameters,nggps_sfc_properties,&
                         nggps_tbd_ddt)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(MPAS_Clock_type),intent(in):: clock

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_atm_input
 type(mpas_pool_type),intent(inout):: nggps_sfc_properties
 type(mpas_pool_type),intent(inout):: nggps_dynamic_parameters
 type(mpas_pool_type),intent(inout):: nggps_tbd_ddt

!local pointers:
 logical,pointer:: config_do_restart

 character(len=StrKIND),pointer:: initial_date

 integer,pointer:: nCellsSolve
 integer,pointer:: nVertLevels,nSoilLevels
 integer,pointer:: nOznLevels,nOznPlCoeff
 integer,pointer:: nMtgwd
 integer,pointer:: numP2d,numP3d
 integer,pointer:: maxran,maxsub
 integer,dimension(:),pointer:: isltyp,ivgtyp
 integer,dimension(:),pointer:: nlonsCell

 real(kind=PHYSKIND),pointer:: fhour
 real(kind=RKIND),dimension(:),pointer:: rdzw
 real(kind=PHYSKIND),dimension(:),pointer:: oro,oro_uf
 real(kind=PHYSKIND),dimension(:),pointer:: tisfc,tsea
 real(kind=PHYSKIND),dimension(:),pointer:: zorl
 real(kind=PHYSKIND),dimension(:),pointer:: latCell,coslatCell,sinlatCell

 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 real(kind=PHYSKIND),dimension(:,:),pointer:: rannum_tank

!local variables:
 type(MPAS_Time_Type)  :: initTime

 logical:: hybrid

 integer:: nCellsTotal
 integer:: iCell,ierr,n
 integer:: hour,day,month,year
 integer:: nshoc_2d,nshoc_3d

!local variables for calculation of random numbers:
 logical:: random_clds

 integer:: iiCell,nn,nrc,k,kk
 integer:: lonrbm,iseedl,maxrs
 integer,dimension(:), allocatable:: indxr

 real(kind=PHYSKIND),parameter:: thousnd = 1000._PHYSKIND
 real(kind=PHYSKIND):: wrk(1)
 real(kind=PHYSKIND),dimension(:), allocatable   :: rannum

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine physics_init NGGPS:'

 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)
 call mpas_pool_get_dimension(mesh,'nMtvr'      ,nMtgwd     )
 call mpas_pool_get_dimension(mesh,'numP2d'     ,numP2d     )
 call mpas_pool_get_dimension(mesh,'numP3d'     ,numP3d     )

 call mpas_pool_get_array(mesh,'initial_date',initial_date)
 call mpas_pool_get_array(mesh,'rdzw',rdzw)
 write(0,*) '--- initial_date = ', trim(initial_date)

 call mpas_pool_get_array(mesh,'nlonsCell'  ,nlonsCell )
 call mpas_pool_get_array(mesh,'latCellPhys',latCell   )
 call mpas_pool_get_array(mesh,'coslatCell' ,coslatCell)
 call mpas_pool_get_array(mesh,'sinlatCell' ,sinlatCell)
 call mpas_pool_get_array(mesh,'zgrid'      ,zgrid     )

 call mpas_pool_get_array(nggps_dynamic_parameters,'forecast_hour',fhour)

 call mpas_pool_get_array(nggps_sfc_properties,'isltyp',isltyp)
 call mpas_pool_get_array(nggps_sfc_properties,'ivgtyp',ivgtyp)
 call mpas_pool_get_array(nggps_sfc_properties,'oro'   ,oro   )
 call mpas_pool_get_array(nggps_sfc_properties,'oro_uf',oro_uf)
 call mpas_pool_get_array(nggps_sfc_properties,'tisfc' ,tisfc )
 call mpas_pool_get_array(nggps_sfc_properties,'tsea'  ,tsea  )
 call mpas_pool_get_array(nggps_sfc_properties,'zorl'  ,zorl  )

!--- initialization of east-north directions to convert u-tendencies from cell centers to cell edges:
 call init_dirs_forphys(mesh)

!--- retrieve the initial starting date and initialize the array idate: 
 call mpas_set_time(curr_time=initTime,dateTimeString=initial_date,ierr=ierr)
 call mpas_get_time(curr_time=initTime,YYYY=year,MM=month,DD=day,H=hour)
 idate(1) = hour
 idate(2) = month
 idate(3) = day
 idate(4) = year
 write(0,*) '... year  =', idate(4)
 write(0,*) '... month =', idate(2)
 write(0,*) '... day   =', idate(3)
 write(0,*) '... hour  =', idate(1)

!--- initialization of dimensions needed in the NGGPS physics:
 lonr     = nCellsSolve
 levr     = nVertLevels
 levs     = nVertLevels
 lsoil    = nSoilLevels
 levozp   = nOznLevels
 pl_coeff = nOznPlCoeff
 nmtvr    = nMtgwd

 me = dminfo % my_proc_id

!--- initialization of the sigma levels needed in subroutine cld_init to define the top layer with PBL clouds:
!--- note that this will need to be changed as MPAS uses a height-coordinate system.
 if(.not.allocated(si)) allocate(si(levr+1))
 si(levr+1) = 0._kind_phys
 do k = 1, levr
    si(k) = float(levr+1-k) / float(levr)
 enddo

!--- initialization of the local variables needed to run the SHOC parameterization. Although we do not run SHOC,
!--- nshoc_2d and nshoc_3d are initialized locally since they are needed to initialize ntot2d and ntot3d in the
!--- subroutine nuopc_phy_init. for now, we set nshoc_3d and nshoc_2d equal to zero:
 nshoc_3d = 0
 nshoc_2d = 0

!--- check the values for num_p3d and num_p2d depending on the choice of the cloud microphysics scheme
!--- (zhao-carr or ferrier). check that num_p3d and num_p2d are as numP3d and numP2d, otherwise, set
!--- numP3d and numP2d to num_p3d and num_p2d and recompile MPAS:
 npdf3d = 0
 if(zhao_mic) then
    num_p3d = 4
    num_p2d = 3
    if(pdfcld) then
       npdf3d = 3
    endif
 else
    num_p3d = 3
    num_p2d = 1
 endif

 if(numP3d /= num_p3d .or. numP2d /= num_p2d) then
    write(0,*) '... numP3d  =', numP3d
    write(0,*) '... numP2d  =', numP2d
    write(0,*) '... num_p3d =', num_p3d
    write(0,*) '... num_p2d =', num_p2d
    write(0,*) '... reset numP3d (numP2d) to num_p3d (num_p2d), and recompile MPAS:'
    call mpas_dmpar_abort(dminfo)
 endif

!--- initialization of ntot2d and ntot3d:
 ntot3d = num_p3d + nshoc_3d + npdf3d
 ntot2d = num_p2d + nshoc_2d

!--- initialization of some logicals and parameters depending on the choice of convection scheme
!--- (ras or sas) for inputs in the type "model_parameters" in nuopc_physics.F:
 old_monin = .false.
 sashal    = .false.
 mstrat    = .false.
 if(ras) then
    old_monin  = .true.
    mstrat     = .true.
    sashal     = .false.
    ccwf(1)    = 0.0_kind_phys
    ccwf(2)    = 1.0_kind_phys
    dlqf(1)    = 0.5_kind_phys
    dlqf(2)    = 0.5_kind_phys
    ctei_rm(1) = 0.60_kind_phys
    ctei_rm(2) = 0.23_kind_phys
 elseif(newsas) then
    old_monin  = .false.
    mstrat     = .false.
    sashal     = .true.
    ccwf(1)    = 1.0_kind_phys
    ccwf(2)    = 1.0_kind_phys
    dlqf(1)    = 0.0_kind_phys
    dlqf(2)    = 0.0_kind_phys
    ctei_rm(1) = 10.0_kind_phys
    ctei_rm(2) = 10.0_kind_phys
 endif

!--- initialization of some logicals and parameters depending on the choice of convection scheme
!--- (ras or sas) for inputs in the type "tbd_ddt" in nuopc_physics.F:
 if(ras) then
    xkzm_m     = 0.25_kind_phys
    xkzm_h     = 0.25_kind_phys 
    xkzm_s     = 1.0_kind_phys
    psautco(1) = 6.0e-4_kind_phys
    psautco(2) = 3.0e-4_kind_phys
    prautco(1) = 6.0e-4_kind_phys
    prautco(2) = 3.0e-4_kind_phys
    wminco(1)  = 5.0e-5_kind_phys
    wminco(2)  = 2.5e-5_kind_phys
    evpco      = 2.0e-5_kind_phys
 elseif(newsas) then
    xkzm_m     = 1.0_kind_phys
    xkzm_h     = 1.0_kind_phys 
    xkzm_s     = 1.0_kind_phys
    psautco(1) = 6.0e-4_kind_phys
    psautco(2) = 3.0e-4_kind_phys
    prautco(1) = 1.0e-4_kind_phys
    prautco(2) = 1.0e-4_kind_phys
    wminco(1)  = 1.0e-5_kind_phys
    wminco(2)  = 1.0e-5_kind_phys
    evpco      = 2.0e-5_kind_phys
 endif

!--- initialization of random number needed for ras or old version of sas:
 random_clds = .false.
 if(ras .or. cal_pre) then
    random_clds = .true.
 elseif(newsas) then
    random_clds = .false.
 endif

 write(0,*)
 write(0,*) '--- begin random numbers:'
 seed0 = 0
 if(.not. newsas .or. cal_pre) then
    if(random_clds) then

       call mpas_pool_get_dimension(mesh,'maxran',maxran)
       call mpas_pool_get_dimension(mesh,'maxsub',maxsub)
       call mpas_pool_get_array(nggps_tbd_ddt,'rannum_tank',rannum_tank)

       maxrs  = maxran/maxsub
       lonrbm = nCellsSolve / maxsub
       if(.not.allocated(rannum)) allocate(rannum(nCellsSolve*maxrs))

       seed0 = idate(1) + idate(2) + idate(3) + idate(4)
       call random_setseed(seed0)
       call random_number(wrk)
       seed0 = seed0 + nint(wrk(1)*thousnd)
       write(0,*) '--- seed0  = ', seed0
       write(0,*) '--- wrk    = ', wrk
       write(0,*) '--- maxran = ', maxran
       write(0,*) '--- maxrs  = ', maxrs
       write(0,*) '--- maxsub = ', maxsub
       write(0,*) '--- lonrbm = ', lonrbm
       write(0,*) '--- nCellsSolve = ',nCellsSolve

       iseedl = seed0 + 1
       call random_setseed(iseedl)
       call random_number(rannum)

       do nrc = 1, maxrs
          nn = (nrc-1)*nCellsSolve
          do k = 1, maxsub
             kk = k - 1
             do iCell = 1, nCellsSolve
                iiCell = kk*lonrbm + iCell
                if (iiCell > nCellsSolve) iiCell = iiCell - nCellsSolve
                rannum_tank(nrc+kk*maxrs,iCell) = rannum(iiCell+nn)
!               if(nrc+kk*maxrs == maxran) &
!                  write(0,301) nrc+kk*maxrs,iCell,iiCell,iiCell+nn
             enddo
          enddo
       enddo
       if(allocated(rannum)) deallocate(rannum)
    endif
    write(0,*) '--- max rannum_tank = ',maxval(rannum_tank(1:maxran,1:nCellsSolve))
    write(0,*) '--- min rannum_tank = ',minval(rannum_tank(1:maxran,1:nCellsSolve))
 endif
 301 format(4(1x,i9))

!--- initialization of nrcm:
 if(ras) then
    if(fix_ncld_hr) then
       nrcm = min(nrcmax,levs-1)
    else
       nrcm = min(nrcmax,levs-1)
    endif
 else
    nrcm = 2
 endif

!--- initialization of dxmaxin, dxminin, and dxinvin:
 dxmaxin = log(1.0/(max_lon*max_lat))
 dxminin = log(1.0/(min_lon*min_lat))
 dxinvin = 1.0 / (dxmaxin-dxminin)
 write(0,*) '--- end dxinvin:'

!---
 flipv = .false.
 if(iflip == 1) flipv = .true.

!---
 sas_shal = .false.
 sas_shal = sashal .and. (.not.ras)

!--- initialization of the NGGPS physics:
 write(0,*) '--- enter subroutine nuopc_phys_init:'
 call nuopc_phys_init( &
             mdl_parameters , ntcw         , ncld         , ntoz        , ntrac            , levs         , &
             me             , lsoil        , lsm          , nmtvr       , nrcm             , levozp       , &
             lonr           , latr         , jcap         , num_p3d     , num_p2d          , npdf3d       , &
             pl_coeff       , ncw          , crtrh        , cdmbgwd     , ccwf             , dlqf         , &
             ctei_rm        , cgwf         , prslrd0      , ras         , pre_rad          , ldiag3d      , &
             lgocart        , lssav_cpl    , flipv        , old_monin   , cnvgwd           , shal_cnv     , &
             sashal         , newsas       , cal_pre      , mom4ice     , mstrat           , trans_trac   , &
             nst_fcst       , moist_adj    , thermodyn_id , sfcpress_id , gen_coord_hybrid , levr         , &
             lsidea         , pdfcld       , shcnvcw      , redrag      , hybedmf          , dspheat      , &
             dxmaxin        , dxminin      , dxinvin      , cscnv       , nctp             , ntke         , &
             do_shoc        , shocaftcnv   , ntot3d       , ntot2d      , si               , ictm         , &
             isol           , ico2         , iaer         , ialb        , iems             , iovr_sw      , &
             iovr_lw        , isubc_sw     , isubc_lw     , shoc_cld    , sas_shal         , crick_proof  , &
             ccnorm         , norad_precip , idate        , iflip       , nlunit)

!deallocate the si array:
 if(allocated(si)) deallocate(si)

!initialization of the forecast hour:
 if(.not. config_do_restart) then
    fhour = 0._PHYSKIND
 endif

!initialization of ozone climatology:
 call init_o3climatology(dminfo,mesh,nggps_atm_input,ntoz)

!initialization of surface roughness. for now, we initialize and update the surface roughness
!length using the monthly-mean climatological surface roughness length zorl12m (see call to
!subroutine physics_gfs_update_surface):
!if(.not. config_do_restart) then
!   do iCell = 1, nCellsSolve
!      if(ivgtyp(iCell) >= 1 .and. ivgtyp(iCell) <= 13) then
!         zorl(iCell) = z0_sib(ivgtyp(iCell)) * 0.01
!      elseif(ivgtyp(iCell) == 0) then
!         zorl(iCell) = 0.01
!      endif
!   enddo
!endif

!initialization of some additional surface variables (needed in mpas_atmphys_nggps_sfc_properties):
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       oro(iCell)    = zgrid(1,iCell)
       oro_uf(iCell) = oro(iCell)
!      tisfc(iCell)  = tsea(iCell)
    enddo
 endif

!initialization of some additional mesh variables (needed in mpas_atmphys_nggps_dynamic_parameters):
 call mpas_dmpar_sum_int(dminfo,nCellsSolve,nCellsTotal)
 
 if(.not. config_do_restart) then
    do iCell = 1, nCellsSolve
       sinlatCell(iCell) = sin(latCell(iCell))
       coslatCell(iCell) = cos(latCell(iCell))
       nlonsCell(iCell)  = nCellsTotal
    enddo 
 endif

!initialization of all physics function tables:
 call gfuncphys

 write(0,*) '--- end subroutine physics_init.'
 write(0,*)

 end subroutine physics_init

!=================================================================================================================
 subroutine init_dirs_forphys(mesh)
!=================================================================================================================

!inout arguments:
 type(mpas_pool_type),intent(in):: mesh

!local pointers:
 integer,pointer:: nCells
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north

!local variables:
 integer:: iCell

!---------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'east'   ,east   )
 call mpas_pool_get_array(mesh,'north'  ,north  )

!Compute unit vectors in east and north directions for each cell:
 do iCell = 1, nCells

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

!=================================================================================================================
 subroutine r3_normalize(ax, ay, az)
!=================================================================================================================
!normalizes the vector (ax, ay, az)

 real (kind=RKIND), intent(inout) :: ax, ay, az
 real (kind=RKIND) :: mi

!-----------------------------------------------------------------------------------------------------------------

 mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
 ax = ax * mi
 ay = ay * mi
 az = az * mi

 end subroutine r3_normalize

!=================================================================================================================
 end module mpas_atmphys_init
!=================================================================================================================


