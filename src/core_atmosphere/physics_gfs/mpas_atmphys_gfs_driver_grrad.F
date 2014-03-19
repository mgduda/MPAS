!==================================================================================================
 module mpas_atmphys_gfs_driver_grrad
 use mpas_dmpar
 use mpas_kind_types
 use mpas_grid_types
 use mpas_timer

 use mpas_atmphys_gfs_constants
 use mpas_atmphys_gfs_interface
 use mpas_atmphys_gfs_vars
 use mpas_atmphys_manager, only: gmt,year,month,day,hour,utc_h,julday,curr_julday
 use mpas_atmphys_utilities

!GFS physics:
 use mersenne_twister, only: random_setseed,random_index,random_stat
 use module_radiation_driver
 use module_radiation_astronomy
 use module_radlw_parameters,only: nbdlw
 use module_radsw_parameters,only: nbdsw

 implicit none
 private
 public:: allocate_radiation,   &
          deallocate_radiation, & 
          mpas_atmphys_grrad

!local variables for astronomy:
!integer,parameter:: latGFS = 1
 integer,dimension(:),allocatable:: lonsperlar,global_lats_r
 real(kind=RKIND),dimension(:),allocatable  :: xlon
 real(kind=RKIND),dimension(:,:),allocatable:: sinlat_r2,coslat_r2,xlon_r2
 real(kind=RKIND),dimension(:,:),allocatable:: coszen_r2,coszdg_r2


!>\brief MPAS top driver for the GFS radiation codes.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-01-27.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Moved the calculation of fhlwr and fhswr (frequencies of long- and short-wave radiation 
!>      calls per hour to subroutine GFS_initialize in mpas_atmphys_gfs_init.F, to compute the
!>      variables nslwr and nsswr  (number of long- and short-wave radiation calls per hour). 
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-19. 


 contains

!==================================================================================================
 subroutine allocate_radiation(mesh)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh

!--------------------------------------------------------------------------------------------------

 if(.not.allocated(slmsk)     ) allocate(slmsk(1:mesh%nCellsSolve)     )
 if(.not.allocated(hprime_rad)) allocate(hprime_rad(1:mesh%nCellsSolve))
 if(.not.allocated(fluxr)     ) allocate(fluxr(1:mesh%nCellsSolve,nfxr))
 if(.not.allocated(fcice)     ) allocate(fcice(1:mesh%nCellsSolve,1:mesh%nVertLevels)  )
 if(.not.allocated(frain)     ) allocate(frain(1:mesh%nCellsSolve,1:mesh%nVertLevels)  )
 if(.not.allocated(rrime)     ) allocate(rrime(1:mesh%nCellsSolve,1:mesh%nVertLevels)  )
 if(.not.allocated(htrlw)     ) allocate(htrlw(1:mesh%nCellsSolve,1:mesh%nVertLevels)  )
 if(.not.allocated(htrsw)     ) allocate(htrsw(1:mesh%nCellsSolve,1:mesh%nVertLevels)  )
 if(.not.allocated(cldcov)    ) allocate(cldcov(1:mesh%nCellsSolve,1:mesh%nVertLevels) )

 if(.not.allocated(htrlwb)    ) allocate(htrlwb(1:mesh%nCellsSolve,1:mesh%nVertLevels,nbdlw))
 if(.not.allocated(htrswb)    ) allocate(htrswb(1:mesh%nCellsSolve,1:mesh%nVertLevels,nbdsw))

 if(.not.allocated(lonsperlar)   ) allocate(lonsperlar(1:latr   ))
 if(.not.allocated(global_lats_r)) allocate(global_lats_r(1:latr))

 if(.not.allocated(sinlat_r2)    ) allocate(sinlat_r2(1:mesh%nCellsSolve,1:latr))
 if(.not.allocated(coslat_r2)    ) allocate(coslat_r2(1:mesh%nCellsSolve,1:latr))
 if(.not.allocated(xlon_r2)      ) allocate(xlon_r2(1:mesh%nCellsSolve,1:latr)  )
 if(.not.allocated(coszen_r2)    ) allocate(coszen_r2(1:mesh%nCellsSolve,1:latr))
 if(.not.allocated(coszdg_r2)    ) allocate(coszdg_r2(1:mesh%nCellsSolve,1:latr))

 end subroutine allocate_radiation

!==================================================================================================
 subroutine deallocate_radiation
!==================================================================================================

 if(allocated(slmsk)     ) deallocate(slmsk     )
 if(allocated(hprime_rad)) deallocate(hprime_rad)
 if(allocated(fluxr)     ) deallocate(fluxr     )
 if(allocated(fcice)     ) deallocate(fcice     )
 if(allocated(frain)     ) deallocate(frain     )
 if(allocated(rrime)     ) deallocate(rrime     )
 if(allocated(htrlw)     ) deallocate(htrlw     )
 if(allocated(htrsw)     ) deallocate(htrsw     )
 if(allocated(cldcov)    ) deallocate(cldcov    )
 if(allocated(htrlwb)    ) deallocate(htrlwb    )
 if(allocated(htrswb)    ) deallocate(htrswb    )

 if(allocated(lonsperlar)   ) deallocate(lonsperlar   )
 if(allocated(global_lats_r)) deallocate(global_lats_r)
 if(allocated(sinlat_r2)    ) deallocate(sinlat_r2    )
 if(allocated(coslat_r2)    ) deallocate(coslat_r2    )
 if(allocated(xlon_r2)      ) deallocate(xlon_r2      )
 if(allocated(coszen_r2)    ) deallocate(coszen_r2    )
 if(allocated(coszdg_r2)    ) deallocate(coszdg_r2    )

 end subroutine deallocate_radiation

!==================================================================================================
 subroutine radiation_from_MPAS(mesh,sfc_input,diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(sfc_input_type),intent(inout)   :: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables:
 integer:: iCell,iLat,k

 real(kind=RKIND),parameter :: dxmax=-16.118095651, &
       dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
 real(kind=RKIND):: work1,work2

!--------------------------------------------------------------------------------------------------

 do iCell = 1, mesh%nCellsSolve
    hprime_rad(iCell) = sfc_input%orogwd%array(1,iCell) 
    slmsk(iCell)      = sfc_input%landmask%array(iCell)
    if(sfc_input%xice%array(iCell) .gt. 0._RKIND) slmsk(iCell) = 2._RKIND
 enddo

!... Local variables needed in the call to astronomy:
 do iLat = 1, latr
    global_lats_r(iLat) = Latr
    lonsperlar(iLat)    = mesh%nCellsSolve
 enddo
 do iCell = 1, mesh%nCellsSolve
 do iLat  = 1, latr
    sinlat_r2(iCell,iLat) = sin(mesh%latCell%array(iCell))
    coslat_r2(iCell,iLat) = cos(mesh%latCell%array(iCell))
    xlon_r2(iCell,iLat)   = mesh%lonCell%array(iCell)
    coszen_r2(iCell,iLat) = 0._RKIND
    coszdg_r2(iCell,iLat) = 0._RKIND
 enddo
 enddo

!... Initialization for the Ferrier cloud microphysics scheme:
 if(num_p3d == 3) then
    do iCell = 1, mesh%nCellsSolve
       do k = 1, mesh%nVertLevels
          fcice(iCell,k) = diag_physics%phy_f3d%array(1,k,iCell)
          frain(iCell,k) = diag_physics%phy_f3d%array(2,k,iCell)
          rrime(iCell,k) = diag_physics%phy_f3d%array(3,k,iCell)
       enddo
    !... weighs range of minimum large ice fraction: 
       do iLat = 1, latr
          work1 = (log(coslat_r2(iCell,iLat)/(mesh%nCellsSolve))-dxmin)*dxinv
          work1 = max(0.0, min(1.0,work1))
          work2 = flgmin(1) * work1 + flgmin(2) * work2
       enddo
!      diag_physics%flgmin_2d%array(iCell) = work2
       diag_physics%flgmin_2d%array(iCell) = 0._RKIND
    enddo
 else
    do iCell = 1, mesh%nCellsSolve
       do k = 1, mesh%nVertLevels
          fcice(iCell,k) = 0._RKIND
          frain(iCell,k) = 0._RKIND
          rrime(iCell,k) = 0._RKIND
       enddo
       diag_physics%flgmin_2d%array(iCell) = 0._RKIND
    enddo
 endif

!... Initialization of output arrays in grrad:
 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       htrlw(iCell,k)  = tend_physics%rthratenlw%array(k,iCell)
       htrsw(iCell,k)  = tend_physics%rthratensw%array(k,iCell)
       cldcov(iCell,k) = diag_physics%cldfrac   %array(k,iCell)
    enddo
    fluxr(iCell,:) = 0._RKIND
!   diag_physics%sfcnsw % array(iCell) = 0._RKIND
!   diag_physics%sfcdsw % array(iCell) = 0._RKIND
!   diag_physics%sfcdlw % array(iCell) = 0._RKIND   
!   diag_physics%sfcalb % array(iCell) = 0._RKIND
!   diag_physics%sfcemis% array(iCell) = 0._RKIND
!   diag_physics%tsflw  % array(iCell) = 0._RKIND
 enddo

 end subroutine radiation_from_MPAS

!==================================================================================================
 subroutine radiation_to_MPAS(mesh,diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables:
 integer:: iCell,k

!--------------------------------------------------------------------------------------------------

 do iCell = 1, mesh%nCellsSolve
    do k = 1, mesh%nVertLevels
       tend_physics%rthratenlw%array(k,iCell) = htrlw(iCell,k)
       tend_physics%rthratensw%array(k,iCell) = htrsw(iCell,k)
       diag_physics%cldfrac   %array(k,iCell) = cldcov(iCell,k)
    enddo
 enddo
 write(0,*)
 write(0,*) '--- end subroutine radiation_to_MPAS:'
 write(0,*) '--- max rthratenlw = ', maxval(tend_physics%rthratenlw%array(:,1:mesh%nCellsSolve))
 write(0,*) '--- min rthratenlw = ', minval(tend_physics%rthratenlw%array(:,1:mesh%nCellsSolve))
 write(0,*) '--- max rthratensw = ', maxval(tend_physics%rthratensw%array(:,1:mesh%nCellsSolve))
 write(0,*) '--- min rthratensw = ', minval(tend_physics%rthratensw%array(:,1:mesh%nCellsSolve))
 write(0,*)

 end subroutine radiation_to_MPAS

!==================================================================================================
 subroutine mpas_atmphys_grrad(itimestep,xtime_s,dminfo,mesh,sfc_input,diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mesh_type),intent(in):: mesh

 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(sfc_input_type),intent(inout)   :: sfc_input
 type(diag_physics_type),intent(inout):: diag_physics
 type(tend_physics_type),intent(inout):: tend_physics

!local variables:
 logical:: lprnt,sas_shal

 integer:: nCellsSolve,nVertLevels
 integer:: i1,i2,k
 integer:: icwp
 integer:: me
 integer:: iCell,iLat

 real(kind=RKIND):: solcon

!local variables for time manager:
 integer:: iddat,ihdat,imdat,iydat
 integer:: idat(8),jdat(8),jdate(4)
 real(kind=RKIND):: rinc(5)
 real(kind=RKIND):: raddt

!local variables for 2-d random seedd generator:
!type(random_stat):: stat
 integer, parameter:: ipsdlim = 1.0e8 ! upper limit for random seeds

 integer:: ipseed,ipsd0
 integer:: numrdm(lonr*latr,2)

!local variables for radiation updates:
 integer:: iaeros

 real(kind=RKIND):: xtime_h

!local variables that are not used, and for which I could not find any initialization:
 integer:: k1oz,k2oz
 real(kind=RKIND):: facoz
 real(kind=RKIND),dimension(mesh%nVertLevels+1):: sigma

!--------------------------------------------------------------------------------------------------
 write(0,100) itimestep

!formats:
 100 format(/,' --- enter subroutine mpas_atmphys_grrad ',i6) 
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

!... initialization of local radiation arrays:
 call radiation_from_MPAS(mesh,sfc_input,diag_physics,tend_physics)

!... date conversion:
 fhour    = utc_h

 jdate(1) = hour
 jdate(2) = month
 jdate(3) = day
 jdate(4) = year

 idat(:)  = 0
 idat(1)  = jdate(4)
 idat(2)  = jdate(2)
 idat(3)  = jdate(3)
 idat(4)  = 0
 idat(5)  = jdate(1)

 jdat(:)  = 0
 rinc(:)  = 0._RKIND
 rinc(2)  = fhour
 call w3movdat(rinc,idat,jdat) 

!... initialization of local variables:
 nCellsSolve = mesh % nCellsSolve
 nVertLevels = mesh % nVertLevels
!if(dminfo % my_proc_id == IO_NODE) me = 0
 me = dminfo % my_proc_id
 sas_shal = sashal .and. (.not. ras)

 if(ntcw > 0) icwp = 1

 lprnt = .false.

!... initialization of local sigma levels (I believe that they are not used anywhere):
 do k = 1, nVertLevels+1
    sigma(k) = 1.0_RKIND
 enddo

!... initialization of radiation:
!write(0,*) 'iflip    =', iflip
!write(0,*) 'ictm     =', ictm
!write(0,*) 'isol     =', isol
!write(0,*) 'ico2     =', ico2
!write(0,*) 'iaer     =', iaer
!write(0,*) 'ialb     =', ialb
!write(0,*) 'iems     =', iems
!write(0,*) 'icwp     =', icwp
!write(0,*) 'num_p3d  =', num_p3d
!write(0,*) 'isubc_sw =', isubc_sw
!write(0,*) 'isubc_lw =', isubc_lw
!write(0,*) 'me       =', me
 call radinit(sigma,nVertLevels,iflip,idat,jdat,ictm,isol,ico2,iaer,ialb,iems,icwp,num_p3d, &
              isubc_sw,isubc_lw,iovr_sw,iovr_lw,me)

!... initialization for climatological ozone (not used):
 if(ntoz .le. 0) &
    call physics_error_fatal('--- climatological O3 option no longer supported. Reset ntoz to 2.')
 k1oz  = 1
 k2oz  = 1
 facoz = 1._RKIND

!... calculation of astronomy for shortwave radiation code:
 lsswr = l_radtsw
 lslwr = l_radtlw
 raddt = min(dt_radtsw,dt_radtlw)

 write(0,*) 'fhlwr =', fhlwr
 write(0,*) 'dtlw  =', dt_radtlw
 write(0,*) 'lslwr =', lslwr
 write(0,*)
 write(0,*) 'fhswr =', fhswr
 write(0,*) 'dtsw  =', dt_radtsw
 write(0,*) 'lsswr =', lsswr

!... calculate the solar zenith angle:
 call astronomy(lonsperlar,global_lats_r,sinlat_r2,coslat_r2,xlon_r2,fhswr,jdat, &
            nCellsSolve,latr,latr,latr,lsswr,me,solcon,slag,sdec,cdec,coszen_r2, &
            coszdg_r2)
!write(0,*)
!write(0,*) 'solcon =', solcon
!write(0,*) 'slag   =', slag
!write(0,*) 'sdec   =', sdec
!write(0,*) 'cdec   =', cdec

!write(0,*)
 do iCell = 1, nCellsSolve
 do iLat  = 1, latr
    diag_physics%coszen%array(iCell) = coszen_r2(iCell,iLat)
    diag_physics%coszdg%array(iCell) = coszdg_r2(iCell,iLat)
!   if(abs(diag_physics%coszen%array(iCell)) .gt. 0._RKIND) &
!      write(0,102) iCell,diag_physics%coszen%array(iCell),diag_physics%coszdg%array(iCell)
 enddo
 enddo

!... generate 2-d random seed array for sub-grid clouds-radiation interactions:
!if(isubc_lw == 2 .or. isubc_sw ==2 ) then
!   ipseed = mod(nint(100.0*sqrt(fhour*3600)),ipsdlim) + 1 + ipsd0
!   call random_setseed(ipseed,stat)
!   write(0,*) '--- ipsd0  =', ipsd0
!   write(0,*) '--- ipseed =', ipseed
!   write(0,*)
!   call random_index(ipsdlim,numrdm,stat)
!endif
!stop

!... CALL GFS PHYSICS DRIVER FOR ALL PHYSICS, EXCEPT RADIATION:
!write(0,*)
!write(0,*) '--- before call to grrad:'
!write(0,*) '... itimestep    = ', itimestep
!write(0,*) '... nCellsSolve  = ', nCellsSolve
!write(0,*) '... nVertLevels  = ', nVertLevels
!write(0,*) '... norad_precip = ', norad_precip
!write(0,*) '... crick_proof  = ', crick_proof
!write(0,*) '... ccnorm       = ', ccnorm
!write(0,*) '... lssav        = ', lssav
!write(0,*) '... sas_shal     = ', sas_shal
!write(0,*) '... lslwr        = ', lslwr
!write(0,*) '... lsswr        = ', lsswr
!write(0,*) '... iflip        = ', iflip
!write(0,*) '... me           = ', me
!write(0,*) '... lprnt        = ', lprnt
!write(0,*) '... ipt          = ', nCellsSolve
!write(0,*) '... np3d         = ', num_p3d
!write(0,*) '... nfxr         = ', nfxr
!write(0,*) '... ntcw         = ', ntcw
!write(0,*) '... ncld         = ', ncld
!write(0,*) '... ntoz         = ', ntoz
!write(0,*) '... ntrac        = ', ntrac
!write(0,*) '... k1oz         = ', k1oz
!write(0,*) '... k2oz         = ', k2oz
!write(0,*) '... facoz        = ', facoz
!write(0,*) '... iovrlw       = ', iovr_lw
!write(0,*) '... iovrsw       = ', iovr_sw
!write(0,*) '... dtlw         = ', dt_radtlw
!write(0,*) '... dtsw         = ', dt_radtsw
!write(0,*) '... solcon       = ', solcon
!write(0,*)
!iCell = nCellsSolve
!do k = 1, nVertLevels+1
!   write(0,103) iCell,k,prsi_rad(iCell,k)
!enddo
!write(0,*)
!do k = 1, nVertLevels
!   write(0,103) iCell,k,prsl_rad(iCell,k),prslk(iCell,k),tgrs(iCell,k),qvrs(iCell,k), & 
!      qgrs(iCell,k,1),qgrs(iCell,k,2),qgrs(iCell,k,3)
!enddo
!write(0,*)
!do k = 1, nVertLevels
!   write(0,103) iCell,k,vvel(iCell,k)
!enddo
!write(0,*)
!do iCell = 1, nCellsSolve
!   write(0,102) iCell,slmsk(iCell),mesh%lonCell%array(iCell),mesh%latCell%array(iCell), &
!         sfc_input%skintemp%array(iCell),sfc_input%snow%array(iCell),sfc_input%snoalb%array(iCell), &
!         sfc_input%facsf%array(iCell),sfc_input%facwf%array(iCell)
!enddo
!write(0,*)
!do iCell = 1, nCellsSolve
!   write(0,102) iCell,diag_physics%zorl%array(iCell),     &
!         diag_physics%alnsf%array(iCell),diag_physics%alnwf%array(iCell), &
!         diag_physics%alvsf%array(iCell),diag_physics%alvwf%array(iCell)
!enddo
!write(0,*)

 call grrad( &
!... inputs:
  kdt          = itimestep,    & ! done.
  ix           = nCellsSolve,  & ! done.
  im           = nCellsSolve,  & ! done.
  lm           = nVertLevels,  & ! done.
  norad_precip = norad_precip, & ! done.
  crick_proof  = crick_proof,  & ! done.
  ccnorm       = ccnorm,       & ! done.
  lssav        = lssav,        & ! done.
  sashal       = sas_shal,     & ! done.
  lslwr        = lslwr,        & ! done.
  lsswr        = lsswr,        & ! done.
  iflip        = iflip,        & ! done. 
  me           = me,           & ! done.
  lprnt        = lprnt,        & ! done.
  ipt          = nCellsSolve,  & ! done.
  np3d         = num_p3d,      & ! done.
  nfxr         = nfxr,         & ! done.
  ntcw         = ntcw,         & ! done.
  ncld         = ncld,         & ! done.
  ntoz         = ntoz,         & ! done.
  ntrac        = ntrac,        & ! done.
  k1oz         = k1oz,         & ! done.
  k2oz         = k2oz,         & ! done.
  facoz        = facoz,        & ! done.
  iovrlw       = iovr_lw,      & ! done.
  iovrsw       = iovr_sw,      & ! done.
  dtlw         = dt_radtlw,    & ! done.
  dtsw         = dt_radtsw,    & ! done.
  solcon       = solcon,       & ! done.
  prsi         = prsi_rad,     & ! done.
  prsl         = prsl_rad,     & ! done.
  prslk        = prslk,        & 
  tgrs         = tgrs,         & ! done.
  qgrs         = qvrs,         & ! done.
  oz           = qgrs,         &
  vvl          = vvel,         &
  slmsk        = slmsk,        & ! done.
  xlon         = mesh%lonCell%array(1:mesh%nCellsSolve),           & ! done.
  xlat         = mesh%latCell%array(1:mesh%nCellsSolve),           & ! done.
  tsfc         = sfc_input%skintemp%array(1:mesh%nCellsSolve),     & ! done.
  snowd        = sfc_input%snow%array(1:mesh%nCellsSolve),         & ! done.
  sncovr       = sfc_input%snowc%array(1:mesh%nCellsSolve),        & ! done.
  snoalb       = sfc_input%snoalb%array(1:mesh%nCellsSolve),       & ! done.
 !hprim        = sfc_input%orogwd%array(1,1:mesh%nCellsSolve),     & ! done.
  hprim        = hprime_rad,                                       & ! done.
  facsf        = sfc_input%facsf%array(1:mesh%nCellsSolve),        & ! done.
  facwf        = sfc_input%facwf%array(1:mesh%nCellsSolve),        & ! done.
  fice         = sfc_input%xice%array(1:mesh%nCellsSolve),         & ! done.
  alvsf        = diag_physics%alvsf%array(1:mesh%nCellsSolve),     & ! done.
  alnsf        = diag_physics%alnsf%array(1:mesh%nCellsSolve),     & ! done.
  alvwf        = diag_physics%alvwf%array(1:mesh%nCellsSolve),     & ! done.
  alnwf        = diag_physics%alnwf%array(1:mesh%nCellsSolve),     & ! done.
  zorl         = diag_physics%zorl%array(1:mesh%nCellsSolve),      & ! done.
  coszen       = diag_physics%coszen%array(1:mesh%nCellsSolve),    & ! done.
  coszdg       = diag_physics%coszdg%array(1:mesh%nCellsSolve),    & ! done.
  cv           = diag_physics%cv%array(1:mesh%nCellsSolve),        & ! done.
  cvb          = diag_physics%cvb%array(1:mesh%nCellsSolve),       & ! done.
  cvt          = diag_physics%cvt%array(1:mesh%nCellsSolve),       & ! done.
  tisfc        = diag_physics%tisfc%array(1:mesh%nCellsSolve),     & ! done.
  flgmin       = diag_physics%flgmin_2d%array(1:mesh%nCellsSolve), & ! done.
  icsdsw       = diag_physics%icsdsw%array(1:mesh%nCellsSolve),    & ! done.
  icsdlw       = diag_physics%icsdlw%array(1:mesh%nCellsSolve),    & ! done.
  fcice        = fcice,  &                                           ! done.
  frain        = frain,  &                                           ! done.
  rrime        = rrime,  &                                           ! done.
!... inouts:
  fluxr        = fluxr,  &
!... outputs:
  htrswb       = htrswb, &
  htrlwb       = htrlwb, &
  htrsw        = htrsw,  &                                           ! done.
  htrlw        = htrlw,  &                                           ! done.
  cldcov       = cldcov, &                                           ! done.
  sfcnsw       = diag_physics%sfcnsw%array(1:nCellsSolve),         & ! done.
  sfcdsw       = diag_physics%sfcdsw%array(1:nCellsSolve),         & ! done.
  sfcdlw       = diag_physics%sfcdlw%array(1:nCellsSolve),         & ! done.
  tsflw        = diag_physics%tsflw%array(1:nCellsSolve),          & ! done.
  sfalb        = diag_physics%sfcalb%array(1:nCellsSolve),         & ! done.
  semis        = diag_physics%sfcemis%array(1:nCellsSolve)         & ! done. 
           )
!write(0,*)
!do iCell = 1,mesh%nCellsSolve
!do k = 1, mesh%nVertLevels
!   write(0,103) iCell,k,htrlw(iCell,k),htrsw(iCell,k)
!enddo
!write(0,*)
!enddo
 write(0,*) '--- end subroutine grrad:'

!copy all arrays back to MPAS grid:
 call radiation_to_MPAS(mesh,diag_physics,tend_physics)
 write(0,*) '--- end subroutine mpas_atmphys_grrad'

!formats:
 102 format(i8,10(1x,e15.8))
 103 format(i8,i3,10(1x,e15.8))

 end subroutine mpas_atmphys_grrad

!==================================================================================================
 end module mpas_atmphys_gfs_driver_grrad
!==================================================================================================
