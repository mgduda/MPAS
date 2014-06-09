!=================================================================================================================
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
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine allocate_radiation(mesh)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 if(.not.allocated(slmsk)     ) allocate(slmsk(1:nCellsSolve)     )
 if(.not.allocated(hprime_rad)) allocate(hprime_rad(1:nCellsSolve))
 if(.not.allocated(fluxr)     ) allocate(fluxr(1:nCellsSolve,nfxr))
 if(.not.allocated(fcice)     ) allocate(fcice(1:nCellsSolve,1:nVertLevels) )
 if(.not.allocated(frain)     ) allocate(frain(1:nCellsSolve,1:nVertLevels) )
 if(.not.allocated(rrime)     ) allocate(rrime(1:nCellsSolve,1:nVertLevels) )
 if(.not.allocated(htrlw)     ) allocate(htrlw(1:nCellsSolve,1:nVertLevels) )
 if(.not.allocated(htrsw)     ) allocate(htrsw(1:nCellsSolve,1:nVertLevels) )
 if(.not.allocated(cldcov)    ) allocate(cldcov(1:nCellsSolve,1:nVertLevels))

 if(.not.allocated(htrlwb)    ) allocate(htrlwb(1:nCellsSolve,1:nVertLevels,nbdlw))
 if(.not.allocated(htrswb)    ) allocate(htrswb(1:nCellsSolve,1:nVertLevels,nbdsw))

 if(.not.allocated(lonsperlar)   ) allocate(lonsperlar(1:latr   ))
 if(.not.allocated(global_lats_r)) allocate(global_lats_r(1:latr))

 if(.not.allocated(sinlat_r2)    ) allocate(sinlat_r2(1:nCellsSolve,1:latr))
 if(.not.allocated(coslat_r2)    ) allocate(coslat_r2(1:nCellsSolve,1:latr))
 if(.not.allocated(xlon_r2)      ) allocate(xlon_r2(1:nCellsSolve,1:latr)  )
 if(.not.allocated(coszen_r2)    ) allocate(coszen_r2(1:nCellsSolve,1:latr))
 if(.not.allocated(coszdg_r2)    ) allocate(coszdg_r2(1:nCellsSolve,1:latr))

 end subroutine allocate_radiation

!=================================================================================================================
 subroutine deallocate_radiation
!=================================================================================================================

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

!=================================================================================================================
 subroutine radiation_from_MPAS(mesh,sfc_input,diag_physics,tend_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: xice
 real(kind=RKIND),dimension(:),pointer    :: flgmin_2d
 real(kind=RKIND),dimension(:),pointer    :: sfcnsw,sfcdsw,sfcdlw,sfcalb,sfcemis,tsflw
 real(kind=RKIND),dimension(:,:),pointer  :: orogwd,rthratenlw,rthratensw,cldfrac
 real(kind=RKIND),dimension(:,:,:),pointer:: phy_f3d

!local variables:
 integer:: iCell,iLat,k

 real(kind=RKIND),parameter :: dxmax=-16.118095651, &
       dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)
 real(kind=RKIND):: work1,work2

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)

 call mpas_pool_get_array(sfc_input   ,'landmask' ,landmask )
 call mpas_pool_get_array(sfc_input   ,'orogwd'   ,orogwd   )
 call mpas_pool_get_array(sfc_input   ,'xice'     ,xice     )
 call mpas_pool_get_array(diag_physics,'phy_f3d'  ,phy_f3d  )
 call mpas_pool_get_array(diag_physics,'flgmin_2d',flgmin_2d)
 call mpas_pool_get_array(diag_physics,'sfcnsw'   ,sfcnsw   )
 call mpas_pool_get_array(diag_physics,'sfcdsw'   ,sfcdsw   )
 call mpas_pool_get_array(diag_physics,'sfcdlw'   ,sfcdlw   )
 call mpas_pool_get_array(diag_physics,'sfcalb'   ,sfcalb   )
 call mpas_pool_get_array(diag_physics,'sfcemis'  ,sfcemis  )
 call mpas_pool_get_array(diag_physics,'tsflw'    ,tsflw    )
 call mpas_pool_get_array(diag_physics,'cldfrac'  ,cldfrac  )

 call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
 call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

 do iCell = 1, nCellsSolve
    hprime_rad(iCell) = orogwd(1,iCell) 
    slmsk(iCell)      = landmask(iCell)
    if(xice(iCell) .gt. 0._RKIND) slmsk(iCell) = 2._RKIND
 enddo

!... Local variables needed in the call to astronomy:
 do iLat = 1, latr
    global_lats_r(iLat) = Latr
    lonsperlar(iLat)    = nCellsSolve
 enddo
 do iCell = 1, nCellsSolve
 do iLat  = 1, latr
    sinlat_r2(iCell,iLat) = sin(latCell(iCell))
    coslat_r2(iCell,iLat) = cos(latCell(iCell))
    xlon_r2(iCell,iLat)   = lonCell(iCell)
    coszen_r2(iCell,iLat) = 0._RKIND
    coszdg_r2(iCell,iLat) = 0._RKIND
 enddo
 enddo

!... Initialization for the Ferrier cloud microphysics scheme:
 if(num_p3d == 3) then
    do iCell = 1, nCellsSolve
       do k = 1, nVertLevels
          fcice(iCell,k) = phy_f3d(1,k,iCell)
          frain(iCell,k) = phy_f3d(2,k,iCell)
          rrime(iCell,k) = phy_f3d(3,k,iCell)
       enddo
    !... weighs range of minimum large ice fraction: 
       do iLat = 1, latr
          work1 = (log(coslat_r2(iCell,iLat)/nCellsSolve)-dxmin)*dxinv
          work1 = max(0.0, min(1.0,work1))
          work2 = flgmin(1) * work1 + flgmin(2) * work2
       enddo
!      flgmin_2d(iCell) = work2
       flgmin_2d(iCell) = 0._RKIND
    enddo
 else
    do iCell = 1, nCellsSolve
       do k = 1, nVertLevels
          fcice(iCell,k) = 0._RKIND
          frain(iCell,k) = 0._RKIND
          rrime(iCell,k) = 0._RKIND
       enddo
       flgmin_2d(iCell) = 0._RKIND
    enddo
 endif

!... Initialization of output arrays in grrad:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       htrlw(iCell,k)  = rthratenlw(k,iCell)
       htrsw(iCell,k)  = rthratensw(k,iCell)
       cldcov(iCell,k) = cldfrac(k,iCell)
    enddo
    fluxr(iCell,:) = 0._RKIND
!   sfcnsw(iCell)  = 0._RKIND
!   sfcdsw(iCell)  = 0._RKIND
!   sfcdlw(iCell)  = 0._RKIND   
!   sfcalb(iCell)  = 0._RKIND
!   sfcemis(iCell) = 0._RKIND
!   tsflw(iCell)   = 0._RKIND
 enddo

 end subroutine radiation_from_MPAS

!=================================================================================================================
 subroutine radiation_to_MPAS(mesh,diag_physics,tend_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw,cldfrac

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(diag_physics,'cldfrac'   ,cldfrac   )
 call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
 call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       rthratenlw(k,iCell) = htrlw(iCell,k)
       rthratensw(k,iCell) = htrsw(iCell,k)
       cldfrac(k,iCell)    = cldcov(iCell,k)
    enddo
 enddo
 write(0,*)
 write(0,*) '--- end subroutine radiation_to_MPAS:'
 write(0,*) '--- max rthratenlw = ', maxval(rthratenlw(:,1:nCellsSolve))
 write(0,*) '--- min rthratenlw = ', minval(rthratenlw(:,1:nCellsSolve))
 write(0,*) '--- max rthratensw = ', maxval(rthratensw(:,1:nCellsSolve))
 write(0,*) '--- min rthratensw = ', minval(rthratensw(:,1:nCellsSolve))
 write(0,*)

 end subroutine radiation_to_MPAS

!=================================================================================================================
 subroutine mpas_atmphys_grrad(itimestep,xtime_s,dminfo,mesh,sfc_input,diag_physics,tend_physics)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,dimension(:),pointer:: icsdlw,icsdsw

 real(kind=RKIND),dimension(:),pointer:: lonCell,latCell
 real(kind=RKIND),dimension(:),pointer:: skintemp,snow,snowc,snoalb,xice,flgmin_2d
 real(kind=RKIND),dimension(:),pointer:: facsf,facwf,alvsf,alnsf,alvwf,alnwf
 real(kind=RKIND),dimension(:),pointer:: zorl,coszen,coszdg,cv,cvb,cvt,tisfc
 real(kind=RKIND),dimension(:),pointer:: sfcnsw,sfcdsw,sfcdlw,tsflw,sfcalb,sfcemis


!local variables:
 logical:: lprnt,sas_shal

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
 real(kind=RKIND),dimension(:),allocatable:: sigma

!-----------------------------------------------------------------------------------------------------------------
 write(0,100) itimestep

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'latCell',latCell)

 call mpas_pool_get_array(sfc_input   ,'skintemp' ,skintemp )
 call mpas_pool_get_array(sfc_input   ,'snow'     ,snow     )
 call mpas_pool_get_array(sfc_input   ,'snowc'    ,snowc    )
 call mpas_pool_get_array(sfc_input   ,'snoalb'   ,snoalb   )
 call mpas_pool_get_array(sfc_input   ,'facsf'    ,facsf    )
 call mpas_pool_get_array(sfc_input   ,'facwf'    ,facwf    )
 call mpas_pool_get_array(sfc_input   ,'xice'     ,xice     )
 call mpas_pool_get_array(diag_physics,'alvsf'    ,alvsf    ) 
 call mpas_pool_get_array(diag_physics,'alnsf'    ,alnsf    )
 call mpas_pool_get_array(diag_physics,'alvwf'    ,alvwf    )
 call mpas_pool_get_array(diag_physics,'alnwf'    ,alnwf    )
 call mpas_pool_get_array(diag_physics,'coszen'   ,coszen   ) 
 call mpas_pool_get_array(diag_physics,'coszdg'   ,coszdg   )
 call mpas_pool_get_array(diag_physics,'cv'       ,cv       )
 call mpas_pool_get_array(diag_physics,'cvb'      ,cvb      )
 call mpas_pool_get_array(diag_physics,'cvt'      ,cvt      )
 call mpas_pool_get_array(diag_physics,'icsdlw'   ,icsdlw   )
 call mpas_pool_get_array(diag_physics,'icsdsw'   ,icsdsw   )
 call mpas_pool_get_array(diag_physics,'sfcalb'   ,sfcalb   )
 call mpas_pool_get_array(diag_physics,'sfcemis'  ,sfcemis  )
 call mpas_pool_get_array(diag_physics,'zorl'     ,zorl     )
 call mpas_pool_get_array(diag_physics,'tisfc'    ,tisfc    )
 call mpas_pool_get_array(diag_physics,'flgmin_2d',flgmin_2d)
 call mpas_pool_get_array(diag_physics,'sfcnsw'   ,sfcnsw   )
 call mpas_pool_get_array(diag_physics,'sfcdsw'   ,sfcdsw   )
 call mpas_pool_get_array(diag_physics,'sfcdlw'   ,sfcdlw   )
 call mpas_pool_get_array(diag_physics,'tsflw'    ,tsflw    )

 if(.not.allocated(sigma)) allocate(sigma(nVertLevels+1))

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
    coszen(iCell) = coszen_r2(iCell,iLat)
    coszdg(iCell) = coszdg_r2(iCell,iLat)
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

!... CALL GFS LONG- AND SHORT-WAVE RADIATION CODES:
 call grrad( &
!... inputs:
  kdt          = itimestep                , ix           = nCellsSolve             , &
  im           = nCellsSolve              , lm           = nVertLevels             , & 
  norad_precip = norad_precip             , crick_proof  = crick_proof             , &
  ccnorm       = ccnorm                   , lssav        = lssav                   , &
  sashal       = sas_shal                 , lslwr        = lslwr                   , &
  lsswr        = lsswr                    , iflip        = iflip                   , & 
  me           = me                       , lprnt        = lprnt                   , &
  ipt          = nCellsSolve              , np3d         = num_p3d                 , &
  nfxr         = nfxr                     , ntcw         = ntcw                    , &
  ncld         = ncld                     , ntoz         = ntoz                    , &
  ntrac        = ntrac                    , k1oz         = k1oz                    , &
  k2oz         = k2oz                     , facoz        = facoz                   , &
  iovrlw       = iovr_lw                  , iovrsw       = iovr_sw                 , &
  dtlw         = dt_radtlw                , dtsw         = dt_radtsw               , &
  solcon       = solcon                   , prsi         = prsi_rad                , &
  prsl         = prsl_rad                 , prslk        = prslk                   , & 
  tgrs         = tgrs                     , qgrs         = qvrs                    , &
  oz           = qgrs                     , vvl          = vvel                    , &
  slmsk        = slmsk                    , xlon         = lonCell(1:nCellsSolve)  , &
  xlat         = latCell(1:nCellsSolve)   , tsfc         = skintemp(1:nCellsSolve) , &
  snowd        = snow(1:nCellsSolve)      , sncovr       = snowc(1:nCellsSolve)    , &
  snoalb       = snoalb(1:nCellsSolve)    , hprim        = hprime_rad              , &
  facsf        = facsf(1:nCellsSolve)     , facwf        = facwf(1:nCellsSolve)    , &
  fice         = xice(1:nCellsSolve)      , alvsf        = alvsf(1:nCellsSolve)    , &
  alnsf        = alnsf(1:nCellsSolve)     , alvwf        = alvwf(1:nCellsSolve)    , &
  alnwf        = alnwf(1:nCellsSolve)     , zorl         = zorl(1:nCellsSolve)     , &
  coszen       = coszen(1:nCellsSolve)    , coszdg       = coszdg(1:nCellsSolve)   , &
  cv           = cv(1:nCellsSolve)        , cvb          = cvb(1:nCellsSolve)      , &
  cvt          = cvt(1:nCellsSolve)       , tisfc        = tisfc(1:nCellsSolve)    , &
  flgmin       = flgmin_2d(1:nCellsSolve) , icsdsw       = icsdsw(1:nCellsSolve)   , &
  icsdlw       = icsdlw(1:nCellsSolve)    , fcice        = fcice                   , &
  frain        = frain                    , rrime        = rrime                   , & 
!... inouts:
  fluxr        = fluxr ,                                                             &
!... outputs:
  htrswb       = htrswb                   , htrlwb       = htrlwb                  , &
  htrsw        = htrsw                    , htrlw        = htrlw                   , &
  cldcov       = cldcov                   , sfcnsw       = sfcnsw(1:nCellsSolve)   , &
  sfcdsw       = sfcdsw(1:nCellsSolve)    , sfcdlw       = sfcdlw(1:nCellsSolve)   , &
  tsflw        = tsflw(1:nCellsSolve)     , sfalb        = sfcalb(1:nCellsSolve)   , &
  semis        = sfcemis(1:nCellsSolve)                                              &
           )
 write(0,*) '--- end subroutine grrad:'

!copy all arrays back to MPAS grid:
 call radiation_to_MPAS(mesh,diag_physics,tend_physics)
 write(0,*) '--- end subroutine mpas_atmphys_grrad'

 if(allocated(sigma)) deallocate(sigma)

!formats:
 102 format(i8,10(1x,e15.8))
 103 format(i8,i3,10(1x,e15.8))

 end subroutine mpas_atmphys_grrad

!=================================================================================================================
 end module mpas_atmphys_gfs_driver_grrad
!=================================================================================================================
