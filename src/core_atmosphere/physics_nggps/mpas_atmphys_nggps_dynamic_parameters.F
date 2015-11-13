! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_dynamic_parameters
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_manager,only: gmt,year,month,day,hour,minute,second,msecond,utc_h,julday,curr_julday
 use mpas_atmphys_utilities
 use mpas_atmphys_nggps_vars

 use nuopc_physics
 use mersenne_twister
 use module_radiation_driver

 implicit none
 private
 public:: mpas_dynamic_parameters

 type,extends(dynamic_parameters):: mpas_dynamic_parameters

       integer,dimension(:),pointer:: nlons_l  => null()
       integer,dimension(:),pointer:: icsdsw_l => null()!auxiliary cloud control arrays passed to main radiation.
       integer,dimension(:),pointer:: icsdlw_l => null()!if isubcsw/isubclw (input to init) are set to 2, these
                                                        !arrays contain provided random seeds for sub-column
                                                        !clouds generators.

       real(kind=kind_phys),dimension(:),pointer:: xlon_l   => null()!grid longitude in radians.
       real(kind=kind_phys),dimension(:),pointer:: xlat_l   => null()!grid latitude in radians.
       real(kind=kind_phys),dimension(:),pointer:: sinlat_l => null()!sine of the grids corresponding latitudes.
       real(kind=kind_phys),dimension(:),pointer:: coslat_l => null()!cosine of the grids corresponding latitudes.

    contains
       procedure:: mpas_setrad    => mpas_dyn_param_setrad
       procedure:: mpas_setphys   => mpas_dyn_param_setphys
       procedure:: allocate_mem   => mpas_dyn_param_allocate
       procedure:: deallocate_mem => mpas_dyn_param_deallocate

 end type

 real(kind=RKIND):: cdec,sdec,slag

 real(kind=RKIND),parameter:: cons_24 = 24._RKIND
 real(kind=RKIND),parameter:: cons_99 = 99._RKIND

 contains


!=================================================================================================================
 subroutine mpas_dyn_param_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_dynamic_parameters),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine mpas_dyn_param_allocate:'

 call mpas_pool_get_dimension(mesh,"nCellsSolve",nCellsSolve)

 if(.not.associated(this%nlons_l) ) allocate(this%nlons_l(1:nCellsSolve) )
 if(.not.associated(this%xlon_l)  ) allocate(this%xlon_l(1:nCellsSolve)  )
 if(.not.associated(this%xlat_l)  ) allocate(this%xlat_l(1:nCellsSolve)  )
 if(.not.associated(this%sinlat_l)) allocate(this%sinlat_l(1:nCellsSolve))
 if(.not.associated(this%coslat_l)) allocate(this%coslat_l(1:nCellsSolve))
 if(.not.associated(this%icsdsw_l)) allocate(this%icsdsw_l(1:nCellsSolve))
 if(.not.associated(this%icsdlw_l)) allocate(this%icsdlw_l(1:nCellsSolve))

 write(0,*) '--- end subroutine mpas_dyn_param_allocate'

 end subroutine mpas_dyn_param_allocate

!=================================================================================================================
 subroutine mpas_dyn_param_deallocate(this)
!=================================================================================================================

 class(mpas_dynamic_parameters),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine mpas_dyn_param_deallocate:'

 if(associated(this%nlons_l) ) deallocate(this%nlons_l )
 if(associated(this%xlon_l)  ) deallocate(this%xlon_l  )
 if(associated(this%xlat_l)  ) deallocate(this%xlat_l  )
 if(associated(this%sinlat_l)) deallocate(this%sinlat_l)
 if(associated(this%coslat_l)) deallocate(this%coslat_l)
 if(associated(this%icsdsw_l)) deallocate(this%icsdsw_l)
 if(associated(this%icsdlw_l)) deallocate(this%icsdlw_l)

 write(0,*) '--- end subroutine mpas_dyn_param_deallocate'

 end subroutine mpas_dyn_param_deallocate

!=================================================================================================================
 subroutine mpas_dyn_param_setrad(this,configs,mesh,itimestep,xtime_s)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 class(mpas_dynamic_parameters),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve

 real(kind=RKIND),pointer:: config_dt
 real(kind=RKIND),dimension(:),pointer:: lonCell,latCell

 integer,dimension(:),pointer:: icsdlw_l,icsdsw_l
 real(kind=kind_phys),dimension(:),pointer:: xlon_l,xlat_l,sinlat_l,coslat_l

!local variables:
 logical:: lprnt,lssav,lslwr,lsswr
 integer:: iCell,ipt,im,ix,kdt,n
 integer:: me,latr,lonr

 real(kind=RKIND):: deltim,dtlw,dtsw,fhour,xtime_h
 real(kind=RKIND):: solcon,solhr
!real(kind=RKIND):: slag,sdec,cdec,solcon,solhr

!variables used for random number generator:
 type(random_stat):: stat

 integer,parameter:: ipsdlim = 1.0e8 !upper limit for random seeds.
 integer:: ipseed
 integer,dimension(:),allocatable  :: numrdm
 integer,dimension(:,:,:),allocatable:: ixseed

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_config(configs,'config_dt',config_dt)

 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'latCell',latCell)

 icsdlw_l => this%icsdlw_l
 icsdsw_l => this%icsdsw_l
 xlon_l   => this%xlon_l
 xlat_l   => this%xlat_l
 sinlat_l => this%sinlat_l
 coslat_l => this%coslat_l

!... pointers to model_parameters:
 me       = mdl_parameters%me
 lonr     = mdl_parameters%lonr
 latr     = mdl_parameters%latr
 idat(:)  = mdl_parameters%idat(:)

!... date conversion:
 xtime_h = xtime_s/3600._RKIND
 write(0,100) itimestep
 write(0,101) itimestep,year,month,day,hour,julday,curr_julday,utc_h,xtime_h
 write(0,*) '--- idat(1) = ', idat(1)
 write(0,*) '--- idat(2) = ', idat(2)
 write(0,*) '--- idat(3) = ', idat(3)
 write(0,*) '--- idat(4) = ', idat(4)
 write(0,*) '--- idat(5) = ', idat(5)
 write(0,*) '--- idat(6) = ', idat(6)
 write(0,*) '--- idat(7) = ', idat(7)
 write(0,*) '--- idat(8) = ', idat(8)
 write(0,*)

 jdate(:) = 0
 jdate(1) = year
 jdate(2) = month
 jdate(3) = day
 jdate(4) = utc_h
 jdate(5) = hour
 jdate(6) = minute
 jdate(7) = second
 jdate(8) = msecond

 fhour = utc_h
 solhr = mod(fhour+jdate(5),cons_24)

 write(0,*) '--- jdate(1) = ', jdate(1)
 write(0,*) '--- jdate(2) = ', jdate(2)
 write(0,*) '--- jdate(3) = ', jdate(3)
 write(0,*) '--- jdate(4) = ', jdate(4)
 write(0,*) '--- jdate(5) = ', jdate(5)
 write(0,*) '--- jdate(6) = ', jdate(6)
 write(0,*) '--- jdate(7) = ', jdate(7)
 write(0,*) '--- jdate(8) = ', jdate(8)
 write(0,*)
 write(0,*) '--- fhour    = ', fhour
 write(0,*) '--- solhr    = ', solhr
 write(0,*)

 lssav  = .true.
 lprnt  = .false.
 lslwr  = l_radtlw
 lsswr  = l_radtsw

 im     = nCellsSolve
 ix     = nCellsSolve
 kdt    = itimestep
 ipt    = nCellsSolve

 deltim = config_dt
 dtlw   = dt_radtlw
 dtsw   = dt_radtsw

 do iCell = 1, nCellsSolve
    xlon_l(iCell)   = lonCell(iCell)
    xlat_l(iCell)   = latCell(iCell)
    sinlat_l(iCell) = sin(latCell(iCell))
    coslat_l(iCell) = cos(latCell(iCell))
 enddo

!... for now:
 write(0,*) '--- isubc_lw  = ', isubc_lw
 write(0,*) '--- isubc_sw  = ', isubc_sw
 write(0,*) '--- ipsd0     = ', ipsd0

 if(.not.allocated(numrdm)) allocate(numrdm(2*nCellsSolve))
 if(.not.allocated(ixseed)) allocate(ixseed(1:nCellsSolve,1,2))

 ipseed = mod(nint(100.0*sqrt(fhour*3600)),ipsdlim) + 1 + ipsd0
 write(0,*) '--- ipseed    = ', ipseed
 call random_setseed(ipseed,stat)
 write(0,*) '--- end random_setseed:'
 call random_index(ipsdlim,numrdm,stat)
 write(0,*) '--- end random_index:'
 do n = 1, 2
    do iCell = 1, nCellsSolve
       ixseed(iCell,1,n) = numrdm(iCell+(latr-1)*lonr+(n-1)*latr)
    enddo
 enddo

 if(isubc_lw == 2 .or. isubc_sw == 2) then
    do iCell = 1, nCellsSolve
       icsdsw_l(iCell) = ixseed(iCell,1,1)
       icsdlw_l(iCell) = ixseed(iCell,1,2)
    enddo
 else
    do iCell = 1, nCellsSolve
       icsdlw_l(iCell) = 0
       icsdsw_l(iCell) = 0
    enddo
 endif

 if(allocated(numrdm)) deallocate(numrdm)
 if(allocated(ixseed)) deallocate(ixseed)

!... updates stuff:
 write(0,*) '--- enter subroutine radupdate:'
 call radupdate(idat,jdate,dtsw,deltim,lsswr,me,slag,sdec,cdec,solcon)
 write(0,*) '--- lsswr  = ', lsswr
 write(0,*) '--- slag   = ', slag
 write(0,*) '--- sdec   = ', sdec
 write(0,*) '--- cdec   = ', cdec
 write(0,*) '--- solcon = ', solcon
 write(0,*) '--- end subroutine radupdate'

!... fills variables and arrays contains in dynamic_parameters:
 write(0,*) '--- enter subroutine dyn_param_setrad:'
 call this%setrad( &
           xlon_l , xlat_l   , sinlat_l , coslat_l , solhr , ix    , im    , kdt   , jdate , &
           solcon , icsdsw_l , icsdlw_l , dtlw     , dtsw  , lsswr , lslwr , lssav , ipt   , &
           lprnt  , deltim   ,  slag    , sdec     , cdec)
 write(0,*) '--- end subroutine dyn_param_setrad' 

!... formats:
 100 format(' --- enter subroutine mpas_dyn_param_setrad: ',i6)
 101 format('     ITIMESTEP   =', i6   ,/, &
            '     YEAR        =', i6   ,/, &
            '     MONTH       =', i6   ,/, &
            '     DAY         =', i6   ,/, &
            '     HOUR        =', i6   ,/, &
            '     JULDAY      =', i6   ,/, &
            '     CURR_JULDAY =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     XTIME_H     =', f16.9,/)

 write(0,*) '--- end subroutine mpas_dyn_param_setrad'
 write(0,*)

 end subroutine mpas_dyn_param_setrad

!=================================================================================================================
 subroutine mpas_dyn_param_setphys(this,configs,mesh,itimestep,xtime_s)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!inout argumentsL
 class(mpas_dynamic_parameters),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve

 real(kind=RKIND),pointer:: config_dt
 real(kind=RKIND),dimension(:),pointer:: lonCell,latCell

 integer,dimension(:),pointer:: nlons_l
 real(kind=RKIND),dimension(:),pointer:: xlon_l,xlat_l,sinlat_l,coslat_l

!local variables:
 logical:: lscca,lssav,lsswr

 integer:: iCell,im,ix,kdt,lat,nnp,nsswr

 real(kind=RKIND),parameter:: tol = 0.01
 real(kind=RKIND):: clstp,dtf,dtp,fhour,fhswr,xtime_h
 real(kind=RKIND):: solhr
!real(kind=RKIND):: slag,sdec,cdec,solhr

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_config(configs,'config_dt'  ,config_dt  )

 call mpas_pool_get_array(mesh,'lonCell',lonCell)
 call mpas_pool_get_array(mesh,'latCell',latCell)

 xlon_l   => this%xlon_l
 xlat_l   => this%xlat_l
 sinlat_l => this%sinlat_l
 coslat_l => this%coslat_l
 nlons_l  => this%nlons_l

!... date conversion:
 xtime_h = xtime_s/3600._RKIND
 write(0,100) itimestep
 write(0,101) itimestep,year,month,day,hour,julday,curr_julday,utc_h,xtime_h

 fhour = utc_h
 solhr = mod(fhour+hour,cons_24)
 lsswr = l_radtsw

!... initialization of variable clstp needed in the calculation of bottom and top pressures for
!... convective clouds (Ken Campana diagnostics). set switch for saving convective clouds:
 fhswr = dt_radtsw / 3600._RKIND !frequency of shortwave radiation calls per hour.
 nsswr=nint(fhswr*3600._RKIND/config_dt)
 if(nsswr.le.0.or.abs(nsswr-fhswr*3600._RKIND/config_dt).gt.tol) &
    call physics_error_fatal('--- Error in mpas_dyn_param_setphys: nsswr < 0')

 lscca = mod(itimestep,nsswr) == 0
 if(lscca.and.l_radtsw) then
    clstp = 1100.+min(fhswr,fhour,cons_99)  !initialize,accumulate,convert
 elseif(lscca) then
    clstp = 0100.+min(fhswr,fhour,cons_99)  !accumulate,convert
 elseif(lsswr) then
    clstp = 1100.                           !initialize,accumulate
 else
    clstp = 0100.                           !accumulate
 endif
!write(0,*) '--- fhswr = ', fhswr
!write(0,*) '--- nsswr = ', nsswr
!write(0,*) '--- lscca = ', lscca
!write(0,*) '--- clstp = ', clstp

 do iCell = 1, nCellsSolve
    xlon_l(iCell)   = lonCell(iCell)
    xlat_l(iCell)   = latCell(iCell)
    sinlat_l(iCell) = sin(latCell(iCell))
    coslat_l(iCell) = cos(latCell(iCell))
    nlons_l(iCell)  = nCellsSolve
 enddo

 im    = nCellsSolve
 ix    = nCellsSolve
 nnp   = 1
 lat   = 1
 kdt   = itimestep
 dtf   = config_dt
 dtp   = config_dt
 lssav = .true.

!cdec = this%cdec
!sdec = this%sdec
!slag = this%slag

 write(0,*)
 write(0,*) '--- fhour  = ', fhour
 write(0,*) '--- solhr  = ', solhr
 write(0,*) '--- slag   = ', slag
 write(0,*) '--- sdec   = ', sdec
 write(0,*) '--- cdec   = ', cdec

!... fills variables and arrays contains in dynamic_parameters:
 write(0,*) '--- enter subroutine dyn_param_setphys:'
 call this%setphys (&
           xlon_l , xlat_l , sinlat_l , coslat_l , solhr , ix      , im    , kdt  , lssav , &
           lat    , dtp    , dtf      , clstp    , nnp   , nlons_l , fhour , slag , sdec  , &
           cdec)
 write(0,*) '--- end subroutine dyn_param_setphys'

!... formats:
 100 format(' --- enter subroutine mpas_dyn_param_setphys: ',i6)
 101 format('     ITIMESTEP   =', i6   ,/, &
            '     YEAR        =', i6   ,/, &
            '     MONTH       =', i6   ,/, &
            '     DAY         =', i6   ,/, &
            '     HOUR        =', i6   ,/, &
            '     JULDAY      =', i6   ,/, &
            '     CURR_JULDAY =', f16.9,/, &
            '     UTC_H       =', f16.9,/, &
            '     XTIME_H     =', f16.9,/)

 write(0,*) '--- end subroutine mpas_dyn_param_setphys'
 write(0,*)

 end subroutine mpas_dyn_param_setphys

!=================================================================================================================
 end module mpas_atmphys_nggps_dynamic_parameters
!=================================================================================================================
