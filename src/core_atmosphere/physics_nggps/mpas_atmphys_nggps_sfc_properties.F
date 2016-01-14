! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_sfc_properties
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_vars
 use machine
 use nuopc_physics

 implicit none
 private
 public:: mpas_sfc_properties

 type,extends(sfc_properties):: mpas_sfc_properties

    real(kind=kind_phys),dimension(:),pointer  :: rsltyp => null()
    real(kind=kind_phys),dimension(:),pointer  :: rvgtyp => null()
    real(kind=kind_phys),dimension(:,:),pointer:: hprime => null()

    contains
       procedure:: mpas_setphys   => mpas_sfc_prop_setphys
       procedure:: mpas_setrad    => mpas_sfc_prop_setrad
       procedure:: mpas_write     => mpas_sfc_prop_write
       procedure:: allocate_mem   => mpas_sfc_prop_allocate
       procedure:: deallocate_mem => mpas_sfc_prop_deallocate 

 end type

!modifications to sourcecode:
!----------------------------
! * revised the input of the static fields for the parameterization of the gravity wave drag over orography.
!   Laura D. Fowler (laura@ucar.edu) / 2015-12-09.
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.

 contains


!=================================================================================================================
 subroutine mpas_sfc_prop_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_sfc_properties),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve
 integer,pointer:: nMtvr

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_sfc_prop_allocate:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nMtvr'      ,nMtvr)

 if(.not.associated(this%rvgtyp)) allocate(this%rvgtyp(1:nCellsSolve)        )
 if(.not.associated(this%rsltyp)) allocate(this%rsltyp(1:nCellsSolve)        )
 if(.not.associated(this%hprim) ) allocate(this%hprim(1:nCellsSolve)         )
 if(.not.associated(this%hprime)) allocate(this%hprime(1:nCellsSolve,1:nMtvr))

!write(0,*) '--- end subroutine mpas_sfc_prop_allocate'

 end subroutine mpas_sfc_prop_allocate

!=================================================================================================================
 subroutine mpas_sfc_prop_deallocate(this)
!=================================================================================================================

 class(mpas_sfc_properties),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_sfc_prop_deallocate:'

 if(associated(this%rvgtyp)) deallocate(this%rvgtyp)
 if(associated(this%rsltyp)) deallocate(this%rsltyp)
 if(associated(this%hprim) ) deallocate(this%hprim )
 if(associated(this%hprime)) deallocate(this%hprime)

!write(0,*) '--- end subroutine mpas_surface_deallocate'

 end subroutine mpas_sfc_prop_deallocate

!=================================================================================================================
 subroutine mpas_sfc_prop_setphys(this,mesh,nggps_sfc_properties)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_sfc_properties

!inout arguments:
 class(mpas_sfc_properties),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nMtvr
 integer,dimension(:),pointer:: landmask,ivgtyp,isltyp

 real(kind=RKIND),dimension(:),pointer:: slope,shdmin,shdmax,snoalb,tg3,slmsk,vfrac,uustar,oro,oro_uf,   &
                                         hice,fice,tisfc,tsfc,snwdph,weasd,sncovr,zorl,canopy,ffmm,ffhh, &
                                         f10m,t2m,q2m
 real(kind=RKIND),dimension(:),pointer:: var2d,con,oa1,oa2,oa3,oa4,ol1,ol2,ol3,ol4,theoro, &
                                         gamoro,sigoro,elvoro

!local nggps pointers:
 real(kind=kind_phys),dimension(:),pointer  :: rsltyp,rvgtyp
 real(kind=kind_phys),dimension(:,:),pointer:: hprime

!local variables:
 integer:: iCell,n

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_sfc_prop_setphys:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nMtvr'      ,nMtvr      )

 call mpas_pool_get_array(nggps_sfc_properties,'landmask',landmask)
 call mpas_pool_get_array(nggps_sfc_properties,'ivgtyp'  ,ivgtyp  )
 call mpas_pool_get_array(nggps_sfc_properties,'isltyp'  ,isltyp  )
 call mpas_pool_get_array(nggps_sfc_properties,'slope'   ,slope   )
 call mpas_pool_get_array(nggps_sfc_properties,'shdmin'  ,shdmin  )
 call mpas_pool_get_array(nggps_sfc_properties,'shdmax'  ,shdmax  )
 call mpas_pool_get_array(nggps_sfc_properties,'snoalb'  ,snoalb  )
 call mpas_pool_get_array(nggps_sfc_properties,'tg3'     ,tg3     )
 call mpas_pool_get_array(nggps_sfc_properties,'slmsk'   ,slmsk   )
 call mpas_pool_get_array(nggps_sfc_properties,'vfrac'   ,vfrac   )
 call mpas_pool_get_array(nggps_sfc_properties,'uustar'  ,uustar  )
 call mpas_pool_get_array(nggps_sfc_properties,'oro'     ,oro     )
 call mpas_pool_get_array(nggps_sfc_properties,'oro_uf'  ,oro_uf  )
 call mpas_pool_get_array(nggps_sfc_properties,'hice'    ,hice    )
 call mpas_pool_get_array(nggps_sfc_properties,'fice'    ,fice    )
 call mpas_pool_get_array(nggps_sfc_properties,'tisfc'   ,tisfc   )
 call mpas_pool_get_array(nggps_sfc_properties,'tsfc'    ,tsfc    )
 call mpas_pool_get_array(nggps_sfc_properties,'snwdph'  ,snwdph  )
 call mpas_pool_get_array(nggps_sfc_properties,'weasd'   ,weasd   )
 call mpas_pool_get_array(nggps_sfc_properties,'sncovr'  ,sncovr  )
 call mpas_pool_get_array(nggps_sfc_properties,'zorl'    ,zorl    )
 call mpas_pool_get_array(nggps_sfc_properties,'canopy'  ,canopy  )
 call mpas_pool_get_array(nggps_sfc_properties,'ffmm'    ,ffmm    )
 call mpas_pool_get_array(nggps_sfc_properties,'ffhh'    ,ffhh    )
 call mpas_pool_get_array(nggps_sfc_properties,'f10m'    ,f10m    )
 call mpas_pool_get_array(nggps_sfc_properties,'t2m'     ,t2m     )
 call mpas_pool_get_array(nggps_sfc_properties,'q2m'     ,q2m     )

 call mpas_pool_get_array(nggps_sfc_properties,'var2d'   ,var2d   )
 call mpas_pool_get_array(nggps_sfc_properties,'con'     ,con     )
 call mpas_pool_get_array(nggps_sfc_properties,'oa1'     ,oa1     )
 call mpas_pool_get_array(nggps_sfc_properties,'oa2'     ,oa2     )
 call mpas_pool_get_array(nggps_sfc_properties,'oa3'     ,oa3     )
 call mpas_pool_get_array(nggps_sfc_properties,'oa4'     ,oa4     )
 call mpas_pool_get_array(nggps_sfc_properties,'ol1'     ,ol1     )
 call mpas_pool_get_array(nggps_sfc_properties,'ol2'     ,ol2     )
 call mpas_pool_get_array(nggps_sfc_properties,'ol3'     ,ol3     )
 call mpas_pool_get_array(nggps_sfc_properties,'ol4'     ,ol4     )
 call mpas_pool_get_array(nggps_sfc_properties,'theoro'  ,theoro  )
 call mpas_pool_get_array(nggps_sfc_properties,'gamoro'  ,gamoro  )
 call mpas_pool_get_array(nggps_sfc_properties,'sigoro'  ,sigoro  )
 call mpas_pool_get_array(nggps_sfc_properties,'elvoro'  ,elvoro  )

!... nggps arrays:
 rsltyp => this%rsltyp
 rvgtyp => this%rvgtyp
 hprime => this%hprime

 do iCell = 1, nCellsSolve
    slmsk(iCell) = real(landmask(iCell),kind_phys)
    if(slmsk(iCell) == 0._kind_phys .and. fice(iCell) .gt. 0.) slmsk(iCell) = 2._kind_phys

    slope(iCell)  = min(slope(iCell),9._RKIND)
    rsltyp(iCell) = real(isltyp(iCell),kind_phys)
    rvgtyp(iCell) = real(ivgtyp(iCell),kind_phys)
 enddo

!... fills the array hprime needed for the gravity wave drag over orography:
 if(nMtvr .ge. 6) then
    do iCell = 1, nCellsSolve
       hprime(iCell,1)  = var2d(iCell)
       hprime(iCell,2)  = con(iCell)
       hprime(iCell,3)  = oa1(iCell)
       hprime(iCell,4)  = oa2(iCell)
       hprime(iCell,5)  = oa3(iCell)
       hprime(iCell,6)  = oa4(iCell)

       if(nMtvr .ge. 10) then
          hprime(iCell,7)  = ol1(iCell)
          hprime(iCell,8)  = ol2(iCell)
          hprime(iCell,9)  = ol3(iCell)
          hprime(iCell,10) = ol4(iCell)
          if(nMtvr .eq. 14) then
             hprime(iCell,11) = theoro(iCell)
             hprime(iCell,12) = gamoro(iCell)
             hprime(iCell,13) = sigoro(iCell)
             hprime(iCell,14) = elvoro(iCell)
          endif
       endif
    enddo
 else
    do iCell = 1, nCellsSolve
       hprime(iCell,:) = 0._kind_phys
    enddo
 endif

!write(0,*) '--- before call to subroutine setphys: nMtvr = ', nMtvr
!do iCell = 1, 25
!   write(0,203) iCell,(hprime(iCell,n),n=1,10)
!enddo
!write(0,*)
!do iCell = 1, nCellsSolve
!   if(slmsk(iCell) .eq. 1._kind_phys) write(0,203) iCell,(hprime(iCell,n),n=11,14)
!enddo
!stop

 call this%setphys( &
           hprime , slope  , shdmin , shdmax , snoalb , tg3    , slmsk , vfrac ,  &
           rvgtyp , rsltyp , uustar , oro    , oro_uf , hice   , fice  , tisfc ,  &
           tsfc   , snwdph , weasd  , sncovr , zorl   , canopy , ffmm  , ffhh  ,  &
           f10m   , t2m    , q2m)

!do iCell = 1, nCellsSolve
!   if(slmsk(iCell) == 1._kind_phys) then
!      write(0,203) iCell,this%slope(iCell),this%shdmin(iCell),this%shdmax(iCell), this%snoalb(iCell),&
!         this%tg3(iCell),this%slmsk(iCell),this%vfrac(iCell),this%vtype(iCell),this%stype(iCell)
!   endif
!   if(this%fice(iCell) > 0._kind_phys) &
!      write(0,203) iCell,this%fice(iCell),this%tsfc(iCell),this%tisfc(iCell)
!   if(slmsk(iCell) == 1._kind_phys) &
!      write(0,203) iCell,this%oro(iCell),this%oro_uf(iCell)
!   if(weasd(iCell) > 0._kind_phys) &
!      write(0,203) iCell,this%sncovr(iCell),this%weasd(iCell),this%snowd(iCell)
!enddo

!formats:
 203 format(i6,10(1x,e15.8))

!write(0,*) '--- end subroutine mpas_sfc_prop_setphys'

 end subroutine mpas_sfc_prop_setphys

!=================================================================================================================
 subroutine mpas_sfc_prop_setrad(this,mesh,nggps_sfc_properties)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_sfc_properties

!inout arguments:
 class(mpas_sfc_properties),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:),pointer:: snoalb,slmsk,fice,tisfc,tsfc,snwdph,sncovr,zorl,facsf,facwf, &
                                         alnsf,alnwf,alvsf,alvwf
 real(kind=RKIND),dimension(:),pointer:: var2d

!local nggps pointers:
 real(kind=kind_phys),dimension(:),pointer:: hprim

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_sfc_prop_setrad:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(nggps_sfc_properties,'landmask',landmask)
 call mpas_pool_get_array(nggps_sfc_properties,'snoalb'  ,snoalb  )
 call mpas_pool_get_array(nggps_sfc_properties,'slmsk'   ,slmsk   )
 call mpas_pool_get_array(nggps_sfc_properties,'fice'    ,fice    )
 call mpas_pool_get_array(nggps_sfc_properties,'snwdph'  ,snwdph  )
 call mpas_pool_get_array(nggps_sfc_properties,'sncovr'  ,sncovr  )
 call mpas_pool_get_array(nggps_sfc_properties,'tsfc'    ,tsfc    )
 call mpas_pool_get_array(nggps_sfc_properties,'zorl'    ,zorl    )
 call mpas_pool_get_array(nggps_sfc_properties,'tisfc'   ,tisfc   )
 call mpas_pool_get_array(nggps_sfc_properties,'alnsf'   ,alnsf   )
 call mpas_pool_get_array(nggps_sfc_properties,'alnwf'   ,alnwf   )
 call mpas_pool_get_array(nggps_sfc_properties,'alvsf'   ,alvsf   )
 call mpas_pool_get_array(nggps_sfc_properties,'alvwf'   ,alvwf   )
 call mpas_pool_get_array(nggps_sfc_properties,'facsf'   ,facsf   )
 call mpas_pool_get_array(nggps_sfc_properties,'facwf '  ,facwf   )
 call mpas_pool_get_array(nggps_sfc_properties,'var2d'  , var2d   )

 hprim => this%hprim

 do iCell = 1, nCellsSolve
    slmsk(iCell) = real(landmask(iCell),kind_phys)
    if(slmsk(iCell) == 0._kind_phys .and. fice(iCell) .gt. 0.) slmsk(iCell) = 2._kind_phys

    !... standard deviation of orography:
    hprim(iCell) = var2d(iCell)
 enddo
 203 format(i6,10(1x,e15.8))

 call this%setrad( &
           slmsk , tsfc  , snwdph , sncovr , snoalb , zorl  , hprim , fice , &
           tisfc , alvsf , alnsf  , alvwf  , alnwf  , facsf , facwf )

!write(0,*) '--- end subroutine mpas_sfc_prop_setrad'

 end subroutine mpas_sfc_prop_setrad

!=================================================================================================================
 subroutine mpas_sfc_prop_write(this,nggps_sfc_properties)
!=================================================================================================================

!input arguments:
 class(mpas_sfc_properties),intent(in):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_sfc_properties

!local variables:
 real(kind=RKIND),dimension(:),pointer  :: alnsf,alvsf,alnwf,alvwf,slmsk,vfrac,uustar,hice,fice,tisfc,tsfc, &
                                           snwdph,weasd,sncovr,zorl,canopy,ffmm,ffhh,f10m,t2m,q2m

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_sfc_prop_write:'

 call mpas_pool_get_array(nggps_sfc_properties,'slmsk'   ,slmsk   )
 call mpas_pool_get_array(nggps_sfc_properties,'vfrac'   ,vfrac   )
 call mpas_pool_get_array(nggps_sfc_properties,'uustar'  ,uustar  )
 call mpas_pool_get_array(nggps_sfc_properties,'hice'    ,hice    )
 call mpas_pool_get_array(nggps_sfc_properties,'fice'    ,fice    )
 call mpas_pool_get_array(nggps_sfc_properties,'tisfc'   ,tisfc   )
 call mpas_pool_get_array(nggps_sfc_properties,'tsfc'    ,tsfc    )
 call mpas_pool_get_array(nggps_sfc_properties,'snwdph'  ,snwdph  )
 call mpas_pool_get_array(nggps_sfc_properties,'weasd'   ,weasd   )
 call mpas_pool_get_array(nggps_sfc_properties,'sncovr'  ,sncovr  )
 call mpas_pool_get_array(nggps_sfc_properties,'zorl'    ,zorl    )
 call mpas_pool_get_array(nggps_sfc_properties,'canopy'  ,canopy  )
 call mpas_pool_get_array(nggps_sfc_properties,'ffmm'    ,ffmm    )
 call mpas_pool_get_array(nggps_sfc_properties,'ffhh'    ,ffhh    )
 call mpas_pool_get_array(nggps_sfc_properties,'f10m'    ,f10m    )
 call mpas_pool_get_array(nggps_sfc_properties,'t2m'     ,t2m     )
 call mpas_pool_get_array(nggps_sfc_properties,'q2m'     ,q2m     )
 call mpas_pool_get_array(nggps_sfc_properties,'alnsf'   ,alnsf   )
 call mpas_pool_get_array(nggps_sfc_properties,'alvsf'   ,alvsf   )
 call mpas_pool_get_array(nggps_sfc_properties,'alnwf'   ,alnwf   )
 call mpas_pool_get_array(nggps_sfc_properties,'alvwf'   ,alvwf   )

!... back to Registry:
 slmsk  => this%slmsk
 vfrac  => this%vfrac
 uustar => this%uustar
 hice   => this%hice
 fice   => this%fice
 tisfc  => this%tisfc
 tsfc   => this%tsfc
 snwdph => this%snowd
 weasd  => this%weasd
 sncovr => this%sncovr
 zorl   => this%zorl
 canopy => this%canopy
 ffmm   => this%ffmm
 ffhh   => this%ffhh
 f10m   => this%f10m
 t2m    => this%t2m
 q2m    => this%q2m
 alnsf  => this%alnsf
 alvsf  => this%alvsf
 alnwf  => this%alnwf
 alvwf  => this%alvwf

!write(0,*) '--- end subroutine mpas_sfc_prop_write'

 end subroutine mpas_sfc_prop_write

!=================================================================================================================
 end module mpas_atmphys_nggps_sfc_properties
!=================================================================================================================
