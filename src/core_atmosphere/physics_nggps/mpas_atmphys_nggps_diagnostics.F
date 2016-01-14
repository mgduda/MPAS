! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_diagnostics
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_vars
 use machine
 use nuopc_physics
 use module_radlw_parameters, only: topflw_type
 use module_radsw_parameters, only: topfsw_type

 implicit none
 private
 public:: mpas_diagnostics

 type,extends(diagnostics):: mpas_diagnostics
    type(topflw_type),dimension(:),pointer:: topflw_l => null()
    type(topfsw_type),dimension(:),pointer:: topfsw_l => null()

    real(kind=kind_phys),dimension(:,:),pointer:: fluxr_l  => null()
    real(kind=kind_phys),dimension(:,:),pointer:: dswcmp_l => null()
    real(kind=kind_phys),dimension(:,:),pointer:: uswcmp_l => null()

    contains
       procedure:: mpas_setphys   => mpas_diagnostics_setphys
       procedure:: mpas_setrad    => mpas_diagnostics_setrad
       procedure:: mpas_write     => mpas_diagnostics_write
       procedure:: allocate_mem   => mpas_diagnostics_allocate
       procedure:: deallocate_mem => mpas_diagnostics_deallocate 

 end type

!modifications to sourcecode:
!----------------------------
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.
 
 contains


!=================================================================================================================
 subroutine mpas_diagnostics_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_diagnostics),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nOznPlCoeff,nFxr

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_diagnostics_allocate:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)
 call mpas_pool_get_dimension(mesh,'nFxr'       ,nFxr       )

 if(.not.associated(this%topflw_l)) allocate(this%topflw_l(1:nCellsSolve)                    )
 if(.not.associated(this%topfsw_l)) allocate(this%topfsw_l(1:nCellsSolve)                    )
 if(.not.associated(this%fluxr_l) ) allocate(this%fluxr_l(1:nCellsSolve,nFxr)                )
 if(.not.associated(this%dswcmp_l)) allocate(this%dswcmp_l(1:nCellsSolve,4)                  )
 if(.not.associated(this%uswcmp_l)) allocate(this%uswcmp_l(1:nCellsSolve,4)                  ) 

 if(.not.allocated(dt3dt)         ) allocate(dt3dt(1:nCellsSolve,1:nVertLevels,6)            ) 
 if(.not.allocated(dq3dt)         ) allocate(dq3dt(1:nCellsSolve,1:nVertLevels,5+nOznPlCoeff)) 
 if(.not.allocated(du3dt)         ) allocate(du3dt(1:nCellsSolve,1:nVertLevels,4)            )
 if(.not.allocated(dv3dt)         ) allocate(dv3dt(1:nCellsSolve,1:nVertLevels,4)            )
 if(.not.allocated(dqdt_v)        ) allocate(dqdt_v(1:nCellsSolve,1:nVertLevels)             )

!write(0,*) '--- end subroutine mpas_diagnostics_allocate'

 end subroutine mpas_diagnostics_allocate

!=================================================================================================================
 subroutine mpas_diagnostics_deallocate(this)
!=================================================================================================================

 class(mpas_diagnostics),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_diagnostics_deallocate:'

 if(associated(this%topflw_l)) deallocate(this%topflw_l)
 if(associated(this%topfsw_l)) deallocate(this%topfsw_l)
 if(associated(this%fluxr_l) ) deallocate(this%fluxr_l )
 if(associated(this%dswcmp_l)) deallocate(this%dswcmp_l)
 if(associated(this%uswcmp_l)) deallocate(this%uswcmp_l) 

 if(allocated(dt3dt)         ) deallocate(dt3dt        ) 
 if(allocated(dq3dt)         ) deallocate(dq3dt        ) 
 if(allocated(du3dt)         ) deallocate(du3dt        )
 if(allocated(dv3dt)         ) deallocate(dv3dt        )
 if(allocated(dqdt_v)        ) deallocate(dqdt_v       )

!write(0,*) '--- end subroutine mpas_diagnostics_deallocate'

 end subroutine mpas_diagnostics_deallocate

!=================================================================================================================
 subroutine mpas_diagnostics_setphys(this,mesh,nggps_diagnostics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_diagnostics

!inout arguments:
 class(mpas_diagnostics),intent(inout):: this

!local pointers:
 integer, pointer:: nCellsSolve,nVertLevels
 real(kind=RKIND),dimension(:),pointer:: srunoff,evbsa,evcwa,snohfa,transa,sbsnoa,snowca,soilm,tmpmin,tmpmax,   &
                                         dusfc,dvsfc,dtsfc,dqsfc,totprcp,gflux,dlwsfc,ulwsfc,suntim,runoff,ep,  &
                                         cldwrk,dugwd,dvgwd,psmean,cnvprcp,spfhmin,spfhmax,rain,rainc,u10m,     &
                                         v10m,zlvl,psurf,hpbl,pwat,t1,q1,u1,v1,chh,cmm,dlwsfci,ulwsfci,dswsfci, &
                                         uswsfci,dusfci,dvsfci,dtsfci,dqsfci,gfluxi,epi,smcwlt2,smcref2,wet1,sr

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_diagnostics_setphys:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_diagnostics,'srunoff',srunoff)
 call mpas_pool_get_array(nggps_diagnostics,'evbsa'  ,evbsa  )
 call mpas_pool_get_array(nggps_diagnostics,'evcwa'  ,evcwa  )
 call mpas_pool_get_array(nggps_diagnostics,'snohfa' ,snohfa )
 call mpas_pool_get_array(nggps_diagnostics,'transa' ,transa )
 call mpas_pool_get_array(nggps_diagnostics,'sbsnoa' ,sbsnoa )
 call mpas_pool_get_array(nggps_diagnostics,'snowca' ,snowca )
 call mpas_pool_get_array(nggps_diagnostics,'soilm'  ,soilm  )
 call mpas_pool_get_array(nggps_diagnostics,'tmpmin' ,tmpmin )
 call mpas_pool_get_array(nggps_diagnostics,'tmpmax' ,tmpmax )
 call mpas_pool_get_array(nggps_diagnostics,'dusfc'  ,dusfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dvsfc'  ,dvsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dtsfc'  ,dtsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dqsfc'  ,dqsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'totprcp',totprcp)
 call mpas_pool_get_array(nggps_diagnostics,'gflux'  ,gflux  )
 call mpas_pool_get_array(nggps_diagnostics,'dlwsfc' ,dlwsfc )
 call mpas_pool_get_array(nggps_diagnostics,'ulwsfc' ,ulwsfc )
 call mpas_pool_get_array(nggps_diagnostics,'suntim' ,suntim )
 call mpas_pool_get_array(nggps_diagnostics,'runoff' ,runoff )
 call mpas_pool_get_array(nggps_diagnostics,'ep'     ,ep     )
 call mpas_pool_get_array(nggps_diagnostics,'cldwrk' ,cldwrk )
 call mpas_pool_get_array(nggps_diagnostics,'dugwd'  ,dugwd  )
 call mpas_pool_get_array(nggps_diagnostics,'dvgwd'  ,dvgwd  )
 call mpas_pool_get_array(nggps_diagnostics,'psmean' ,psmean )
 call mpas_pool_get_array(nggps_diagnostics,'cnvprcp',cnvprcp)
 call mpas_pool_get_array(nggps_diagnostics,'spfhmin',spfhmin)
 call mpas_pool_get_array(nggps_diagnostics,'spfhmax',spfhmax)
 call mpas_pool_get_array(nggps_diagnostics,'rain'   ,rain   )
 call mpas_pool_get_array(nggps_diagnostics,'rainc'  ,rainc  )
 call mpas_pool_get_array(nggps_diagnostics,'u10m'   ,u10m   )
 call mpas_pool_get_array(nggps_diagnostics,'v10m'   ,v10m   )
 call mpas_pool_get_array(nggps_diagnostics,'zlvl'   ,zlvl   )
 call mpas_pool_get_array(nggps_diagnostics,'psurf'  ,psurf  )
 call mpas_pool_get_array(nggps_diagnostics,'hpbl'   ,hpbl   )
 call mpas_pool_get_array(nggps_diagnostics,'pwat'   ,pwat   )
 call mpas_pool_get_array(nggps_diagnostics,'t1'     ,t1     )
 call mpas_pool_get_array(nggps_diagnostics,'q1'     ,q1     )
 call mpas_pool_get_array(nggps_diagnostics,'u1'     ,u1     )
 call mpas_pool_get_array(nggps_diagnostics,'v1'     ,v1     )
 call mpas_pool_get_array(nggps_diagnostics,'chh'    ,chh    )
 call mpas_pool_get_array(nggps_diagnostics,'cmm'    ,cmm    )
 call mpas_pool_get_array(nggps_diagnostics,'dlwsfci',dlwsfci)
 call mpas_pool_get_array(nggps_diagnostics,'ulwsfci',ulwsfci)
 call mpas_pool_get_array(nggps_diagnostics,'dswsfci',dswsfci)
 call mpas_pool_get_array(nggps_diagnostics,'uswsfci',uswsfci)
 call mpas_pool_get_array(nggps_diagnostics,'dusfci' ,dusfci )
 call mpas_pool_get_array(nggps_diagnostics,'dvsfci' ,dvsfci )
 call mpas_pool_get_array(nggps_diagnostics,'dtsfci' ,dtsfci )
 call mpas_pool_get_array(nggps_diagnostics,'dqsfci' ,dqsfci )
 call mpas_pool_get_array(nggps_diagnostics,'gfluxi' ,gfluxi )
 call mpas_pool_get_array(nggps_diagnostics,'epi'    ,epi    )
 call mpas_pool_get_array(nggps_diagnostics,'smcwlt2',smcwlt2)
 call mpas_pool_get_array(nggps_diagnostics,'smcref2',smcref2)
 call mpas_pool_get_array(nggps_diagnostics,'wet1'   ,wet1   )
 call mpas_pool_get_array(nggps_diagnostics,'sr'     ,sr     )

!initialization of local physics tendencies:
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       dqdt_v(iCell,k)  = 0._kind_phys
       dt3dt(iCell,k,:) = 0._kind_phys
       dq3dt(iCell,k,:) = 0._kind_phys
       du3dt(iCell,k,:) = 0._kind_phys
       dv3dt(iCell,k,:) = 0._kind_phys
    enddo
 enddo

 call this%setphys( &
           srunoff , evbsa   , evcwa   , snohfa  , transa  , sbsnoa  , &
           snowca  , soilm   , tmpmin  , tmpmax  , dusfc   , dvsfc   , &
           dtsfc   , dqsfc   , totprcp , gflux   , dlwsfc  , ulwsfc  , &
           suntim  , runoff  , ep      , cldwrk  , dugwd   , dvgwd   , &
           psmean  , cnvprcp , spfhmin , spfhmax , rain    , rainc   , &
           dt3dt   , dq3dt   , du3dt   , dv3dt   , dqdt_v  , u10m    , &
           v10m    , zlvl    , psurf   , hpbl    , pwat    , t1      , &
           q1      , u1      , v1      , chh     , cmm     , dlwsfci , &
           ulwsfci , dswsfci , uswsfci , dusfci  , dvsfci  , dtsfci  , &
           dqsfci  , gfluxi  , epi     , smcwlt2 , smcref2 , wet1    , &
           sr)

!write(0,*) '--- end subroutine mpas_diagnostics_setphys'

 end subroutine mpas_diagnostics_setphys

!=================================================================================================================
 subroutine mpas_diagnostics_setrad(this,mesh,nggps_diagnostics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_diagnostics

!inout arguments:
 class(mpas_diagnostics),intent(inout):: this

!local pointers:
 integer,pointer:: nFxr,nCellsSolve
 real(kind=RKIND),dimension(:),pointer  :: toplwup,toplwupc,topswup,topswupc,topswdn
 real(kind=RKIND),dimension(:,:),pointer:: fluxr,dswcmp,uswcmp

 type(topflw_type),dimension(:),pointer:: topflw_l
 type(topfsw_type),dimension(:),pointer:: topfsw_l
 real(kind=kind_phys),dimension(:,:),pointer:: fluxr_l,dswcmp_l,uswcmp_l

!local variables:
 integer:: iCell,n

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_diagnostics_setrad:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nFxr'       ,nFxr       )

 call mpas_pool_get_array(nggps_diagnostics,'toplwup' ,toplwup )
 call mpas_pool_get_array(nggps_diagnostics,'toplwupc',toplwupc)
 call mpas_pool_get_array(nggps_diagnostics,'topswup' ,topswup )
 call mpas_pool_get_array(nggps_diagnostics,'topswupc',topswupc)
 call mpas_pool_get_array(nggps_diagnostics,'topswdn' ,topswdn )
 call mpas_pool_get_array(nggps_diagnostics,'fluxr'   ,fluxr   )
 call mpas_pool_get_array(nggps_diagnostics,'dswcmp'  ,dswcmp  )
 call mpas_pool_get_array(nggps_diagnostics,'uswcmp'  ,uswcmp  )

 fluxr_l  => this%fluxr_l
 dswcmp_l => this%dswcmp_l
 uswcmp_l => this%uswcmp_l
 topflw_l => this%topflw_l
 topfsw_l => this%topfsw_l

 do iCell = 1, nCellsSolve
    topflw_l(iCell)%upfxc = toplwup(iCell)
    topflw_l(iCell)%upfx0 = toplwupc(iCell)
    topfsw_l(iCell)%upfxc = topswup(iCell)
    topfsw_l(iCell)%upfx0 = topswupc(iCell)
    topfsw_l(iCell)%dnfxc = topswdn(iCell)
    do n = 1, nFxr
       fluxr_l(iCell,n) = fluxr(n,iCell)
    enddo
    do n = 1, 4
       dswcmp_l(iCell,n) = dswcmp(n,iCell)
       uswcmp_l(iCell,n) = uswcmp(n,iCell)
    enddo
 enddo

 call this%setrad(nFxr,fluxr_l,topfsw_l,topflw_l,dswcmp_l,uswcmp_l)
       
!write(0,*) '--- end subroutine mpas_diagnostics_setrad'
!write(0,*)

 end subroutine mpas_diagnostics_setrad

!=================================================================================================================
 subroutine mpas_diagnostics_write(this,mesh,nggps_diagnostics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)   :: mesh
 class(mpas_diagnostics),intent(in):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_diagnostics

!local pointers:
 integer,pointer:: nCellsSolve
 real(kind=RKIND),dimension(:),pointer:: srunoff,evbsa,evcwa,snohfa,transa,sbsnoa,snowca,soilm,tmpmin,tmpmax,   &
                                         dusfc,dvsfc,dtsfc,dqsfc,totprcp,gflux,dlwsfc,ulwsfc,suntim,runoff,ep,  &
                                         cldwrk,dugwd,dvgwd,psmean,cnvprcp,spfhmin,spfhmax,rain,rainc,u10m,     &
                                         v10m,zlvl,psurf,hpbl,pwat,t1,q1,u1,v1,chh,cmm,dlwsfci,ulwsfci,dswsfci, &
                                         uswsfci,dusfci,dvsfci,dtsfci,dqsfci,gfluxi,epi,smcwlt2,smcref2,wet1,sr

!local variables:
 integer:: iCell,n

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_diagnostics_write:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(nggps_diagnostics,'srunoff'   ,srunoff   )
 call mpas_pool_get_array(nggps_diagnostics,'evbsa'     ,evbsa     )
 call mpas_pool_get_array(nggps_diagnostics,'evcwa'     ,evcwa     )
 call mpas_pool_get_array(nggps_diagnostics,'snohfa'    ,snohfa    )
 call mpas_pool_get_array(nggps_diagnostics,'transa'    ,transa    )
 call mpas_pool_get_array(nggps_diagnostics,'sbsnoa'    ,sbsnoa    )
 call mpas_pool_get_array(nggps_diagnostics,'snowca'    ,snowca    )
 call mpas_pool_get_array(nggps_diagnostics,'soilm'     ,soilm     )
 call mpas_pool_get_array(nggps_diagnostics,'tmpmin'    ,tmpmin    )
 call mpas_pool_get_array(nggps_diagnostics,'tmpmax'    ,tmpmax    )
 call mpas_pool_get_array(nggps_diagnostics,'dusfc'     ,dusfc     )
 call mpas_pool_get_array(nggps_diagnostics,'dvsfc'     ,dvsfc     )
 call mpas_pool_get_array(nggps_diagnostics,'dtsfc'     ,dtsfc     )
 call mpas_pool_get_array(nggps_diagnostics,'dqsfc'     ,dqsfc     )
 call mpas_pool_get_array(nggps_diagnostics,'totprcp'   ,totprcp   )
 call mpas_pool_get_array(nggps_diagnostics,'gflux'     ,gflux     )
 call mpas_pool_get_array(nggps_diagnostics,'dlwsfc'    ,dlwsfc    )
 call mpas_pool_get_array(nggps_diagnostics,'ulwsfc'    ,ulwsfc    )
 call mpas_pool_get_array(nggps_diagnostics,'suntim'    ,suntim    )
 call mpas_pool_get_array(nggps_diagnostics,'runoff'    ,runoff    )
 call mpas_pool_get_array(nggps_diagnostics,'ep'        ,ep        )
 call mpas_pool_get_array(nggps_diagnostics,'cldwrk'    ,cldwrk    )
 call mpas_pool_get_array(nggps_diagnostics,'dugwd'     ,dugwd     )
 call mpas_pool_get_array(nggps_diagnostics,'dvgwd'     ,dvgwd     )
 call mpas_pool_get_array(nggps_diagnostics,'psmean'    ,psmean    )
 call mpas_pool_get_array(nggps_diagnostics,'cnvprcp'   ,cnvprcp   )
 call mpas_pool_get_array(nggps_diagnostics,'spfhmin'   ,spfhmin   )
 call mpas_pool_get_array(nggps_diagnostics,'spfhmax'   ,spfhmax   )
 call mpas_pool_get_array(nggps_diagnostics,'rain'      ,rain      )
 call mpas_pool_get_array(nggps_diagnostics,'rainc'     ,rainc     )
 call mpas_pool_get_array(nggps_diagnostics,'u10m'      ,u10m      )
 call mpas_pool_get_array(nggps_diagnostics,'v10m'      ,v10m      )
 call mpas_pool_get_array(nggps_diagnostics,'zlvl'      ,zlvl      )
 call mpas_pool_get_array(nggps_diagnostics,'psurf'     ,psurf     )
 call mpas_pool_get_array(nggps_diagnostics,'hpbl'      ,hpbl      )
 call mpas_pool_get_array(nggps_diagnostics,'pwat'      ,pwat      )
 call mpas_pool_get_array(nggps_diagnostics,'t1'        ,t1        )
 call mpas_pool_get_array(nggps_diagnostics,'q1'        ,q1        )
 call mpas_pool_get_array(nggps_diagnostics,'u1'        ,u1        )
 call mpas_pool_get_array(nggps_diagnostics,'v1'        ,v1        )
 call mpas_pool_get_array(nggps_diagnostics,'chh'       ,chh       )
 call mpas_pool_get_array(nggps_diagnostics,'cmm'       ,cmm       )
 call mpas_pool_get_array(nggps_diagnostics,'dlwsfci'   ,dlwsfci   )
 call mpas_pool_get_array(nggps_diagnostics,'ulwsfci'   ,ulwsfci   )
 call mpas_pool_get_array(nggps_diagnostics,'dswsfci'   ,dswsfci   )
 call mpas_pool_get_array(nggps_diagnostics,'uswsfci'   ,uswsfci   )
 call mpas_pool_get_array(nggps_diagnostics,'dusfci'    ,dusfci    )
 call mpas_pool_get_array(nggps_diagnostics,'dvsfci'    ,dvsfci    )
 call mpas_pool_get_array(nggps_diagnostics,'dtsfci'    ,dtsfci    )
 call mpas_pool_get_array(nggps_diagnostics,'dqsfci'    ,dqsfci    )
 call mpas_pool_get_array(nggps_diagnostics,'gfluxi'    ,gfluxi    )
 call mpas_pool_get_array(nggps_diagnostics,'epi'       ,epi       )
 call mpas_pool_get_array(nggps_diagnostics,'smcwlt2'   ,smcwlt2   )
 call mpas_pool_get_array(nggps_diagnostics,'smcref2'   ,smcref2   )
 call mpas_pool_get_array(nggps_diagnostics,'wet1'      ,wet1      )
 call mpas_pool_get_array(nggps_diagnostics,'sr'        ,sr        )

!... back to Registry:
 srunoff => this%srunoff
 evbsa   => this%evbsa
 evcwa   => this%evcwa
 snohfa  => this%snohfa
 transa  => this%transa
 sbsnoa  => this%sbsnoa
 snowca  => this%snowca
 soilm   => this%soilm
 tmpmin  => this%tmpmin
 tmpmax  => this%tmpmax
 dusfc   => this%dusfc
 dvsfc   => this%dvsfc
 dtsfc   => this%dtsfc
 dqsfc   => this%dqsfc
 totprcp => this%totprcp
 gflux   => this%gflux
 dlwsfc  => this%dlwsfc
 ulwsfc  => this%ulwsfc
 suntim  => this%suntim
 runoff  => this%runoff
 ep      => this%ep
 cldwrk  => this%cldwrk
 dugwd   => this%dugwd
 dvgwd   => this%dvgwd
 psmean  => this%psmean
 cnvprcp => this%cnvprcp
 spfhmin => this%spfhmin
 spfhmax => this%spfhmax
 rain    => this%rain
 rainc   => this%rainc
 u10m    => this%u10m
 v10m    => this%v10m
 zlvl    => this%zlvl
 psurf   => this%psurf
 hpbl    => this%hpbl
 pwat    => this%pwat
 t1      => this%t1
 q1      => this%q1
 u1      => this%u1
 v1      => this%v1
 chh     => this%chh
 cmm     => this%cmm
 dlwsfci => this%dlwsfci
 ulwsfci => this%ulwsfci
 dswsfci => this%dswsfci
 uswsfci => this%uswsfci
 dusfci  => this%dusfci
 dvsfci  => this%dvsfci
 dtsfci  => this%dtsfci
 dqsfci  => this%dqsfci
 gfluxi  => this%gfluxi
 epi     => this%epi
 smcwlt2 => this%smcwlt2
 smcref2 => this%smcref2
 wet1    => this%wet1
 sr      => this%sr

!write(0,*) '--- end subroutine mpas_diagnostics_write'

 end subroutine mpas_diagnostics_write

!=================================================================================================================
 end module mpas_atmphys_nggps_diagnostics
!=================================================================================================================
