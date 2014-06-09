! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_update
 use mpas_grid_types

 implicit none
 private
 public:: physics_update, &
          physics_gfs_update_diagnostics


!>\brief update diagnostics.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-01-28.
!>
!>\details
!>
!> subroutines in mpas_atmphys_update:
!> -----------------------------------
!> physics_update                : not used.
!> physics_gfs_update_diagnostics: reset all accumulated diagnostics to zero.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * added initialization of variable dkh (accumulated vertical diffusion coefficient).
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains
 
!=================================================================================================================
 subroutine physics_update(domain,dt)
!=================================================================================================================

!input arguments:
 type(domain_type),intent(in):: domain
 real(kind=RKIND),intent(in) :: dt 

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!-----------------------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- begin physics_update:'
 write(0,*) '--- end physics_update.'

 end subroutine physics_update

!=================================================================================================================
 subroutine physics_gfs_update_diagnostics(mesh,diag_physics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout)::diag_physics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

 real(kind=RKIND),dimension(:),pointer  :: psmean,suntim,dlwsfc,ulwsfc,gflux,evbsa,evcwa,transa,    &
                                           sbsnoa,snowca,snohfa,ep,gtmp2m,gu10m,gv10m,gustar,gzorl, &
                                           goro,dusfc,dvsfc,dtsfc,dqsfc,gpblh,dugwd,dvgwd,cldwrk,   &
                                           cnvprcp,totprcp,runoff,srunoff,gsoil,tmpmax,tmpmin,      &
                                           spfhmax,spfhmin
 real(kind=RKIND),dimension(:,:),pointer:: dkh

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------
 write(0,*) '--- enter subroutine physics_gfs_update_diagnostics:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(diag_physics,'psmean' ,psmean )
 call mpas_pool_get_array(diag_physics,'suntim' ,suntim )
 call mpas_pool_get_array(diag_physics,'dlwsfc' ,dlwsfc )
 call mpas_pool_get_array(diag_physics,'ulwsfc' ,ulwsfc )
 call mpas_pool_get_array(diag_physics,'gflux'  ,gflux  )
 call mpas_pool_get_array(diag_physics,'evbsa'  ,evbsa  )
 call mpas_pool_get_array(diag_physics,'evcwa'  ,evcwa  )
 call mpas_pool_get_array(diag_physics,'transa' ,transa )
 call mpas_pool_get_array(diag_physics,'sbsnoa' ,sbsnoa )
 call mpas_pool_get_array(diag_physics,'snowca' ,snowca )
 call mpas_pool_get_array(diag_physics,'snohfa' ,snohfa )
 call mpas_pool_get_array(diag_physics,'ep'     ,ep     )
 call mpas_pool_get_array(diag_physics,'gtmp2m' ,gtmp2m )
 call mpas_pool_get_array(diag_physics,'gu10m'  ,gu10m  )
 call mpas_pool_get_array(diag_physics,'gv10m'  ,gv10m  )
 call mpas_pool_get_array(diag_physics,'gustar' ,gustar )
 call mpas_pool_get_array(diag_physics,'gzorl'  ,gzorl  )
 call mpas_pool_get_array(diag_physics,'goro'   ,goro   )
 call mpas_pool_get_array(diag_physics,'dusfc'  ,dusfc  )
 call mpas_pool_get_array(diag_physics,'dvsfc'  ,dvsfc  )
 call mpas_pool_get_array(diag_physics,'dtsfc'  ,dtsfc  )
 call mpas_pool_get_array(diag_physics,'dqsfc'  ,dqsfc  )
 call mpas_pool_get_array(diag_physics,'gpblh'  ,gpblh  )
 call mpas_pool_get_array(diag_physics,'dugwd'  ,dugwd  )
 call mpas_pool_get_array(diag_physics,'dvgwd'  ,dvgwd  )
 call mpas_pool_get_array(diag_physics,'cldwrk' ,cldwrk )
 call mpas_pool_get_array(diag_physics,'cnvprcp',cnvprcp)
 call mpas_pool_get_array(diag_physics,'totprcp',totprcp)
 call mpas_pool_get_array(diag_physics,'runoff' ,runoff )
 call mpas_pool_get_array(diag_physics,'srunoff',srunoff)
 call mpas_pool_get_array(diag_physics,'gsoil'  ,gsoil  )
 call mpas_pool_get_array(diag_physics,'tmpmax' ,tmpmax )
 call mpas_pool_get_array(diag_physics,'tmpmin' ,tmpmin )
 call mpas_pool_get_array(diag_physics,'spfhmax',spfhmax)
 call mpas_pool_get_array(diag_physics,'spfhmin',spfhmin)
 call mpas_pool_get_array(diag_physics,'dkh'    ,dkh    )

 do iCell = 1, nCellsSolve
    psmean(iCell)  = 0._RKIND
    suntim(iCell)  = 0._RKIND
    dlwsfc(iCell)  = 0._RKIND
    ulwsfc(iCell)  = 0._RKIND
    gflux(iCell)   = 0._RKIND
    evbsa(iCell)   = 0._RKIND
    evcwa(iCell)   = 0._RKIND
    transa(iCell)  = 0._RKIND
    sbsnoa(iCell)  = 0._RKIND
    snowca(iCell)  = 0._RKIND
    snohfa(iCell)  = 0._RKIND
    ep(iCell)      = 0._RKIND
    gtmp2m(iCell)  = 0._RKIND
    gu10m(iCell)   = 0._RKIND
    gv10m(iCell)   = 0._RKIND
    gustar(iCell)  = 0._RKIND
    gzorl(iCell)   = 0._RKIND
    goro(iCell)    = 0._RKIND
    dusfc(iCell)   = 0._RKIND
    dvsfc(iCell)   = 0._RKIND
    dtsfc(iCell)   = 0._RKIND
    dqsfc(iCell)   = 0._RKIND
    gpblh(iCell)   = 0._RKIND
    dugwd(iCell)   = 0._RKIND
    dvgwd(iCell)   = 0._RKIND
    cldwrk(iCell)  = 0._RKIND
    cnvprcp(iCell) = 0._RKIND
    totprcp(iCell) = 0._RKIND
    runoff(iCell)  = 0._RKIND
    srunoff(iCell) = 0._RKIND
    gsoil(iCell)   = 0._RKIND
    tmpmax(iCell)  = 0._RKIND
    tmpmin(iCell)  = 0._RKIND
    spfhmax(iCell) = 0._RKIND
    spfhmin(iCell) = 0._RKIND

    do k = 1, nVertLevels
       dkh(k,iCell) = 0._RKIND
    enddo

 enddo

 write(0,*) '--- end subroutine physics_gfs_update_diagnostics:'

 end subroutine physics_gfs_update_diagnostics

!=================================================================================================================
 end module mpas_atmphys_update
!=================================================================================================================
