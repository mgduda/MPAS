! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_update
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 implicit none
 private
 public:: physics_update, &
          physics_nggps_update_diagnostics

!modifications to sourcecode:
!----------------------------
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.

 contains
 

!=================================================================================================================
 subroutine physics_update(domain,dt)
!=================================================================================================================

!input arguments:
 type(domain_type),intent(in):: domain
 real(kind=R8KIND),intent(in) :: dt 

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!-----------------------------------------------------------------------------------------------------------------

!write(0,*)
!write(0,*) '--- begin physics_update NGGPS:'
!write(0,*) '--- end physics_update NGGPS.'

 end subroutine physics_update

!=================================================================================================================
 subroutine physics_nggps_update_diagnostics(mesh,nggps_diagnostics)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_diagnostics

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

 real(kind=R8KIND),dimension(:),pointer:: psmean,suntim,dlwsfc,ulwsfc,gflux,evbsa,evcwa,transa, &
                                        sbsnoa,snowca,snohfa,ep,dusfc,dvsfc,dtsfc,dqsfc,dugwd, &  
                                        dvgwd,cldwrk,cnvprcp,totprcp,runoff,srunoff,tmpmax,    &
                                        tmpmin,spfhmax,spfhmin

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------
!write(0,*) '--- enter subroutine physics_nggps_update_diagnostics:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_diagnostics,'psmean' ,psmean )
 call mpas_pool_get_array(nggps_diagnostics,'suntim' ,suntim )
 call mpas_pool_get_array(nggps_diagnostics,'dlwsfc' ,dlwsfc )
 call mpas_pool_get_array(nggps_diagnostics,'ulwsfc' ,ulwsfc ) 
 call mpas_pool_get_array(nggps_diagnostics,'gflux'  ,gflux  )
 call mpas_pool_get_array(nggps_diagnostics,'evbsa'  ,evbsa  )
 call mpas_pool_get_array(nggps_diagnostics,'evcwa'  ,evcwa  )
 call mpas_pool_get_array(nggps_diagnostics,'transa' ,transa )
 call mpas_pool_get_array(nggps_diagnostics,'sbsnoa' ,sbsnoa )
 call mpas_pool_get_array(nggps_diagnostics,'snowca' ,snowca )
 call mpas_pool_get_array(nggps_diagnostics,'snohfa' ,snohfa )
 call mpas_pool_get_array(nggps_diagnostics,'ep'     ,ep     )
 call mpas_pool_get_array(nggps_diagnostics,'dusfc'  ,dusfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dvsfc'  ,dvsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dtsfc'  ,dtsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dqsfc'  ,dqsfc  )
 call mpas_pool_get_array(nggps_diagnostics,'dugwd'  ,dugwd  )
 call mpas_pool_get_array(nggps_diagnostics,'dvgwd'  ,dvgwd  )
 call mpas_pool_get_array(nggps_diagnostics,'cldwrk' ,cldwrk )
 call mpas_pool_get_array(nggps_diagnostics,'cnvprcp',cnvprcp)
 call mpas_pool_get_array(nggps_diagnostics,'totprcp',totprcp)
 call mpas_pool_get_array(nggps_diagnostics,'runoff' ,runoff )
 call mpas_pool_get_array(nggps_diagnostics,'srunoff',srunoff)
 call mpas_pool_get_array(nggps_diagnostics,'tmpmax' ,tmpmax )
 call mpas_pool_get_array(nggps_diagnostics,'tmpmin' ,tmpmin )
 call mpas_pool_get_array(nggps_diagnostics,'spfhmax',spfhmax)
 call mpas_pool_get_array(nggps_diagnostics,'spfhmin',spfhmin)

 do iCell = 1, nCellsSolve
    psmean(iCell)  = 0._R8KIND
    suntim(iCell)  = 0._R8KIND
    dlwsfc(iCell)  = 0._R8KIND
    ulwsfc(iCell)  = 0._R8KIND
    gflux(iCell)   = 0._R8KIND
    evbsa(iCell)   = 0._R8KIND
    evcwa(iCell)   = 0._R8KIND
    transa(iCell)  = 0._R8KIND
    sbsnoa(iCell)  = 0._R8KIND
    snowca(iCell)  = 0._R8KIND
    snohfa(iCell)  = 0._R8KIND
    ep(iCell)      = 0._R8KIND
    dusfc(iCell)   = 0._R8KIND
    dvsfc(iCell)   = 0._R8KIND
    dtsfc(iCell)   = 0._R8KIND
    dqsfc(iCell)   = 0._R8KIND
    dugwd(iCell)   = 0._R8KIND
    dvgwd(iCell)   = 0._R8KIND
    cldwrk(iCell)  = 0._R8KIND
    cnvprcp(iCell) = 0._R8KIND
    totprcp(iCell) = 0._R8KIND
    runoff(iCell)  = 0._R8KIND
    srunoff(iCell) = 0._R8KIND
    tmpmax(iCell)  = 0._R8KIND
    tmpmin(iCell)  = 0._R8KIND
    spfhmax(iCell) = 0._R8KIND
    spfhmin(iCell) = 0._R8KIND
 enddo

!write(0,*) '--- end subroutine physics_nggps_update_diagnostics'
!write(0,*)

 end subroutine physics_nggps_update_diagnostics

!=================================================================================================================
 end module mpas_atmphys_update
!=================================================================================================================
