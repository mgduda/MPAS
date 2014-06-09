! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_todynamics
 use mpas_dmpar
 use mpas_grid_types

 use mpas_atmphys_gfs_constants, only: R_d,R_v,degrad

 implicit none
 private
 public:: physics_addtend


!>\brief interface between the physics parameterizations and the non-hydrostatic dynamical core.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-05-01.
!>
!>\details
!>
!> subroutines in mpas_atmphys_todynamics:
!> ---------------------------------------
!> physics_addtend: add and mass-weigh tendencies before being added to dynamics tendencies.
!> tend_toEdges   : interpolate wind-tendencies from centers to edges of grid-cells.
!>
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains
 
!=================================================================================================================
 subroutine physics_addtend(block,mesh,state,diag,tend,tend_physics,configs,rk_step)
!=================================================================================================================

!input variables:
 type(block_type),pointer:: block
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(in):: configs
 integer, intent(in):: rk_step

!inout variables:
 type(mpas_pool_type),intent(inout):: tend

!local pointers:
 integer,pointer:: nCellsSolve,nEdges,nEdgesSolve,nVertLevels
 integer,pointer:: index_qv,index_qc,index_qo3

!local variables:
 integer:: iCell,iEdge,k,n

 real(kind=RKIND),dimension(:,:),pointer:: mass          ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:),pointer:: mass_edge     ! diag rho_edge
 real(kind=RKIND),dimension(:,:),pointer:: theta_m,qv
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND),dimension(:,:),pointer:: tend_thdiffusion,tend_qvdiffusion,  &
                                           tend_qcdiffusion,tend_qo3diffusion, &
                                           tend_udiffusion ,tend_vdiffusion 
 real(kind=RKIND),dimension(:,:),pointer:: tend_thcudeep,tend_qvcudeep,        &
                                           tend_qccudeep,tend_ucudeep,         &
                                           tend_vcudeep
 real(kind=RKIND),dimension(:,:),pointer:: tend_thcushallow,tend_qvcushallow,  &
                                           tend_qccushallow
 real(kind=RKIND),dimension(:,:),pointer:: tend_thradlw,tend_thradsw
 real(kind=RKIND),dimension(:,:),pointer:: tend_ugwdo,tend_vgwdo
 real(kind=RKIND),dimension(:,:),pointer:: tend_ugwdc,tend_vgwdc

 real(kind=RKIND),dimension(:,:),pointer:: tend_thphysics,tend_qvphysics,  &
                                           tend_qcphysics,tend_qo3physics, &
                                           tend_uphysics ,tend_vphysics

 real(kind=RKIND),dimension(:,:),pointer  :: tend_theta,tend_u
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars
 real(kind=RKIND),dimension(:,:),pointer:: tend_qv,tend_qc,tend_qo3

 real(kind=RKIND),dimension(:,:),allocatable:: tend_uphysics_Edge
 real(kind=RKIND),dimension(:,:),allocatable:: theta,tend_th

!=================================================================================================================

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nEdges'     ,nEdges     )
 call mpas_pool_get_dimension(mesh,'nEdgesSolve',nEdgesSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(state,'theta_m' ,theta_m,1)
 call mpas_pool_get_array(state,'scalars' ,scalars,1)
 call mpas_pool_get_array(state,'rho_zz'  ,mass,2   )
 call mpas_pool_get_array(diag ,'rho_edge',mass_edge)

 call mpas_pool_get_dimension(state,'index_qv' ,index_qv )
 call mpas_pool_get_dimension(state,'index_qc' ,index_qc )
 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)
 qv => scalars(index_qv,:,:)

!PBL and diffusion:
 call mpas_pool_get_array(tend_physics,'udiffusion'  ,tend_udiffusion  )
 call mpas_pool_get_array(tend_physics,'vdiffusion'  ,tend_vdiffusion  )
 call mpas_pool_get_array(tend_physics,'thdiffusion' ,tend_thdiffusion )
 call mpas_pool_get_array(tend_physics,'qvdiffusion' ,tend_qvdiffusion )
 call mpas_pool_get_array(tend_physics,'qcdiffusion' ,tend_qcdiffusion )
 call mpas_pool_get_array(tend_physics,'qo3diffusion',tend_qo3diffusion)

!deep convection and moist convective adjustment:
 call mpas_pool_get_array(tend_physics,'ucudeep' ,tend_ucudeep )
 call mpas_pool_get_array(tend_physics,'vcudeep' ,tend_vcudeep )
 call mpas_pool_get_array(tend_physics,'thcudeep',tend_thcudeep)
 call mpas_pool_get_array(tend_physics,'qvcudeep',tend_qvcudeep)
 call mpas_pool_get_array(tend_physics,'qccudeep',tend_qccudeep)

!shallow convection:
 call mpas_pool_get_array(tend_physics,'thcushallow',tend_thcushallow)
 call mpas_pool_get_array(tend_physics,'qvcushallow',tend_qvcushallow)
 call mpas_pool_get_array(tend_physics,'qccushallow',tend_qccushallow)

!long-wave and short-wave convection:
 call mpas_pool_get_array(tend_physics,'thradlw',tend_thradlw)
 call mpas_pool_get_array(tend_physics,'thradsw',tend_thradsw)

!gravity wave drag over orography and convective gravity wave drag:
 call mpas_pool_get_array(tend_physics,'ugwdo',tend_ugwdo)
 call mpas_pool_get_array(tend_physics,'ugwdc',tend_ugwdc)
 call mpas_pool_get_array(tend_physics,'vgwdo',tend_vgwdo)
 call mpas_pool_get_array(tend_physics,'vgwdc',tend_vgwdc)

!total physics tendencies from driver to gbphys:
 call mpas_pool_get_array(tend_physics,'thphysics' ,tend_thphysics )
 call mpas_pool_get_array(tend_physics,'qvphysics' ,tend_qvphysics )
 call mpas_pool_get_array(tend_physics,'qcphysics' ,tend_qcphysics )
 call mpas_pool_get_array(tend_physics,'qo3physics',tend_qo3physics)
 call mpas_pool_get_array(tend_physics,'uphysics'  ,tend_uphysics  )
 call mpas_pool_get_array(tend_physics,'vphysics'  ,tend_vphysics  )

 call mpas_pool_get_array(tend,'u'      ,tend_u    )
 call mpas_pool_get_array(tend,'theta_m',tend_theta)
 call mpas_pool_get_array(tend,'scalars_tend',tend_scalars)
 tend_qv  => tend_scalars(index_qv,:,:)
 tend_qc  => tend_scalars(index_qc,:,:)
 tend_qo3 => tend_scalars(index_qo3,:,:)

!initialize the tendency for the potential temperature:
 allocate(theta(nVertLevels,nCellsSolve)  )
 allocate(tend_th(nVertLevels,nCellsSolve))
 tend_th      = 0.
 tend_scalars = 0.

!add coupled physics tendency to u tendency:
 allocate(tend_uphysics_Edge(nVertLevels,nEdges))
 tend_uphysics_Edge(:,:) = 0._RKIND
 call tend_toEdges(block,mesh,tend_uphysics,tend_vphysics,tend_uphysics_Edge)
 do iEdge = 1, nEdgesSolve
 do k  = 1, nVertLevels
    tend_u(k,iEdge)=tend_u(k,iEdge)+tend_uphysics_Edge(k,iEdge)*mass_edge(k,iEdge)
 enddo
 enddo
 deallocate(tend_uphysics_Edge)

 do iCell = 1, nCellsSolve
 do k = 1, nVertLevels
    tend_th(k,iCell)  = tend_th(k,iCell)  + tend_thphysics(k,iCell)*mass(k,iCell)
    tend_qv(k,iCell)  = tend_qv(k,iCell)  + tend_qvphysics(k,iCell)*mass(k,iCell)
    tend_qc(k,iCell)  = tend_qc(k,iCell)  + tend_qcphysics(k,iCell)*mass(k,iCell)
    tend_qo3(k,iCell) = tend_qo3(k,iCell) + tend_qo3physics(k,iCell)*mass(k,iCell)
 enddo
 enddo

!if non-hydrostatic core, convert the tendency for the potential temperature to a
!tendency for the modified potential temperature:
 do iCell = 1, nCellsSolve
 do k = 1, nVertLevels
    theta(k,iCell) = theta_m(k,iCell) / (1. + R_v/R_d * qv(k,iCell))
    tend_th(k,iCell) = (1. + R_v/R_d * qv(k,iCell)) * tend_th(k,iCell) &
                    + R_v/R_d * theta(k,iCell) * tend_qv(k,iCell)
    tend_theta(k,iCell) = tend_theta(k,iCell) + tend_th(k,iCell)
 enddo
 enddo
 deallocate(theta)
 deallocate(tend_th)

!if(rk_step .eq. 3) then
    write(0,*)
    write(0,*) '--- enter subroutine physics_addtend:'
    write(0,*) 'max tend_thphysics  = ',maxval(tend_thphysics(:,1:nCellsSolve) )
    write(0,*) 'min tend_thphysics  = ',minval(tend_thphysics(:,1:nCellsSolve) )
    write(0,*) 'max tend_qvphysics  = ',maxval(tend_qvphysics(:,1:nCellsSolve) )
    write(0,*) 'min tend_qvphysics  = ',minval(tend_qvphysics(:,1:nCellsSolve) )
    write(0,*) 'max tend_qcphysics  = ',maxval(tend_qcphysics(:,1:nCellsSolve) )
    write(0,*) 'min tend_qcphysics  = ',minval(tend_qcphysics(:,1:nCellsSolve) )
    write(0,*) 'max tend_qo3physics = ',maxval(tend_qo3physics(:,1:nCellsSolve))
    write(0,*) 'min tend_qo3physics = ',minval(tend_qo3physics(:,1:nCellsSolve))
    write(0,*) '--- end subroutine physics_addtend'
    write(0,*)
!   stop
!endif

 end subroutine  physics_addtend

!=================================================================================================================
 subroutine tend_toEdges(block,mesh,Ux_tend,Uy_tend,U_tend)
!=================================================================================================================

!input arguments:
!----------------
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 real(kind=RKIND),intent(in),dimension(:,:):: Ux_tend,Uy_tend 

!output arguments:
!-----------------
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
!-----------------
 type (field2DReal), pointer :: tempField
 type (field2DReal), target :: tempFieldTarget
 integer:: iCell,iEdge,k,j
 integer,pointer:: nCells,nCellsSolve,nVertLevels
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: edgesOnCell

 real(kind=RKIND),dimension(:,:),pointer:: east,north,edge_normal
 real(kind=RKIND),dimension(:,:),allocatable,target:: Ux_tend_halo,Uy_tend_halo
 
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

 call mpas_pool_get_array(mesh, 'east', east)
 call mpas_pool_get_array(mesh, 'north', north)
 call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'edgeNormalVectors', edge_normal)

 allocate(Ux_tend_halo(nVertLevels,nCells+1))
 allocate(Uy_tend_halo(nVertLevels,nCells+1))

 Ux_tend_halo(:,:) = 0.
 Uy_tend_halo(:,:) = 0.
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       Ux_tend_halo(k,iCell) = Ux_tend(k,iCell)
       Uy_tend_halo(k,iCell) = Uy_tend(k,iCell)
    enddo
 enddo

 tempField => tempFieldTarget 
 tempField % block => block
 tempField % dimSizes(1) = nVertLevels
 tempField % dimSizes(2) = nCellsSolve
 tempField % sendList => block % parinfo % cellsToSend
 tempField % recvList => block % parinfo % cellsToRecv
 tempField % copyList => block % parinfo % cellsToCopy
 tempField % prev => null()
 tempField % next => null()

 tempField % array => Ux_tend_halo
 call mpas_dmpar_exch_halo_field(tempField)
 
 tempField % array => Uy_tend_halo
 call mpas_dmpar_exch_halo_field(tempField)

 U_tend(:,:) = 0.0
 do iCell = 1, nCells
 do j = 1, nEdgesOnCell(iCell)
    iEdge = edgesOnCell(j,iCell)
    do k = 1, nVertLevels
       U_tend(k,iEdge) = U_tend(k,iEdge) & 
                + 0.5 * Ux_tend_halo(k,iCell) * (edge_normal(1,iEdge) * east(1,iCell)   &
                                              +  edge_normal(2,iEdge) * east(2,iCell)   &
                                              +  edge_normal(3,iEdge) * east(3,iCell))  &
                + 0.5 * Uy_tend_halo(k,iCell) * (edge_normal(1,iEdge) * north(1,iCell)  &
                                              +  edge_normal(2,iEdge) * north(2,iCell)  &
                                              +  edge_normal(3,iEdge) * north(3,iCell))
    enddo
 enddo
 enddo

 deallocate(Ux_tend_halo)
 deallocate(Uy_tend_halo)
 
 end subroutine tend_toEdges

!=================================================================================================================
 end module mpas_atmphys_todynamics
!=================================================================================================================
