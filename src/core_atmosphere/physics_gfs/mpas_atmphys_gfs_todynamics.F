! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_todynamics
 use mpas_configure
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


 contains
 
!==================================================================================================
 subroutine physics_addtend(mesh, state, diag, tend, tend_physics, mass, mass_edge, rk_step)
!==================================================================================================

!input variables:
!----------------
 type(mesh_type),intent(in):: mesh
 type(state_type),intent(in):: state
 type(diag_type),intent(in):: diag
 type(tend_physics_type),intent(inout):: tend_physics
 integer, intent(in):: rk_step
 real(kind=RKIND),dimension(:,:),intent(in):: mass
 real(kind=RKIND),dimension(:,:),intent(in):: mass_edge

!inout variables:
!----------------
 type(tend_type),intent(inout):: tend

!local variables:
!----------------
 type(block_type),pointer :: block

 integer:: iCell,iEdge,k,n,nCellsSolve,nEdges,nEdgesSolve,nVertLevels

 real(kind=RKIND),dimension(:,:),pointer:: theta_m,qv

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

 real(kind=RKIND),dimension(:,:),allocatable:: tend_uphysics_Edge
 real(kind=RKIND),dimension(:,:),allocatable:: theta,tend_th

!==================================================================================================
 block => mesh % block
 
 nEdges      = mesh % nEdges
 nCellsSolve = mesh % nCellsSolve
 nEdgesSolve = mesh % nEdgesSolve
 nVertLevels = mesh % nVertLevels

!theta   => diag % theta % array
 theta_m => state % theta_m % array
 qv      => state % scalars % array(state%index_qv,:,:)

!PBL and diffusion:
 tend_udiffusion   => tend_physics % udiffusion   % array
 tend_vdiffusion   => tend_physics % vdiffusion   % array
 tend_thdiffusion  => tend_physics % thdiffusion  % array
 tend_qvdiffusion  => tend_physics % qvdiffusion  % array
 tend_qcdiffusion  => tend_physics % qcdiffusion  % array
 tend_qo3diffusion => tend_physics % qo3diffusion % array

!deep convection and moist convective adjustment:
 tend_ucudeep  => tend_physics % ucudeep  % array
 tend_vcudeep  => tend_physics % vcudeep  % array
 tend_thcudeep => tend_physics % thcudeep % array
 tend_qvcudeep => tend_physics % qvcudeep % array
 tend_qccudeep => tend_physics % qccudeep % array

!shallow convection:
 tend_thcushallow => tend_physics % thcushallow % array
 tend_qvcushallow => tend_physics % qvcushallow % array
 tend_qccushallow => tend_physics % qccushallow % array

!long-wave and short-wave convection:
 tend_thradlw => tend_physics % thradlw % array
 tend_thradsw => tend_physics % thradsw % array

!gravity wave drag over orography and convective gravity wave drag:
 tend_ugwdo   => tend_physics % ugwdo % array
 tend_ugwdc   => tend_physics % ugwdc % array
 tend_vgwdo   => tend_physics % vgwdo % array
 tend_vgwdc   => tend_physics % vgwdc % array 

!total physics tendencies from driver to gbphys:
 tend_thphysics  => tend_physics % thphysics  % array
 tend_qvphysics  => tend_physics % qvphysics  % array
 tend_qcphysics  => tend_physics % qcphysics  % array
 tend_qo3physics => tend_physics % qo3physics % array
 tend_uphysics   => tend_physics % uphysics   % array
 tend_vphysics   => tend_physics % vphysics   % array

 tend_u       => tend % u % array
 tend_theta   => tend % theta_m % array
 tend_scalars => tend % scalars % array

!initialize the tendency for the potential temperature:
 allocate(theta(nVertLevels,nCellsSolve)  )
 allocate(tend_th(nVertLevels,nCellsSolve))
 tend_th      = 0.
 tend_scalars = 0.

!add coupled physics tendency to u tendency:
 allocate(tend_uphysics_Edge(nVertLevels,nEdges))
 tend_uphysics_Edge(:,:) = 0._RKIND
 call tend_toEdges(mesh,tend_uphysics,tend_vphysics,tend_uphysics_Edge)
 do iEdge = 1, nEdgesSolve
 do k  = 1, nVertLevels
    tend_u(k,iEdge)=tend_u(k,iEdge)+tend_uphysics_Edge(k,iEdge)*mass_edge(k,iEdge)
 enddo
 enddo
 deallocate(tend_uphysics_Edge)

 do iCell = 1, nCellsSolve
 do k = 1, nVertLevels
    tend_th(k,iCell) = tend_th(k,iCell) + tend_thphysics(k,iCell)*mass(k,iCell)
    tend_scalars(tend%index_qv,k,iCell)  = tend_scalars(tend%index_qv,k,iCell)  &
                                       + tend_qvphysics(k,iCell)*mass(k,iCell)
    tend_scalars(tend%index_qc,k,iCell)  = tend_scalars(tend%index_qc,k,iCell)  &
                                       + tend_qcphysics(k,iCell)*mass(k,iCell)
    tend_scalars(tend%index_qo3,k,iCell) = tend_scalars(tend%index_qo3,k,iCell) &
                                       + tend_qo3physics(k,iCell)*mass(k,iCell)
 enddo
 enddo

!if non-hydrostatic core, convert the tendency for the potential temperature to a
!tendency for the modified potential temperature:
 do iCell = 1, nCellsSolve
 do k = 1, nVertLevels
    theta(k,iCell) = theta_m(k,iCell) / (1. + R_v/R_d * qv(k,iCell))
    tend_th(k,iCell) = (1. + R_v/R_d * qv(k,iCell)) * tend_th(k,iCell) &
                    + R_v/R_d * theta(k,iCell) * tend_scalars(tend%index_qv,k,iCell)
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

!==================================================================================================
 subroutine tend_toEdges(mesh,Ux_tend,Uy_tend,U_tend)
!==================================================================================================

!input arguments:
!----------------
 type(mesh_type),intent(in):: mesh
 real(kind=RKIND),intent(in),dimension(:,:):: Ux_tend,Uy_tend 

!output arguments:
!-----------------
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
!-----------------
 type(block_type),pointer :: block
 type (field2DReal), pointer :: tempField
 type (field2DReal), target :: tempFieldTarget
 integer:: iCell,iEdge,k,j,nCells,nCellsSolve,nVertLevels
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: edgesOnCell

 real(kind=RKIND),dimension(:,:),pointer:: east,north,edge_normal
 real(kind=RKIND),dimension(:,:),allocatable,target:: Ux_tend_halo,Uy_tend_halo
 
!--------------------------------------------------------------------------------------------------

 block => mesh % block
 
 nCells       = mesh % nCells
 nCellsSolve  = mesh % nCellsSolve
 nVertLevels  = mesh % nVertLevels

 east  => mesh % east  % array
 north => mesh % north % array
 edgesOnCell  => mesh % edgesOnCell  % array
 nEdgesOnCell => mesh % nEdgesOnCell % array
 edge_normal  => mesh % edgeNormalVectors % array

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

!==================================================================================================
 end module mpas_atmphys_todynamics
!==================================================================================================
