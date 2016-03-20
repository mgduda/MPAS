! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_todynamics
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_dmpar

 use mpas_atmphys_nggps_constants, only: R_d,R_v,degrad

 implicit none
 private
 public:: physics_addtend

!modifications to sourcecode:
!----------------------------
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.

 contains
 

!==================================================================================================
 subroutine physics_addtend(block,mesh,state,diag,tend,nggps_physics_tendencies,rk_step,do_full_tend)
!==================================================================================================

 use mpas_atm_dimensions

!input variables:
 type(block_type),intent(in)    :: block
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: nggps_physics_tendencies
 logical, intent(inout) :: do_full_tend

 integer, intent(in):: rk_step

!inout variables:
 type(mpas_pool_type),intent(inout):: tend

!local pointers:
 integer,pointer:: nCells,nCellsSolve,nEdges,nEdgesSolve
 integer,pointer:: index_qv,index_qc,index_qo3

 real(kind=RKIND),dimension(:,:),pointer  :: mass
 real(kind=RKIND),dimension(:,:),pointer  :: mass_edge
 real(kind=PHYSKIND),dimension(:,:),pointer  :: qv
 real(kind=RKIND),dimension(:,:),pointer  :: theta_m
 real(kind=PHYSKIND),dimension(:,:,:),pointer:: scalars

 real(kind=PHYSKIND),dimension(:,:),pointer  :: tend_thdiffusion,tend_qvdiffusion,  &
                                             tend_qcdiffusion,tend_qo3diffusion, &
                                             tend_udiffusion ,tend_vdiffusion 
 real(kind=PHYSKIND),dimension(:,:),pointer  :: tend_thcudeep,tend_qvcudeep,        &
                                             tend_qccudeep,tend_ucudeep,         &
                                             tend_vcudeep
 real(kind=PHYSKIND),dimension(:,:),pointer  :: tend_thcushallow,tend_qvcushallow,  &
                                             tend_qccushallow
 real(kind=PHYSKIND),dimension(:,:),pointer  :: tend_thradlw,tend_thradsw
 real(kind=PHYSKIND),dimension(:,:),pointer  :: tend_ugwdo,tend_vgwdo
 real(kind=PHYSKIND),dimension(:,:),pointer  :: tend_ugwdc,tend_vgwdc

 real(kind=PHYSKIND),dimension(:,:),pointer  :: tend_thphysics,tend_qvphysics,  &
                                             tend_qcphysics,tend_qo3physics, &
                                             tend_uphysics ,tend_vphysics

 real(kind=RKIND),dimension(:,:),pointer  :: tend_theta,tend_u
 real(kind=PHYSKIND),dimension(:,:,:),pointer:: tend_scalars

 real(kind=RKIND),dimension(:,:),pointer:: tend_uphysics_Edge
 real(kind=RKIND),dimension(:,:),allocatable:: theta
 real(kind=RKIND),dimension(:,:),pointer:: tend_th

!local variables:
 integer:: iCell,iEdge,k
 real(kind=RKIND) :: rvord

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells)
 call mpas_pool_get_dimension(mesh,'nEdges'     ,nEdges)
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nEdgesSolve',nEdgesSolve)

 call mpas_pool_get_dimension(state,'index_qv' ,index_qv )
 call mpas_pool_get_dimension(state,'index_qc' ,index_qc )
 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)

 call mpas_pool_get_array(diag ,'rho_edge',mass_edge )
 call mpas_pool_get_array(state,'rho_zz'  ,mass, 2   )
 call mpas_pool_get_array(state,'theta_m' ,theta_m,1 )
 call mpas_pool_get_array(state,'scalars' ,scalars,1 )
 qv => scalars(index_qv,:,:)

!PBL and diffusion:
 call mpas_pool_get_array(nggps_physics_tendencies,'udiffusion'  ,tend_udiffusion  )
 call mpas_pool_get_array(nggps_physics_tendencies,'vdiffusion'  ,tend_vdiffusion  )
 call mpas_pool_get_array(nggps_physics_tendencies,'thdiffusion' ,tend_thdiffusion )
 call mpas_pool_get_array(nggps_physics_tendencies,'qvdiffusion' ,tend_qvdiffusion )
 call mpas_pool_get_array(nggps_physics_tendencies,'qcdiffusion' ,tend_qcdiffusion )
 call mpas_pool_get_array(nggps_physics_tendencies,'qo3diffusion',tend_qo3diffusion)

!deep convection and moist convective adjustment:
 call mpas_pool_get_array(nggps_physics_tendencies,'ucudeep' ,tend_ucudeep )
 call mpas_pool_get_array(nggps_physics_tendencies,'vcudeep' ,tend_vcudeep )
 call mpas_pool_get_array(nggps_physics_tendencies,'thcudeep',tend_thcudeep)
 call mpas_pool_get_array(nggps_physics_tendencies,'qvcudeep',tend_qvcudeep)
 call mpas_pool_get_array(nggps_physics_tendencies,'qccudeep',tend_qccudeep)

!shallow convection:
 call mpas_pool_get_array(nggps_physics_tendencies,'thcushallow',tend_thcushallow)
 call mpas_pool_get_array(nggps_physics_tendencies,'qvcushallow',tend_qvcushallow)
 call mpas_pool_get_array(nggps_physics_tendencies,'qccushallow',tend_qccushallow)

!long-wave and short-wave convection:
 call mpas_pool_get_array(nggps_physics_tendencies,'thradlw',tend_thradlw)
 call mpas_pool_get_array(nggps_physics_tendencies,'thradsw',tend_thradsw)

!gravity wave drag over orography and convective gravity wave drag:
 call mpas_pool_get_array(nggps_physics_tendencies,'ugwdo',tend_ugwdo)
 call mpas_pool_get_array(nggps_physics_tendencies,'ugwdc',tend_ugwdc)
 call mpas_pool_get_array(nggps_physics_tendencies,'vgwdo',tend_vgwdo)
 call mpas_pool_get_array(nggps_physics_tendencies,'vgwdc',tend_vgwdc)

!total physics tendencies from driver to gbphys:
 call mpas_pool_get_array(nggps_physics_tendencies,'tend_th'   ,tend_th )
 call mpas_pool_get_array(nggps_physics_tendencies,'tend_uphysics_Edge' ,tend_uphysics_Edge )
 call mpas_pool_get_array(nggps_physics_tendencies,'thphysics' ,tend_thphysics )
 call mpas_pool_get_array(nggps_physics_tendencies,'qvphysics' ,tend_qvphysics )
 call mpas_pool_get_array(nggps_physics_tendencies,'qcphysics' ,tend_qcphysics )
 call mpas_pool_get_array(nggps_physics_tendencies,'qo3physics',tend_qo3physics)
 call mpas_pool_get_array(nggps_physics_tendencies,'uphysics'  ,tend_uphysics  )
 call mpas_pool_get_array(nggps_physics_tendencies,'vphysics'  ,tend_vphysics  )

!dynamic tendencies:
 call mpas_pool_get_array(tend,'u',tend_u)
 call mpas_pool_get_array(tend,'theta_m',tend_theta)
 call mpas_pool_get_array(tend,'scalars_tend',tend_scalars)

!initialize the tendency for the potential temperature and all scalars due to PBL, deep and
!shallow convection, and longwave and shortwave radiation:

!add coupled physics tendency to u tendency:
 if (do_full_tend) then
    call tend_toEdges(block,mesh,tend_uphysics,tend_vphysics,tend_uphysics_Edge)
 end if

 rvord = R_v / R_d

 allocate(theta(nVertLevels,nCellsSolve)  )
 call physics_addtend_work(nCells, nEdges, nCellsSolve, nEdgesSolve, tend_u, tend_uphysics_Edge, mass_edge, &
                           tend_th, tend_scalars, scalars, index_qv, index_qc, index_qo3, rvord, &
                           tend_thphysics, tend_qvphysics, tend_qcphysics, tend_qo3physics, mass, &
                           theta, theta_m, tend_theta &
                           )
 deallocate(theta)

 if (do_full_tend) then
    do_full_tend = .false.
 end if

 end subroutine  physics_addtend


 subroutine physics_addtend_work(nCells, nEdges, nCellsSolve, nEdgesSolve, tend_u, tend_uphysics_Edge, mass_edge, &
                                 tend_th, tend_scalars, scalars, index_qv, index_qc, index_qo3, rvord, &
                                 tend_thphysics, tend_qvphysics, tend_qcphysics, tend_qo3physics, mass, &
                                 theta, theta_m, tend_theta &
                                )

    use mpas_atm_dimensions

    implicit none

    integer, intent(in) :: nCells, nEdges, nCellsSolve, nEdgesSolve
    integer, intent(in) :: index_qv, index_qc, index_qo3
    real(kind=RKIND), intent(in) :: rvord
    real(kind=RKIND), dimension(nVertLevels,nEdges), intent(inout) :: tend_u
    real(kind=RKIND), dimension(nVertLevels,nEdges), intent(in) :: tend_uphysics_Edge
    real(kind=RKIND), dimension(nVertLevels,nEdges), intent(in) :: mass_edge
    real(kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: tend_th
    real(kind=PHYSKIND), dimension(num_scalars,nVertLevels,nCells), intent(inout) :: tend_scalars
    real(kind=PHYSKIND), dimension(num_scalars,nVertLevels,nCells), intent(in) :: scalars
    real(kind=PHYSKIND), dimension(nVertLevels,nCells), intent(in) :: tend_thphysics
    real(kind=PHYSKIND), dimension(nVertLevels,nCells), intent(in) :: tend_qvphysics
    real(kind=PHYSKIND), dimension(nVertLevels,nCells), intent(in) :: tend_qcphysics
    real(kind=PHYSKIND), dimension(nVertLevels,nCells), intent(in) :: tend_qo3physics
    real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: mass
    real(kind=RKIND), dimension(nVertLevels,nCellsSolve), intent(inout) :: theta
    real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: theta_m
    real(kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: tend_theta

    integer :: iCell, iEdge, k

    do iEdge = 1, nEdgesSolve
       do k  = 1, nVertLevels
          tend_u(k,iEdge)=tend_u(k,iEdge)+tend_uphysics_Edge(k,iEdge)*mass_edge(k,iEdge)
       end do
    end do

    do iCell = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_scalars(index_qv,k,iCell)  = tend_qvphysics(k,iCell)*mass(k,iCell)
          tend_scalars(index_qc,k,iCell)  = tend_qcphysics(k,iCell)*mass(k,iCell)
          tend_scalars(index_qo3,k,iCell) = tend_qo3physics(k,iCell)*mass(k,iCell)
       end do
    end do

    do iCell = 1, nCellsSolve
       do k = 1, nVertLevels
          tend_th(k,iCell)  = tend_thphysics(k,iCell)*mass(k,iCell)
          theta(k,iCell) = theta_m(k,iCell) / (1.0 + rvord * scalars(index_qv,k,iCell))
          tend_th(k,iCell) = (1.0 + rvord * scalars(index_qv,k,iCell)) * tend_th(k,iCell) &
                          + rvord * theta(k,iCell) * tend_scalars(index_qv,k,iCell)
          tend_theta(k,iCell) = tend_theta(k,iCell) + tend_th(k,iCell)
       end do
    end do

 end subroutine physics_addtend_work


!==================================================================================================
 subroutine tend_toEdges(block,mesh,Ux_tend,Uy_tend,U_tend)
!==================================================================================================

 use mpas_atm_dimensions

!input arguments:
!----------------
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 real(kind=PHYSKIND),intent(in),dimension(:,:):: Ux_tend,Uy_tend 

!output arguments:
!-----------------
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
!-----------------
 type (field2DReal), pointer :: tempField
 type (field2DReal), target :: tempFieldTarget
 integer:: iCell,iEdge,k,j
 integer,pointer:: nCells,nCellsSolve
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: edgesOnCell

 real(kind=RKIND),dimension(:,:),pointer:: east,north,edge_normal
 real(kind=RKIND),dimension(:,:),allocatable,target:: Ux_tend_halo,Uy_tend_halo
 
!--------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)

 call mpas_pool_get_array(mesh, 'east', east)
 call mpas_pool_get_array(mesh, 'north', north)
 call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'edgeNormalVectors', edge_normal)

 allocate(Ux_tend_halo(nVertLevels,nCells+1))
 allocate(Uy_tend_halo(nVertLevels,nCells+1))

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
