! DH* 20160818
#ifdef MPAS_ALTOMP_TIME_INTEGRATION
! *DH 20160818
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_todynamics
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_dmpar

 use mpas_atmphys_constants, only: R_d,R_v,degrad

 implicit none
 private
 public:: physics_addtend


!Interface between the physics parameterizations and the non-hydrostatic dynamical core.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
!
! subroutines in mpas_atmphys_todynamics:
! ---------------------------------------
! physics_addtend: add and mass-weigh tendencies before being added to dynamics tendencies.
! tend_toEdges   : interpolate wind-tendencies from centers to edges of grid-cells.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * added calculation of the advective tendency of the potential temperature due to horizontal
!   and vertical advection, and horizontal mixing (diffusion).
!   Laura D. Fowler (birch.mmm.ucar.edu) / 2013-11-19.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.
! * renamed "tiedtke" with "cu_tiedtke".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-22.
! * modified the sourcecode to accomodate the packages "cu_kain_fritsch_in" and "cu_tiedtke_in".
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-24.
! * added the calculation of rthdynten which is the tendency of potential temperature due to horizontal and
!   vertical advections needed in the Grell-Freitas scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-30.
! * added the option bl_mynn for the calculation of the tendency for the cloud ice number concentration.
!   Laura D. Fowler (laura@ucar.edu) / 2016-04-11.


 contains

 
!=================================================================================================================
 subroutine physics_addtend(block, mesh, state, diag, tend, tend_physics, configs, rk_step, dynamics_substep, &
                            theta, tend_th, Ux_tend_halo, Uy_tend_halo, rublten_Edge, rucuten_Edge, &
                            cellStart, cellEnd, edgeStart, edgeEnd, &
                            cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd)
!=================================================================================================================

 use mpas_atm_dimensions

!input variables:
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: configs
 integer, intent(in):: rk_step
 integer, intent(in):: dynamics_substep

 real(kind=RKIND),dimension(:,:),intent(inout):: theta, tend_th
 real(kind=RKIND),dimension(:,:),target,intent(inout):: Ux_tend_halo,Uy_tend_halo
 real(kind=RKIND),dimension(:,:),intent(inout):: rublten_Edge, rucuten_Edge
 integer, intent(in):: cellStart, cellEnd, edgeStart, edgeEnd
 integer, intent(in):: cellSolveStart, cellSolveEnd, edgeSolveStart, edgeSolveEnd

!inout variables:
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: tend
 type(mpas_pool_type),intent(inout):: tend_physics

!local variables:
 character(len=StrKIND), pointer :: config_pbl_scheme, config_convection_scheme, &
                                    config_radt_lw_scheme, config_radt_sw_scheme

 integer:: i,iCell,k,n
 integer,pointer:: index_qv, index_qc, index_qr, index_qi, index_qs, index_qg
 integer,pointer:: index_ni
 integer,pointer:: nCells,nCellsSolve,nEdges,nEdgesSolve

 ! required for tend_toEdges
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: edgesOnCell
 real(kind=RKIND),dimension(:,:),pointer:: east,north,edge_normal

 real(kind=RKIND),dimension(:,:),pointer:: mass          ! time level 2 rho_zz
 real(kind=RKIND),dimension(:,:),pointer:: mass_edge     ! diag rho_edge
 real(kind=RKIND),dimension(:,:),pointer:: theta_m,qv    ! time level 1
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rqvblten,rqcblten, &
                                           rqiblten,rublten,rvblten
 real(kind=RKIND),dimension(:,:),pointer:: rniblten
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten, &
                                           rqrcuten,rqicuten,rqscuten, &
                                           rucuten,rvcuten
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw                                    
 real(kind=RKIND),dimension(:,:),pointer:: rthdynten
 
 real(kind=RKIND),dimension(:,:),pointer:: tend_u_phys !nick
 real(kind=RKIND),dimension(:,:),pointer  :: tend_theta,tend_theta_euler,tend_diabatic,tend_u
 real(kind=RKIND),dimension(:,:,:),pointer:: tend_scalars
 real(kind=RKIND),dimension(:,:),pointer  :: tend_qv,tend_qc,tend_qr,tend_qi,tend_qs,tend_qg
 real(kind=RKIND),dimension(:,:),pointer  :: tend_ni

 real(kind=RKIND):: coeff
 real(kind=RKIND), parameter :: R_v_over_R_d = R_v/R_d

!=================================================================================================================
 call mpas_pool_get_dimension(mesh, 'nCells', nCells)
 call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
 call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
 call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
 call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

 call mpas_pool_get_config(configs, 'config_pbl_scheme', config_pbl_scheme)
 call mpas_pool_get_config(configs, 'config_convection_scheme', config_convection_scheme)
 call mpas_pool_get_config(configs, 'config_radt_lw_scheme', config_radt_lw_scheme)
 call mpas_pool_get_config(configs, 'config_radt_sw_scheme', config_radt_sw_scheme)

 call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
 call mpas_pool_get_array(state, 'scalars', scalars, 1)
 call mpas_pool_get_array(state, 'rho_zz', mass, 2)
 call mpas_pool_get_array(diag , 'rho_edge', mass_edge)

 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_dimension(state,'index_qc',index_qc)
 call mpas_pool_get_dimension(state,'index_qr',index_qr)
 call mpas_pool_get_dimension(state,'index_qi',index_qi)
 call mpas_pool_get_dimension(state,'index_qs',index_qs)
 call mpas_pool_get_dimension(state,'index_qg',index_qg)

 call mpas_pool_get_array(tend,'rt_diabatic_tend',tend_diabatic)

 call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
 call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)

 call mpas_pool_get_array(tend,'u'           , tend_u         )
 call mpas_pool_get_array(tend,'theta_m'     , tend_theta     )
 call mpas_pool_get_array(tend,'theta_euler' ,tend_theta_euler)
 call mpas_pool_get_array(tend,'scalars_tend',tend_scalars    )

 
 call mpas_pool_get_array(diag , 'tend_u_phys', tend_u_phys) !nick
 call mpas_pool_get_array(mesh, 'east', east)
 call mpas_pool_get_array(mesh, 'north', north)
 call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
 call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
 call mpas_pool_get_array(mesh, 'edgeNormalVectors', edge_normal)

 qv => scalars(index_qv,:,:)
 tend_qv => tend_scalars(index_qv,:,:)
 tend_qc => tend_scalars(index_qc,:,:)
 tend_qr => tend_scalars(index_qr,:,:)
 tend_qi => tend_scalars(index_qi,:,:)
 tend_qs => tend_scalars(index_qs,:,:)
 tend_qg => tend_scalars(index_qg,:,:)

!initialize the tendency for the potential temperature and all scalars
!due to PBL, convection, and longwave and shortwave radiation:
 tend_qv(:,cellStart:cellEnd) = 0._RKIND
 tend_qc(:,cellStart:cellEnd) = 0._RKIND
 tend_qr(:,cellStart:cellEnd) = 0._RKIND
 tend_qi(:,cellStart:cellEnd) = 0._RKIND
 tend_qs(:,cellStart:cellEnd) = 0._RKIND
 tend_qg(:,cellStart:cellEnd) = 0._RKIND
 tend_th(:,cellSolveStart:cellSolveEnd) = 0._RKIND

!add coupled tendencies due to PBL processes:
 if(config_pbl_scheme .ne. 'off') then
    call mpas_pool_get_array(tend_physics,'rublten' ,rublten )
    call mpas_pool_get_array(tend_physics,'rvblten' ,rvblten )
    call mpas_pool_get_array(tend_physics,'rthblten',rthblten)
    call mpas_pool_get_array(tend_physics,'rqvblten',rqvblten)
    call mpas_pool_get_array(tend_physics,'rqcblten',rqcblten)
    call mpas_pool_get_array(tend_physics,'rqiblten',rqiblten)

    ! Barriers implemented in tend_toEdges, in particular on exit
    call tend_toEdges(block,mesh,rublten,rvblten,rublten_Edge, &
                      nCells, nCellsSolve, nVertLevels, east, north, edgesOnCell, &
                      nEdgesOnCell, edge_normal, Ux_tend_halo, Uy_tend_halo, &
                      cellStart, cellEnd, edgeStart, edgeEnd, cellSolveStart, cellSolveEnd)

    tend_u_phys(1:nVertLevels,edgeStart:edgeEnd) = rublten_Edge(1:nVertLevels,edgeStart:edgeEnd)
    
    do i = edgeSolveStart, edgeSolveEnd
!DIR$ IVDEP
       do k  = 1, nVertLevels
          tend_u(k,i) = tend_u(k,i) + rublten_Edge(k,i)*mass_edge(k,i)
       enddo
    enddo

    do i = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthblten(k,i)*mass(k,i)
          tend_qv(k,i) = tend_qv(k,i) + rqvblten(k,i)*mass(k,i)
          tend_qc(k,i) = tend_qc(k,i) + rqcblten(k,i)*mass(k,i)
          tend_qi(k,i) = tend_qi(k,i) + rqiblten(k,i)*mass(k,i)
       enddo
    enddo

    pbl_select: select case (config_pbl_scheme)
    
       case("bl_mynn")
          call mpas_pool_get_dimension(state,'index_ni',index_ni)
          tend_ni => tend_scalars(index_ni,:,:)
          tend_ni(:,cellStart:cellEnd) = 0._RKIND

          call mpas_pool_get_array(tend_physics,'rniblten',rniblten)

!$OMP barrier

          do i = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
             do k = 1, nVertLevels
                tend_ni(k,i) = tend_ni(k,i) + rniblten(k,i)*mass(k,i)
             enddo
          enddo

       case default

     end select pbl_select

 endif

!add coupled tendencies due to convection:
 if(config_convection_scheme .ne. 'off') then
    call mpas_pool_get_array(tend_physics,'rthcuten',rthcuten)
    call mpas_pool_get_array(tend_physics,'rqvcuten',rqvcuten)
    call mpas_pool_get_array(tend_physics,'rqccuten',rqccuten)
    call mpas_pool_get_array(tend_physics,'rqicuten',rqicuten)

    if(config_convection_scheme .eq. 'cu_tiedtke') then
       call mpas_pool_get_array(tend_physics,'rucuten' ,rucuten )
       call mpas_pool_get_array(tend_physics,'rvcuten' ,rvcuten )
       
       ! Barriers implemented in tend_toEdges, in particular on exit
       call tend_toEdges(block,mesh,rucuten,rvcuten,rucuten_Edge, &
                         nCells, nCellsSolve, nVertLevels, east, north, edgesOnCell, &
                         nEdgesOnCell, edge_normal, Ux_tend_halo, Uy_tend_halo, &
                         cellStart, cellEnd, edgeStart, edgeEnd, cellSolveStart, cellSolveEnd)

       tend_u_phys(1:nVertLevels,edgeStart:edgeEnd) = tend_u_phys(1:nVertLevels,edgeStart:edgeEnd) + &
                                                    & rucuten_Edge(1:nVertLevels,edgeStart:edgeEnd)

       do i = edgeSolveStart, edgeSolveEnd
!DIR$ IVDEP
          do k  = 1, nVertLevels
             tend_u(k,i)=tend_u(k,i)+rucuten_Edge(k,i)*mass_edge(k,i)
          enddo
       enddo

    endif

    do i = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthcuten(k,i)*mass(k,i)
          tend_qv(k,i) = tend_qv(k,i) + rqvcuten(k,i)*mass(k,i)
          tend_qc(k,i) = tend_qc(k,i) + rqccuten(k,i)*mass(k,i)
          tend_qi(k,i) = tend_qi(k,i) + rqicuten(k,i)*mass(k,i)
       enddo
    enddo

    convection_select: select case(config_convection_scheme)
    
       case('cu_grell_freitas')
          call mpas_pool_get_array(tend_physics, 'rthdynten' , rthdynten)

          do i = cellSolveStart, cellSolveEnd
             do k = 1, nVertLevels
                coeff = 1._RKIND/(1._RKIND + R_v_over_R_d*qv(k,i))
                rthdynten(k,i) = coeff*(tend_theta(k,i)-tend_theta_euler(k,i) &
                               - mass(k,i)*tend_diabatic(k,i))/mass(k,i)
             enddo
          enddo

       case('cu_kain_fritsch')
          call mpas_pool_get_array(tend_physics,'rqrcuten',rqrcuten)
          call mpas_pool_get_array(tend_physics,'rqscuten',rqscuten)
          do i = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
             do k = 1, nVertLevels
                tend_qr(k,i) = tend_qr(k,i) + rqrcuten(k,i)*mass(k,i)
                tend_qs(k,i) = tend_qs(k,i) + rqscuten(k,i)*mass(k,i)
             enddo
          enddo
    
       case('cu_tiedtke')
          call mpas_pool_get_array(tend_physics,'rucuten' ,rucuten )
          call mpas_pool_get_array(tend_physics,'rvcuten' ,rvcuten )

          ! Barriers implemented in tend_toEdges, in particular on exit
          call tend_toEdges(block,mesh,rucuten,rvcuten,rucuten_Edge, &
                            nCells, nCellsSolve, nVertLevels, east, north, edgesOnCell, &
                            nEdgesOnCell, edge_normal, Ux_tend_halo, Uy_tend_halo, &
                            cellStart, cellEnd, edgeStart, edgeEnd, cellSolveStart, cellSolveEnd)

          
          do i = edgeSolveStart, edgeSolveEnd
!DIR$ IVDEP
             do k  = 1, nVertLevels
                tend_u(k,i) = tend_u(k,i) + rucuten_Edge(k,i)*mass_edge(k,i)
             enddo
          enddo

        case default

    end select convection_select

 endif

!add coupled tendencies due to longwave radiation:
 if(config_radt_lw_scheme .ne. 'off') then
    do i = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthratenlw(k,i)*mass(k,i)
       enddo
    enddo
 endif

!add coupled tendencies due to shortwave radiation:
 if(config_radt_sw_scheme .ne. 'off') then
    do i = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
       do k = 1, nVertLevels
          tend_th(k,i) = tend_th(k,i) + rthratensw(k,i)*mass(k,i)
       enddo
    enddo
 endif

!if non-hydrostatic core, convert the tendency for the potential temperature to a
!tendency for the modified potential temperature:
 do i = cellSolveStart, cellSolveEnd
!DIR$ IVDEP
    do k = 1, nVertLevels
       theta(k,i) = theta_m(k,i) / (1. + R_v_over_R_d * qv(k,i))
       tend_th(k,i) = (1. + R_v_over_R_d * qv(k,i)) * tend_th(k,i) &
                       + R_v_over_R_d * theta(k,i) * tend_qv(k,i)
       tend_theta(k,i) = tend_theta(k,i) + tend_th(k,i)
    enddo
 enddo

 if(rk_step .eq. 3) then
 !$OMP barrier
    if(mpas_threading_get_thread_num() .eq. 0) then
       write(0,*)
       write(0,*) '--- enter subroutine physics_addtend:'
       write(0,*) 'max rthblten   = ',maxval(rthblten(:,:))
       write(0,*) 'min rthblten   = ',minval(rthblten(:,:))
       write(0,*) 'max rthcuten   = ',maxval(rthcuten(:,:))
       write(0,*) 'min rthcuten   = ',minval(rthcuten(:,:))
       write(0,*) 'max rthratenlw = ',maxval(rthratenlw(:,:))
       write(0,*) 'min rthratenlw = ',minval(rthratenlw(:,:))
       write(0,*) 'max rthratensw = ',maxval(rthratensw(:,:))
       write(0,*) 'min rthratensw = ',minval(rthratensw(:,:))
       
       write(0,*) 'mean rthblten   = ',sum(rthblten(:,:))/(max(1,size(rthblten(:,:))))
       write(0,*) 'mean rthcuten   = ',sum(rthcuten(:,:))/(max(1,size(rthcuten(:,:))))
       write(0,*) 'mean rthratenlw = ',sum(rthratenlw(:,:))/(max(1,size(rthratenlw(:,:))))
       write(0,*) 'mean rthratensw = ',sum(rthratensw(:,:))/(max(1,size(rthratensw(:,:))))
       
       write(0,*) 'max tend_qv = ',   maxval(tend_qv(:,:))
       write(0,*) 'min tend_qv = ',   minval(tend_qv(:,:))
       write(0,*) 'max tend_qc = ',   maxval(tend_qc(:,:))
       write(0,*) 'min tend_qc = ',   minval(tend_qc(:,:))
       write(0,*) 'max tend_qr = ',   maxval(tend_qr(:,:))
       write(0,*) 'min tend_qr = ',   minval(tend_qr(:,:))
       write(0,*) 'max tend_qi = ',   maxval(tend_qi(:,:))
       write(0,*) 'min tend_qi = ',   minval(tend_qi(:,:))
       write(0,*) 'max tend_qs = ',   maxval(tend_qs(:,:))
       write(0,*) 'min tend_qs = ',   minval(tend_qs(:,:))
       write(0,*) 'max tend_qg = ',   maxval(tend_qg(:,:))
       write(0,*) 'min tend_qg = ',   minval(tend_qg(:,:))
       write(0,*) 'max tend_th = ',   maxval(tend_th(:,:))
       write(0,*) 'min tend_th = ',   minval(tend_th(:,:))
       write(0,*) 'max tend_u     = ',maxval(tend_u(:,:))
       write(0,*) 'min tend_u     = ',minval(tend_u(:,:))
       write(0,*) 'max tend_th    = ',maxval(tend_th(:,:))
       write(0,*) 'min tend_th    = ',minval(tend_th(:,:))
       write(0,*) 'max tend_theta = ',maxval(tend_theta(:,:))
       write(0,*) 'min tend_theta = ',minval(tend_theta(:,:))
       write(0,*) 'max theta      = ',maxval(theta(:,:))
       write(0,*) 'min theta      = ',minval(theta(:,:))
       
       write(0,*) 'mean tend_qv = ',   sum(tend_qv(:,:))/(max(1,size(tend_qv(:,:))))
       write(0,*) 'mean tend_qc = ',   sum(tend_qc(:,:))/(max(1,size(tend_qc(:,:))))
       write(0,*) 'mean tend_qr = ',   sum(tend_qr(:,:))/(max(1,size(tend_qr(:,:))))
       write(0,*) 'mean tend_qi = ',   sum(tend_qi(:,:))/(max(1,size(tend_qi(:,:))))
       write(0,*) 'mean tend_qs = ',   sum(tend_qs(:,:))/(max(1,size(tend_qs(:,:))))
       write(0,*) 'mean tend_qg = ',   sum(tend_qg(:,:))/(max(1,size(tend_qg(:,:))))
       write(0,*) 'mean tend_th = ',   sum(tend_th(:,:))/(max(1,size(tend_th(:,:))))
       write(0,*) 'mean tend_u     = ',sum(tend_u(:,:))/(max(1,size(tend_u(:,:))))
       write(0,*) 'mean tend_th    = ',sum(tend_th(:,:))/(max(1,size(tend_th(:,:))))
       write(0,*) 'mean tend_theta = ',sum(tend_theta(:,:))/(max(1,size(tend_theta(:,:))))
       write(0,*) 'mean theta      = ',sum(theta(:,:))/(max(1,size(theta(:,:))))
       
       write(0,*) '--- end subroutine physics_addtend'
       write(0,*)
    endif
 endif

!formats:
 201 format(2i6,10(1x,e15.8))
 202 format(3i6,10(1x,e15.8))

 end subroutine physics_addtend

!=================================================================================================================
subroutine tend_toEdges(block,mesh,Ux_tend,Uy_tend,U_tend, &
                        nCells, nCellsSolve, nVertLevels, east, north, edgesOnCell, &
                        nEdgesOnCell, edge_normal, Ux_tend_halo, Uy_tend_halo, &
                        cellStart, cellEnd, edgeStart, edgeEnd, cellSolveStart, cellSolveEnd)
!=================================================================================================================

!input arguments:
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: mesh
 real(kind=RKIND),intent(in),dimension(:,:):: Ux_tend,Uy_tend 
 integer,pointer,intent(in):: nCells,nCellsSolve,nVertLevels
 integer,dimension(:),pointer,intent(in)  :: nEdgesOnCell
 integer,dimension(:,:),pointer,intent(in):: edgesOnCell
 real(kind=RKIND),dimension(:,:),pointer,intent(in):: east,north,edge_normal
 integer, intent(in) :: cellStart, cellEnd, edgeStart, edgeEnd
 integer, intent(in) :: cellSolveStart, cellSolveEnd

 !output arguments:
 real(kind=RKIND),dimension(:,:),target,intent(out):: Ux_tend_halo,Uy_tend_halo
 real(kind=RKIND),intent(out),dimension(:,:):: U_tend

!local variables:
 type (field2DReal), pointer :: tempField
 type (field2DReal), target :: tempFieldTarget
 integer:: iCell,iEdge,k,j
 
!-----------------------------------------------------------------------------------------------------------------

 U_tend(:,edgeStart:edgeEnd)       = 0.
 
 do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
    do k = 1, nVertLevels
       Ux_tend_halo(k,iCell) = Ux_tend(k,iCell)
       Uy_tend_halo(k,iCell) = Uy_tend(k,iCell)
    enddo
 enddo
 
 ! Communication block - master only, needs barriers
 ! before and after for data consistency
!$OMP barrier
!$OMP master
 tempField => tempFieldTarget 
 tempField % block => block
 tempField % dimSizes(1) = nVertLevels
 tempField % dimSizes(2) = nCellsSolve
 tempField % sendList => block % parinfo % cellsToSend
 tempField % recvList => block % parinfo % cellsToRecv
 tempField % copyList => block % parinfo % cellsToCopy
 tempField % prev => null()
 tempField % next => null()
 tempField % isActive = .true.

 tempField % array => Ux_tend_halo
 call mpas_dmpar_exch_halo_field(tempField)

 tempField % array => Uy_tend_halo
 call mpas_dmpar_exch_halo_field(tempField)
! DH* 20160818 reactivate those and remove the ones below the do iCell
! loop to parallelise the U_tend calculation for all edges; this results
! in different results for different numbers of MPI tasks due to tiny
! rounding differences when adding up the contribution in different order,
! but gives a significant speedup for the calls to physics_addtend!
! -$- OMP end master
! -$- OMP barrier
! *DH 20160818
 
 ! DH* 20160818
 !do iCell = cellStart,cellEnd
 do iCell = 1, nCells
 ! *DH 20160618
    do j = 1, nEdgesOnCell(iCell)
       iEdge = edgesOnCell(j,iCell)
!DIR$ IVDEP
       do k = 1, nVertLevels
          U_tend(k,iEdge) = U_tend(k,iEdge) & 
                   + 0.5 * Ux_tend_halo(k,iCell) * (edge_normal(1,iEdge) * east(1,iCell)   &
                                                 +  edge_normal(2,iEdge) * east(2,iCell)   &
                                                 +  edge_normal(3,iEdge) * east(3,iCell))  &
                   + 0.5 * Uy_tend_halo(k,iCell) * (edge_normal(1,iEdge) * north(1,iCell)  &
                                                 +  edge_normal(2,iEdge) * north(2,iCell)  &
                                                 +  edge_normal(3,iEdge) * north(3,iCell))
       end do
    end do
 end do
! DH* 20160818 
!$OMP end master
! *DH 20160818
!$OMP barrier
 end subroutine tend_toEdges

!=================================================================================================================
 end module mpas_atmphys_todynamics
!=================================================================================================================
! DH* 20160818
#endif
! *DH 20160818
