   subroutine atm_compute_dyn_tend(tend, state, diag, mesh, configs, nVertLevels, rk_step, dt, &
                                   cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                   cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: state - current model state
   !        mesh - grid metadata
   !        diag - some grid diagnostics
   !
   ! Output: tend - tendencies: tend_u, tend_w, tend_theta and tend_rho
   !                these are all coupled-variable tendencies.
   !         various other quantities in diag: Smagorinsky eddy viscosity
   !                
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(in) :: state
      type (mpas_pool_type), intent(in) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: nVertLevels              ! for allocating stack variables
      integer, intent(in) :: rk_step
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, iq, iAdvCell
      real (kind=RKIND) :: flux, workpv
      real (kind=RKIND) :: edge_sign, pr_scale, r_dc, r_dv, u_mix_scale

      integer, pointer :: nCells, nEdges, nVertices, nCellsSolve, nEdgesSolve, vertexDegree
      integer, pointer :: moist_start, moist_end
      real (kind=RKIND) :: h_mom_eddy_visc4
      real (kind=RKIND) :: h_theta_eddy_visc4
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), pointer ::  fEdge, dvEdge, dcEdge, areaCell, areaTriangle, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, &
                                                    divergence, vorticity, ke, pv_edge, theta_m, rw, tend_rho, &
                                                    rt_diabatic_tend, tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zxu, cqu, & 
                                                    h_divergence, kdiff, edgesOnCell_sign, edgesOnVertex_sign, rw_save, ru_save

      real (kind=RKIND), dimension(:,:), pointer :: rr_save

      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      real (kind=RKIND), dimension(:,:), pointer :: tend_u_euler, tend_w_euler, tend_theta_euler

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnEdge, cellsOnCell, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, latEdge, angleEdge, u_init

      real (kind=RKIND), dimension( nVertLevels + 1 ) :: wduz, wdwz, wdtz, dpzx
      real (kind=RKIND), dimension( nVertLevels ) :: ru_edge_w, q
      real (kind=RKIND) :: theta_turb_flux, w_turb_flux, r

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension(:), pointer :: rdzu, rdzw, fzm, fzp, qv_init
      real (kind=RKIND), dimension(:,:), pointer :: t_init 

      real (kind=RKIND), dimension(:,:), pointer :: pzp, pzm

      real (kind=RKIND), pointer :: cf1, cf2, cf3
      real (kind=RKIND) :: prandtl_inv, r_areaCell


      real (kind=RKIND), pointer :: r_earth
      real (kind=RKIND), dimension(:,:), pointer :: ur_cell, vr_cell

      real (kind=RKIND), parameter :: c_s = 0.125
      real (kind=RKIND), dimension( nVertLevels ) :: d_diag, d_off_diag, flux_arr
      real (kind=RKIND), dimension(:,:), pointer :: defc_a, defc_b

      real(kind=RKIND), dimension(:,:), pointer :: tend_w_pgf, tend_w_buoy

      real (kind=RKIND), pointer :: coef_3rd_order
      logical, pointer :: config_mix_full
      character (len=StrKIND), pointer :: config_horiz_mixing
      real (kind=RKIND), pointer :: config_del4u_div_factor
      real (kind=RKIND), pointer :: config_h_theta_eddy_visc4
      real (kind=RKIND), pointer :: config_h_mom_eddy_visc4
      real (kind=RKIND), pointer :: config_visc4_2dsmag
      real (kind=RKIND), pointer :: config_len_disp

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
                ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0


      call mpas_pool_get_config(mesh, 'sphere_radius', r_earth)
      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_mix_full', config_mix_full)
      call mpas_pool_get_config(configs, 'config_horiz_mixing', config_horiz_mixing)
      call mpas_pool_get_config(configs, 'config_del4u_div_factor', config_del4u_div_factor)
      call mpas_pool_get_config(configs, 'config_h_theta_eddy_visc4', config_h_theta_eddy_visc4)
      call mpas_pool_get_config(configs, 'config_h_mom_eddy_visc4', config_h_mom_eddy_visc4)
      call mpas_pool_get_config(configs, 'config_visc4_2dsmag', config_visc4_2dsmag)
      call mpas_pool_get_config(configs, 'config_len_disp', config_len_disp)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(state, 'u', u, 2)
      call mpas_pool_get_array(state, 'w', w, 2)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)

      call mpas_pool_get_array(diag, 'uReconstructZonal', ur_cell)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', vr_cell)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'rho_p_save', rr_save)
      call mpas_pool_get_array(diag, 'v', v)
      call mpas_pool_get_array(diag, 'kdiff', kdiff)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'divergence', divergence)
      call mpas_pool_get_array(diag, 'vorticity', vorticity)
      call mpas_pool_get_array(diag, 'ke', ke)
      call mpas_pool_get_array(diag, 'pv_edge', pv_edge)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
      call mpas_pool_get_array(diag, 'h_divergence', h_divergence)

      call mpas_pool_get_array(mesh, 'pzp', pzp)
      call mpas_pool_get_array(mesh, 'pzm', pzm)
      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnVertex_sign', edgesOnVertex_sign)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'angleEdge', angleEdge)
      call mpas_pool_get_array(mesh, 'defc_a', defc_a)
      call mpas_pool_get_array(mesh, 'defc_b', defc_b)
      call mpas_pool_get_array(mesh, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(mesh, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(mesh, 'u_init', u_init)
      call mpas_pool_get_array(mesh, 't_init', t_init)
      call mpas_pool_get_array(mesh, 'qv_init', qv_init)

      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)

      call mpas_pool_get_array(tend, 'u', tend_u)
      call mpas_pool_get_array(tend, 'theta_m', tend_theta)
      call mpas_pool_get_array(tend, 'w', tend_w)
      call mpas_pool_get_array(tend, 'rho_zz', tend_rho)
      call mpas_pool_get_array(tend, 'rt_diabatic_tend', rt_diabatic_tend)
      call mpas_pool_get_array(tend, 'u_euler', tend_u_euler)
      call mpas_pool_get_array(tend, 'theta_euler', tend_theta_euler)
      call mpas_pool_get_array(tend, 'w_euler', tend_w_euler)
      call mpas_pool_get_array(tend, 'w_pgf', tend_w_pgf)
      call mpas_pool_get_array(tend, 'w_buoy', tend_w_buoy)

      call mpas_pool_get_array(diag, 'cqw', cqw)
      call mpas_pool_get_array(diag, 'cqu', cqu)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)

      call mpas_pool_get_dimension(state, 'moist_start', moist_start)
      call mpas_pool_get_dimension(state, 'moist_end', moist_end)

      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(mesh, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(mesh, 'adv_coefs', adv_coefs)
      call mpas_pool_get_array(mesh, 'adv_coefs_3rd', adv_coefs_3rd)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)


      prandtl_inv = 1.0_RKIND / prandtl


      if (rk_step == 1) then

!         tend_u_euler(1:nVertLevels,edgeStart:edgeEnd) = 0.0

         ! Smagorinsky eddy viscosity, based on horizontal deformation (in this case on model coordinate surfaces).
         ! The integration coefficients were precomputed and stored in defc_a and defc_b

         do iCell = cellStart,cellEnd
            d_diag(1:nVertLevels) = 0.0
            d_off_diag(1:nVertLevels) = 0.0
            do iEdge=1,nEdgesOnCell(iCell)
               do k=1,nVertLevels
                  d_diag(k)     = d_diag(k)     + defc_a(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                - defc_b(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                  d_off_diag(k) = d_off_diag(k) + defc_b(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                + defc_a(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
               end do
            end do
            do k=1, nVertLevels
               ! here is the Smagorinsky formulation, 
               ! followed by imposition of an upper bound on the eddy viscosity
               kdiff(k,iCell) = (c_s * config_len_disp)**2 * sqrt(d_diag(k)**2 + d_off_diag(k)**2)
               kdiff(k,iCell) = min(kdiff(k,iCell),(0.01*config_len_disp**2)/dt)
            end do
         end do

         h_mom_eddy_visc4   = config_visc4_2dsmag * config_len_disp**3
         h_theta_eddy_visc4 = h_mom_eddy_visc4

      end if



      ! tendency for density.
      ! accumulate total water here for later use in w tendency calculation.


      ! accumulate horizontal mass-flux

      do iCell=cellStart,cellEnd
         h_divergence(1:nVertLevels,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            edge_sign = edgesOnCell_sign(i,iCell) * dvEdge(iEdge)
            do k=1,nVertLevels
               h_divergence(k,iCell) = h_divergence(k,iCell) + edge_sign * ru(k,iEdge)
            end do
         end do
      end do

      ! compute horiontal mass-flux divergence, add vertical mass flux divergence to complete tend_rho

      do iCell = cellStart,cellEnd
!         qtot(1:nVertLevels,iCell) = 0.0
         r = 1.0 / areaCell(iCell)
         do k = 1,nVertLevels
            h_divergence(k,iCell) = h_divergence(k,iCell) * r

            ! qtot should be moved somewhere else - it is fixed over the dynamics steps
            ! but we are touching the scalars needlessly every acoustic step 

!            do iq = moist_start, moist_end
!               qtot(k,iCell) = qtot(k,iCell) + scalars(iq, k, iCell)
!            end do

         end do
      end do    

      !
      ! dp / dz and tend_rho
      !
      ! only needed on first rk_step with pert variables defined a pert from time t
      !
      if(rk_step == 1) then
        do iCell = cellStart,cellEnd
           k = 1
             dpdz(k,iCell) = 0.5*(  pzm(k,iCell)*(pp(k+1,iCell) - pp(k,iCell))   &
                                  + pzp(k,iCell)*(pp(k+2,iCell) - pp(k,iCell)) )
           do k = 2, nVertLevels-1
              dpdz(k,iCell) = 0.5*(  pzp(k,iCell)*(pp(k+1,iCell) - pp(k  ,iCell))  &
                                   + pzm(k,iCell)*(pp(k  ,iCell) - pp(k-1,iCell)) )
           end do
           k = nVertLevels
           dpdz(k,iCell) = 0.5*pzm(k,iCell)*(pp(k,iCell) - pp(k-1,iCell))

          do k = 1,nVertLevels
            tend_rho(k,iCell) = -h_divergence(k,iCell)-rdzw(k)*(rw(k+1,iCell)-rw(k,iCell))
          end do
        end do
      end if

!$OMP BARRIER

      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !

      do iEdge=edgeSolveStart,edgeSolveEnd

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         ! horizontal pressure gradient 

      if(rk_step == 1) then
!  first push to tend_u_euler
         do k=1,nVertLevels
            tend_u_euler(k,iEdge) =  - cqu(k,iEdge)*((pp(k,cell2)-pp(k,cell1))/dcEdge(iEdge)   &
                                        - zxu(k,iEdge)*(dpdz(k,cell1)+dpdz(k,cell2)) ) / (.5*(zz(k,cell2)+zz(k,cell1)))
         end do
      end if

         ! vertical transport of u

         wduz(1) = 0.

         k = 2
         wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2))*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))
         do k=3,nVertLevels-1
            wduz(k) = flux3( u(k-2,iEdge),u(k-1,iEdge),u(k,iEdge),u(k+1,iEdge),0.5*(rw(k,cell1)+rw(k,cell2)), 1.0_RKIND )
         end do
         k = nVertLevels
         wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2))*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))

         wduz(nVertLevels+1) = 0.

         do k=1,nVertLevels
!            tend_u(k,iEdge) = tend_u(k,iEdge) - rdzw(k)*(wduz(k+1)-wduz(k)) 
            tend_u(k,iEdge) = - rdzw(k)*(wduz(k+1)-wduz(k)) !  first use of tend_u
         end do

         ! Next, nonlinear Coriolis term (q) following Ringler et al JCP 2009

         q(:) = 0.0
         do j = 1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(j,iEdge)
            do k=1,nVertLevels
               workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
               q(k) = q(k) + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * rho_edge(k,eoe)
            end do
         end do

         do k=1,nVertLevels

            ! horizontal ke gradient and vorticity terms in the vector invariant formulation
            ! of the horizontal momentum equation
            tend_u(k,iEdge) = tend_u(k,iEdge) + rho_edge(k,iEdge)* (q(k) - (ke(k,cell2) - ke(k,cell1))       &
                                                                 / dcEdge(iEdge))                            &
                                             - u(k,iEdge)*0.5*(h_divergence(k,cell1)+h_divergence(k,cell2)) 
#ifdef CURVATURE
            ! curvature terms for the sphere
            tend_u(k,iEdge) = tend_u(k,iEdge) &
                             - 2.*omega*cos(angleEdge(iEdge))*cos(latEdge(iEdge))  &
                               *rho_edge(k,iEdge)*.25*(w(k,cell1)+w(k+1,cell1)+w(k,cell2)+w(k+1,cell2))          & 
                             - u(k,iEdge)*.25*(w(k+1,cell1)+w(k,cell1)+w(k,cell2)+w(k+1,cell2))                  &
                               *rho_edge(k,iEdge)/r_earth
#endif
         end do
      end do


      !
      !  horizontal mixing for u
      !  mixing terms are integrated using forward-Euler, so this tendency is only computed in the
      !  first Runge-Kutta substep and saved for use in later RK substeps 2 and 3.
      !

      if (rk_step == 1) then

!  dp/dx now already in tend_u_euler
!
!         tend_u_euler(1:nVertLevels,edgeStart:edgeEnd) = 0.0

!$OMP BARRIER

!         do iEdge=edgeSolveStart,edgeSolveEnd
!            cell1 = cellsOnEdge(1,iEdge)
!            cell2 = cellsOnEdge(2,iEdge)
!            vertex1 = verticesOnEdge(1,iEdge)
!            vertex2 = verticesOnEdge(2,iEdge)
!   
!            do k=1,nVertLevels
               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc2 == constant
               !
!               u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
!                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge),0.25*dcEdge(iEdge))
!               u_diffusion = rho_edge(k,iEdge)* 0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * u_diffusion
!               u_diffusion = u_diffusion * meshScalingDel2(iEdge)
!   
!               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + u_diffusion
!            end do
!         end do


         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
         ! First, storage to hold the result from the first del^2 computation.

         delsq_u(1:nVertLevels,edgeStart:edgeEnd) = 0.0

         do iEdge=edgeStart,edgeEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)
            r_dc = 1./dcEdge(iEdge)
            r_dv = 1./max(dvEdge(iEdge), 0.25*dcEdge(iEdge))

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
!               u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
!                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge), 0.25*dcEdge(iEdge))

               u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) * r_dc  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) * r_dv

               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion

               ! include 2nd-orer diffusion here 
               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) &
                                       + rho_edge(k,iEdge)* 0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * u_diffusion * meshScalingDel2(iEdge)

            end do
         end do

!$OMP BARRIER

         do iVertex=vertexStart,vertexEnd
!            delsq_circulation(1:nVertLevels,iVertex) = 0.0
            delsq_vorticity(1:nVertLevels,iVertex) = 0.0
            r = 1.0 / areaTriangle(iVertex)
            do i=1,vertexDegree
               iEdge = edgesOnVertex(i,iVertex)
               edge_sign = r * dcEdge(iEdge) * edgesOnVertex_sign(i,iVertex)
               do k=1,nVertLevels
!                  delsq_circulation(k,iVertex) = delsq_circulation(k,iVertex) + edge_sign * delsq_u(k,iEdge)
                  delsq_vorticity(k,iVertex) = delsq_vorticity(k,iVertex) + edge_sign * delsq_u(k,iEdge)
               end do
            end do
!            r = 1.0 / areaTriangle(iVertex)
!MGD optimize by using delsq_vorticity from the outset
!            do k=1,nVertLevels
!               delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
!            end do
         end do

         do iCell=cellStart,cellEnd
            delsq_divergence(1:nVertLevels,iCell) = 0.0
            r = 1.0 / areaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               edge_sign = r * dvEdge(iEdge) * edgesOnCell_sign(i,iCell)
               do k=1,nVertLevels
                  delsq_divergence(k,iCell) = delsq_divergence(k,iCell) + edge_sign * delsq_u(k,iEdge)
               end do
            end do
!            do k = 1,nVertLevels
!               delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
!            end do
         end do

!$OMP BARRIER

         do iEdge=edgeSolveStart,edgeSolveEnd
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            u_mix_scale = meshScalingDel4(iEdge)*h_mom_eddy_visc4
            r_dc = u_mix_scale*config_del4u_div_factor/dcEdge(iEdge)
            r_dv = u_mix_scale/max(dvEdge(iEdge), 0.25*dcEdge(iEdge))

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
               ! Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor 
               !    relative to the rotational part.  The stability constraint on the divergence component is much less
               !    stringent than the rotational part, and this flexibility may be useful.
               !
!               u_diffusion =  rho_edge(k,iEdge) * ( config_del4u_div_factor * ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
!                           -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) / max(dvEdge(iEdge), 0.25*dcEdge(iEdge)) &
!                                                  )
!
!               u_diffusion = u_diffusion * meshScalingDel4(iEdge)
!               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - h_mom_eddy_visc4 * u_diffusion

               u_diffusion =  rho_edge(k,iEdge) *  ( ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) * r_dc  &
                                                    -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) * r_dv )
               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - u_diffusion


            end do
         end do


      end if ! (rk_step 1 test for computing mixing terms)

!$OMP BARRIER

!  add in mixing for u

      do iEdge=edgeSolveStart,edgeSolveEnd
         do k=1,nVertLevels
            tend_u(k,iEdge) = tend_u(k,iEdge) + tend_u_euler(k,iEdge)
         end do
      end do


!----------- rhs for w


      !
      !  horizontal advection for w
      !

      do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
         tend_w(1:nVertLevels+1,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            edge_sign = edgesOnCell_sign(i,iCell) * dvEdge(iEdge) * 0.5

            do k=2,nVertLevels
               ru_edge_w(k) = fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge)
            end do

            flux_arr(1:nVertLevels) = 0.0

            ! flux_arr stores the value of w at the cell edge used in the horizontal transport

            do j=1,nAdvCellsForEdge(iEdge)
               iAdvCell = advCellsForEdge(j,iEdge)
               do k=2,nVertLevels
                  scalar_weight = adv_coefs(j,iEdge) + coef_3rd_order * sign(1.0_RKIND,ru_edge_w(k)) * adv_coefs_3rd(j,iEdge)
                  flux_arr(k) = flux_arr(k) + scalar_weight * w(k,iAdvCell)
               end do
            end do

            do k=2,nVertLevels
               tend_w(k,iCell) = tend_w(k,iCell) - edgesOnCell_sign(i,iCell) * ru_edge_w(k)*flux_arr(k)
            end do

         end do
      end do

#ifdef CURVATURE
      do iCell = cellSolveStart, cellSolveEnd
         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell) + (rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))*          &
                                      ( (fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))**2.             &
                                       +(fzm(k)*vr_cell(k,iCell)+fzp(k)*vr_cell(k-1,iCell))**2. )/r_earth   &
                                + 2.*omega*cos(latCell(iCell))                                              &
                                       *(fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))                 &
                                       *(rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))

         end do
      end do
#endif


      !
      !  horizontal mixing for w - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      if (rk_step == 1) then

!         do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
!            tend_w_euler(1:nVertLevels+1,iCell) = 0.0
!            do i=1,nEdgesOnCell(iCell)
!               iEdge = edgesOnCell(i,iCell)
!               cell1 = cellsOnEdge(1,iEdge)
!               cell2 = cellsOnEdge(2,iEdge)
!
!               do k=2,nVertLevels
!                  w_turb_flux = 0.25 * (kdiff(k,cell1)+kdiff(k,cell2)+kdiff(k-1,cell1)+kdiff(k-1,cell2))  &
!                                         * (w(k,cell2) - w(k,cell1)) / dcEdge(iEdge)
!                  w_turb_flux = w_turb_flux * meshScalingDel2(iEdge)
!                  flux = 0.5 * dvEdge(iEdge) * (rho_edge(k,iEdge)+rho_edge(k-1,iEdge)) * w_turb_flux
!                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + edgesOnCell_sign(i,iCell) * flux / areaCell(iCell)
!               end do
!           
!            end do
!         end do


!  !OMP BARRIER  why is this openmp barrier here???

         ! del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
         !
         ! First, storage to hold the result from the first del^2 computation.
         !  we copied code from the theta mixing, hence the theta* names.


         delsq_w(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_w_euler(1:nVertLevels+1,iCell) = 0.0
            r_areaCell = 1./areaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)

               edge_sign = 0.5 * r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) / dcEdge(iEdge)

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

              do k=2,nVertLevels

!                  delsq_w(k,iCell) = delsq_w(k,iCell) + 0.5*edge_sign*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))
                  w_turb_flux =  edge_sign*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))
                  delsq_w(k,iCell) = delsq_w(k,iCell) + w_turb_flux
                  w_turb_flux = w_turb_flux * meshScalingDel2(iEdge) * 0.25 * &
                                  (kdiff(k,cell1)+kdiff(k,cell2)+kdiff(k-1,cell1)+kdiff(k-1,cell2))
                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + w_turb_flux
               end do
            end do
         end do

!         do iCell = cellStart,cellEnd
!            r = 1.0 / areaCell(iCell)
!            do k=2,nVertLevels
!               delsq_w(k,iCell) = delsq_w(k,iCell) * r
!            end do
!         end do

!$OMP BARRIER

         do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
!!            tend_theta_euler(1:nVertLevels,iCell) = 0.0  ! why is this here???
            r_areaCell = h_mom_eddy_visc4/areaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)/dcEdge(iEdge)

               do k=2,nVertLevels
                  tend_w_euler(k,iCell) = tend_w_euler(k,iCell) - edge_sign * (delsq_w(k,cell2) - delsq_w(k,cell1))
               end do
           
            end do
         end do

      end if ! horizontal mixing for w computed in first rk_step


! Note for OpenMP parallelization: We could avoid allocating the delsq_w scratch
!   array, and just use the delsq_theta array as was previously done; however,
!   particularly when oversubscribing cores with threads, there is the risk that
!   some threads may reach code further below that re-uses the delsq_theta array, 
!   in which case we would need a barrier somewhere between here and that code 
!   below to ensure correct behavior.


      !
      !  vertical advection, pressure gradient and buoyancy for w
      !

      do iCell=cellSolveStart,cellSolveEnd

         wdwz(1) = 0.0

         k = 2
         wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
         do k=3,nVertLevels-1
            wdwz(k) = flux3( w(k-2,iCell),w(k-1,iCell),w(k,iCell),w(k+1,iCell),0.5*(rw(k,iCell)+rw(k-1,iCell)), 1.0_RKIND )
         end do
         k = nVertLevels
         wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))

         wdwz(nVertLevels+1) = 0.0

      !  Note: next we are also dividing through by the cell area after the horizontal flux divergence

         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell)/areaCell(iCell) -rdzu(k)*(wdwz(k+1)-wdwz(k))
         end do

!  change from rtheta_pp defined as difference from time t.
!  vertical pressure gradient held at time t.  We will pull buoyance term over later

         if(rk_step == 1) then
            do k=2,nVertLevels
              tend_w_euler(k,iCell) = tend_w_euler(k,iCell)- cqw(k,iCell)*rdzu(k)*(pp(k,iCell)-pp(k-1,iCell)) &
                                  - cqw(k,iCell)*(        &
                                    gravity*  &
                                   ( fzm(k)*(rb(k,iCell)*(qtot(k,iCell)) +         &
                                             rr_save(k,iCell)*(1.+qtot(k,iCell)))                  &
                                    +fzp(k)*(rb(k-1,iCell)*(qtot(k-1,iCell))  +  &
                                             rr_save(k-1,iCell)*(1.+qtot(k-1,iCell))) ))
            end do
          end if

      end do

! add in mixing terms for w

      do iCell = cellSolveStart,cellSolveEnd
         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell) + tend_w_euler(k,iCell)
         end do
      end do



!----------- rhs for theta


      !
      !  horizontal advection for theta
      !

      do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
         tend_theta(1:nVertLevels,iCell) = 0.0
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)

            flux_arr(1:nVertLevels) = 0.0

            do j=1,nAdvCellsForEdge(iEdge)
               iAdvCell = advCellsForEdge(j,iEdge)
               do k=1,nVertLevels
                  scalar_weight = adv_coefs(j,iEdge) + coef_3rd_order*sign(1.0_RKIND,ru(k,iEdge))*adv_coefs_3rd(j,iEdge)
                  flux_arr(k) = flux_arr(k) + scalar_weight* theta_m(k,iAdvCell)
               end do
            end do

            do k=1,nVertLevels
               tend_theta(k,iCell) = tend_theta(k,iCell) - edgesOnCell_sign(i,iCell) * ru(k,iEdge) * flux_arr(k)
            end do

         end do
      end do

!  addition to pick up perturbation flux for rtheta_pp equation

      if(rk_step > 1) then
        do iCell=cellSolveStart,cellSolveEnd
          do i=1,nEdgesOnCell(iCell) 
            iEdge = edgesOnCell(i,iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
               flux = edgesOnCell_sign(i,iCell)*dvEdge(iEdge)*(ru_save(k,iEdge)-ru(k,iEdge))
               flux = flux*0.5*(theta_m(k,cell2)+theta_m(k,cell1))
!               tend_theta(k,iCell) = tend_theta(k,iCell)-flux/areaCell(iCell)
               tend_theta(k,iCell) = tend_theta(k,iCell)-flux  ! division by areaCell picked up down below
            end do
          end do
        end do
      end if

      !
      !  horizontal mixing for theta_m - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      if (rk_step == 1) then

! Smagorinsky 2nd order horizontal mixing
!         do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
!            do i=1,nEdgesOnCell(iCell)
!               iEdge = edgesOnCell(i,iCell)
!               cell1 = cellsOnEdge(1,iEdge)
!               cell2 = cellsOnEdge(2,iEdge)
!
!               do k=1,nVertLevels
!                  theta_turb_flux = 0.5*(kdiff(k,cell1)+kdiff(k,cell2))*prandtl_inv  &
!                                        *(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
!                  theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
!                  flux = dvEdge(iEdge) * rho_edge(k,iEdge) * theta_turb_flux
!                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + edgesOnCell_sign(i,iCell) * flux / areaCell(iCell)
!               end do
!
!            end do
!         end do

         delsq_theta(1:nVertLevels,cellStart:cellEnd) = 0.0

         do iCell=cellStart,cellEnd
            tend_theta_euler(1:nVertLevels,iCell) = 0.0
            r_areaCell = 1./areaCell(iCell)
            do i=1,nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i,iCell)
               edge_sign = r_areaCell*edgesOnCell_sign(i,iCell) * dvEdge(iEdge) / dcEdge(iEdge)
               pr_scale = prandtl_inv * meshScalingDel2(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k=1,nVertLevels
!!                  delsq_theta(k,iCell) = delsq_theta(k,iCell) + edge_sign * rho_edge(k,iEdge) * (theta_m(k,cell2) - theta_m(k,cell1))

!  we are computing the Smagorinsky filter at more points than needed here

                  theta_turb_flux = edge_sign*(theta_m(k,cell2) - theta_m(k,cell1))*rho_edge(k,iEdge)
                  delsq_theta(k,iCell) = delsq_theta(k,iCell) + theta_turb_flux
                  theta_turb_flux = theta_turb_flux*0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * pr_scale
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + theta_turb_flux

!!                  theta_turb_flux = 0.5*(kdiff(k,cell1)+kdiff(k,cell2))*prandtl_inv & 
!!                                        *(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
!!                  theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
!!                  flux = dvEdge(iEdge) * rho_edge(k,iEdge) * theta_turb_flux
!!                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + edgesOnCell_sign(i,iCell) * flux / areaCell(iCell)

               end do
            end do
!!            r = 1.0 / areaCell(iCell)
!!            do k=1,nVertLevels
!!               delsq_theta(k,iCell) = delsq_theta(k,iCell) * r
!!            end do
         end do

!$OMP BARRIER

         do iCell=cellSolveStart,cellSolveEnd    ! Technically updating fewer cells than before...
            r_areaCell = h_theta_eddy_visc4*prandtl_inv/areaCell(iCell)
            do i=1,nEdgesOnCell(iCell)

               iEdge = edgesOnCell(i,iCell)
               edge_sign = meshScalingDel4(iEdge)*r_areaCell*dvEdge(iEdge)*edgesOnCell_sign(i,iCell)/dcEdge(iEdge)

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               do k=1,nVertLevels
                  tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) - edge_sign*(delsq_theta(k,cell2) - delsq_theta(k,cell1))
               end do

            end do
         end do

      end if ! theta mixing calculated first rk_step


      !
      !  vertical advection plus diabatic term
      !  Note: we are also dividing through by the cell area after the horizontal flux divergence
      !
      do iCell = cellSolveStart,cellSolveEnd

         wdtz(1) = 0.0

         k = 2
!!         wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))  ! rtheta_pp redefinition
         wdtz(k) =  rw_save(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
         do k=3,nVertLevels-1
            wdtz(k) = flux3( theta_m(k-2,iCell),theta_m(k-1,iCell),theta_m(k,iCell),theta_m(k+1,iCell), rw(k,iCell), coef_3rd_order )
            wdtz(k) =  wdtz(k) + (rw_save(k,icell)-rw(k,iCell))*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))  ! rtheta_pp redefinition
         end do
         k = nVertLevels
!!         wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
         wdtz(k) =  rw_save(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))  ! rtheta_pp redefinition

         wdtz(nVertLevels+1) = 0.0

         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell)/areaCell(iCell) -rdzw(k)*(wdtz(k+1)-wdtz(k))
            tend_theta(k,iCell) = tend_theta(k,iCell) + rho_zz(k,iCell)*rt_diabatic_tend(k,iCell)
         end do
      end do


      do iCell = cellSolveStart,cellSolveEnd
         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell) + tend_theta_euler(k,iCell)
         end do
      end do

   end subroutine atm_compute_dyn_tend
