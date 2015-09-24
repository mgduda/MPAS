   subroutine atm_recover_large_step_variables( state, diag, tend, mesh, configs, dt, ns, rk_step, &
                                       cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd, &
                                       cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd)

      ! reconstitute state variables from acoustic-step perturbation variables 
      ! after the acoustic steps.  The perturbation variables were originally set in
      ! subroutine atm_set_smlstep_pert_variables prior to their acoustic-steps update.
      ! we are also computing a few other state-derived variables here.

      implicit none

      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(inout) :: tend
      type (mpas_pool_type), intent(inout) :: mesh
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: ns, rk_step
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: cellStart, cellEnd, vertexStart, vertexEnd, edgeStart, edgeEnd
      integer, intent(in) :: cellSolveStart, cellSolveEnd, vertexSolveStart, vertexSolveEnd, edgeSolveStart, edgeSolveEnd


      real (kind=RKIND), dimension(:,:), pointer :: wwAvg, rw_save, w, rw, rw_p, rtheta_p, rtheta_pp,   &
                                                    rtheta_p_save, rt_diabatic_tend, rho_p, rho_p_save, &
                                                    rho_pp, rho_zz, rho_base, ruAvg, ru_save, ru_p, u, ru, &
                                                    exner, exner_base, rtheta_base, pressure_p,         &
                                                    zz, theta_m, pressure_b
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, dvEdge, areaCell
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3, zb_cell, zb3_cell
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      integer, dimension(:), pointer :: nEdgesOnCell

      integer :: i, iCell, iEdge, k, cell1, cell2
      integer, pointer :: nVertLevels, nCells, nCellsSolve, nEdges, nEdgesSolve
      real (kind=RKIND) :: rcv, p0, flux
      real (kind=RKIND), pointer :: cf1, cf2, cf3, coef_3rd_order


      call mpas_pool_get_array(diag, 'wwAvg', wwAvg)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_p', rw_p)
      call mpas_pool_get_array(state, 'w', w, 2)

      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(diag, 'rtheta_p_save', rtheta_p_save)
      call mpas_pool_get_array(diag, 'rtheta_pp', rtheta_pp)
      call mpas_pool_get_array(diag, 'rtheta_base', rtheta_base)
      call mpas_pool_get_array(tend, 'rt_diabatic_tend', rt_diabatic_tend)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 2)
      call mpas_pool_get_array(state, 'scalars', scalars, 2)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 2)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rho_p_save', rho_p_save)
      call mpas_pool_get_array(diag, 'rho_pp', rho_pp)
      call mpas_pool_get_array(diag, 'rho_base', rho_base)

      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'ru_p', ru_p)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(state, 'u', u, 2)

      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'exner_base', exner_base)

      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'pressure_base', pressure_b)

      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)
      call mpas_pool_get_array(mesh, 'zb_cell', zb_cell)
      call mpas_pool_get_array(mesh, 'zb3_cell', zb3_cell)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

      call mpas_pool_get_config(configs, 'config_coef_3rd_order', coef_3rd_order)

      rcv = rgas/(cp-rgas)
      p0 = 1.0e+05  ! this should come from somewhere else...


      ! compute new density everywhere so we can compute u from ru.
      ! we will also need it to compute theta_m below

      do iCell=cellStart,cellEnd

         do k = 1, nVertLevels
            rho_p(k,iCell) = rho_p_save(k,iCell) + rho_pp(k,iCell)

            rho_zz(k,iCell) = rho_p(k,iCell) + rho_base(k,iCell)
         end do

         w(1,iCell) = 0.0

         do k = 2, nVertLevels
            wwAvg(k,iCell) = rw_save(k,iCell) + (wwAvg(k,iCell) / real(ns,RKIND))
            rw(k,iCell) = rw_save(k,iCell) + rw_p(k,iCell)

          ! pick up part of diagnosed w from omega - divide by density later
            w(k,iCell) = rw(k,iCell)/(fzm(k)*zz(k,iCell)+fzp(k)*zz(k-1,iCell))
                                      
         end do

         w(nVertLevels+1,iCell) = 0.0

         if (rk_step == 3) then
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p_save(k,iCell) + rtheta_pp(k,iCell) &
                                 - dt * rho_zz(k,iCell) * rt_diabatic_tend(k,iCell)
            end do
         else
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p_save(k,iCell) + rtheta_pp(k,iCell)
            end do
         end if

         do k = 1, nVertLevels
            theta_m(k,iCell) = (rtheta_p(k,iCell) + rtheta_base(k,iCell))/rho_zz(k,iCell)
            exner(k,iCell) = (zz(k,iCell)*(rgas/p0)*(rtheta_p(k,iCell)+rtheta_base(k,iCell)))**rcv
            ! pressure_p is perturbation pressure
            pressure_p(k,iCell) = zz(k,iCell) * rgas * (exner(k,iCell)*rtheta_p(k,iCell)+rtheta_base(k,iCell)  &
                                                          * (exner(k,iCell)-exner_base(k,iCell)))
         end do

      end do

      ! recover time-averaged ruAvg on all edges of owned cells (for upcoming scalar transport).  
      ! we solved for these in the acoustic-step loop.  
      ! we will compute ru and u here also, given we are here, even though we only need them on nEdgesSolve

      ! Avoid FP errors caused by a potential division by zero below by 
      ! initializing the "garbage cell" of rho_zz to a non-zero value
      rho_zz(:,nCells+1) = 1.0

!$OMP BARRIER

      do iEdge=edgeStart,edgeEnd

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 1, nVertLevels
            ruAvg(k,iEdge) = ru_save(k,iEdge) + (ruAvg(k,iEdge) / real(ns,RKIND))
            ru(k,iEdge) = ru_save(k,iEdge) + ru_p(k,iEdge)
            u(k,iEdge) = 2.*ru(k,iEdge)/(rho_zz(k,cell1)+rho_zz(k,cell2))
         end do
      end do

!$OMP BARRIER

      do iCell=cellStart,cellEnd

         !  finish recovering w from (rho*omega)_p.  as when we formed (rho*omega)_p from u and w, we need
         !  to use the same flux-divergence operator as is used for the horizontal theta transport
         !  (See Klemp et al 2003).

         do i=1,nEdgesOnCell(iCell)
            iEdge=edgesOnCell(i,iCell)

            flux = (cf1*ru(1,iEdge) + cf2*ru(2,iEdge) + cf3*ru(3,iEdge))
            w(1,iCell) = w(1,iCell) + edgesOnCell_sign(i,iCell) * &
                                   (zb_cell(1,i,iCell) + sign(1.0_RKIND,flux)*coef_3rd_order*zb3_cell(1,i,iCell))*flux

            do k = 2, nVertLevels
               flux = (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
               w(k,iCell) = w(k,iCell) + edgesOnCell_sign(i,iCell) * &
                                    (zb_cell(k,i,iCell)+sign(1.0_RKIND,flux)*coef_3rd_order*zb3_cell(k,i,iCell))*flux
            end do

         end do

         w(1,iCell) = w(1,iCell)/(cf1*rho_zz(1,iCell)+cf2*rho_zz(2,iCell)+cf3*rho_zz(3,iCell))
         do k = 2, nVertLevels
           w(k,iCell) = w(k,iCell)/(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell))
         end do

      end do

   end subroutine atm_recover_large_step_variables
