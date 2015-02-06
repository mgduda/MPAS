! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module init_atm_cases

   use mpas_kind_types
   use mpas_derived_types
   use mpas_constants
   use mpas_dmpar
   use atm_advection
   use mpas_RBF_interpolation
   use mpas_vector_reconstruction
   use mpas_timer

   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in mpas_timekeeping
   use mpas_timekeeping !, only: MPAS_Time_type, MPAS_TimeInterval_type, MPAS_Clock_type, &
                        !        mpas_set_time, mpas_set_timeInterval, mpas_get_time, operator(+), add_t_ti
   


   contains


   subroutine init_atm_setup_case(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the hydrostatic test case
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain


      integer :: i
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: mesh
      type (mpas_pool_type), pointer :: fg
      type (mpas_pool_type), pointer :: state
      type (mpas_pool_type), pointer :: diag
      type (mpas_pool_type), pointer :: diag_physics

      integer, pointer :: config_init_case

      character(len=StrKIND), pointer :: mminlu

      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, pointer :: nVertLevels


      call mpas_pool_get_config(domain % blocklist % configs, 'config_init_case', config_init_case)

      write(0,*) ' Jablonowski and Williamson baroclinic wave test case '
      if (config_init_case == 1) write(0,*) ' no initial perturbation '
      if (config_init_case == 2) write(0,*) ' initial perturbation included '
      if (config_init_case == 3) write(0,*) ' normal-mode perturbation included '
      block_ptr => domain % blocklist
      do while (associated(block_ptr))

         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', mesh)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', state)
         call mpas_pool_get_subpool(block_ptr % structs, 'diag', diag)

         call mpas_pool_get_dimension(block_ptr % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block_ptr % dimensions, 'nVertLevels', nVertLevels)

         write(0,*) ' calling test case setup '
         call init_atm_case_jw(mesh, nCells, nVertLevels, state, diag, block_ptr % configs, config_init_case)
         call decouple_variables(mesh, nCells, nVertLevels, state, diag)
         write(0,*) ' returned from test case setup '
         block_ptr => block_ptr % next
      end do

   end subroutine init_atm_setup_case

!----------------------------------------------------------------------------------------------------------

   subroutine init_atm_case_jw(mesh, nCells, nVertLevels, state, diag, configs, test_case)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup baroclinic wave test case from Jablonowski and Williamson 2008 (QJRMS)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: configs
      integer, intent(in) :: test_case

      real (kind=RKIND), parameter :: u0 = 35.0
      real (kind=RKIND), parameter :: alpha_grid = 0.  ! no grid rotation

!      real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e

      real (kind=RKIND), parameter :: t0b = 250., t0 = 288., delta_t = 4.8e+05, dtdz = 0.005, eta_t = 0.2
      real (kind=RKIND), parameter :: u_perturbation = 1., pert_radius = 0.1, latitude_pert = 40., longitude_pert = 20.
      real (kind=RKIND), parameter :: theta_c = pii/4.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: k_x = 9.           ! Normal mode wave number

      real (kind=RKIND), dimension(:), pointer :: rdzw, dzu, rdzu, fzm, fzp
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:), pointer :: zgrid, zxu, zz, hx
      real (kind=RKIND), dimension(:,:), pointer :: pressure, ppb, pb, rho_zz, rb, rr, tb, rtb, p, pp, dss, t, rt
      real (kind=RKIND), dimension(:,:), pointer :: u, ru, w, rw, v
      real (kind=RKIND), dimension(:,:), pointer :: rho, theta
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      
!.. initialization of moisture:
      integer, pointer :: index_qv
      real (kind=RKIND), parameter :: rh_max = 0.40 ! Maximum relative humidity
!      real (kind=RKIND), parameter :: rh_max = 0.70 ! Maximum relative humidity
      real (kind=RKIND), dimension(nVertLevels, nCells) :: qsat, relhum
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars
!.. end initialization of moisture.

!  scalars for scaling test
      integer, pointer :: index_p0,index_p1,index_p2,index_p3,index_p4,index_p5,index_p6,index_p7,index_p8,index_p9
      real (kind=RKIND) :: stmp

      integer :: iCell, iCell1, iCell2 , iEdge, vtx1, vtx2, ivtx, i, k, nz, itr, itrp, cell1, cell2
      integer, pointer :: nz1, nCellsSolve, nEdges, maxEdges, nVertices

      !This is temporary variable here. It just need when calculate tangential velocity v.
      integer :: eoe, j
      integer, dimension(:), pointer :: nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnEdge, cellsOnEdge, verticesOnEdge, cellsOnCell
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      real (kind=RKIND) :: flux, fluxk, lat1, lat2, eta_v, r_pert, u_pert, lat_pert, lon_pert, r

      real (kind=RKIND) :: ptop, p0, phi
      real (kind=RKIND) :: lon_Edge

      real (kind=RKIND) :: r_earth, etavs, ztemp, zd, zt, dz, gam, delt, str

      real (kind=RKIND) :: es, qvs, xnutr, znut, ptemp 
      integer :: iter

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: hyai, hybi, znu, znw, znwc, znwv, hyam, hybm
      real (kind=RKIND), dimension(nVertLevels + 1 ) :: znuc, znuv, bn, divh, dpn

      real (kind=RKIND), dimension(nVertLevels + 1 ) :: sh, zw, ah
      real (kind=RKIND), dimension(nVertLevels ) :: zu, dzw, rdzwp, rdzwm
      real (kind=RKIND), dimension(nVertLevels ) :: eta, etav, teta, ppi, tt, temperature_1d

      real (kind=RKIND) :: d1, d2, d3, cof1, cof2, psurf
      real (kind=RKIND), pointer :: cf1, cf2, cf3

      !  storage for (lat,z) arrays for zonal velocity calculation

      logical, parameter :: rebalance = .true.
      integer, parameter :: nlat=721
      real (kind=RKIND), dimension(nVertLevels) :: flux_zonal
      real (kind=RKIND), dimension(nVertLevels + 1, nlat) :: zgrid_2d
      real (kind=RKIND), dimension(nVertLevels, nlat) :: u_2d, pp_2d, rho_2d, qv_2d, etavs_2d, zz_2d
      real (kind=RKIND), dimension(nVertLevels, nlat) :: p_2d, pb_2d, ppb_2d, rr_2d, rb_2d, tb_2d, rtb_2d
      real (kind=RKIND), dimension(nVertLevels, nlat-1) :: zx_2d 
      real (kind=RKIND), dimension(nlat) :: lat_2d
      real (kind=RKIND) :: dlat, hx_1d
      real (kind=RKIND) :: z_edge, z_edge3, d2fdx2_cell1, d2fdx2_cell2

!      logical, parameter :: moisture = .true.
      logical, parameter :: moisture = .false.

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:), pointer :: xEdge, yEdge, zEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, latVertex, lonVertex, latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: fEdge, fVertex

      real (kind=RKIND), pointer :: sphere_radius
      logical, pointer :: on_a_sphere
      real (kind=RKIND), pointer :: config_coef_3rd_order
      integer, pointer :: config_theta_adv_order
      integer, pointer :: config_init_case


      call mpas_pool_get_config(configs, 'config_init_case', config_init_case)
      call mpas_pool_get_config(configs, 'config_coef_3rd_order', config_coef_3rd_order)
      call mpas_pool_get_config(configs, 'config_theta_adv_order', config_theta_adv_order)


      !
      ! Scale all distances and areas from a unit sphere to one with radius sphere_radius
      !
      call mpas_pool_get_array(mesh, 'xCell', xCell)
      call mpas_pool_get_array(mesh, 'yCell', yCell)
      call mpas_pool_get_array(mesh, 'zCell', zCell)
      call mpas_pool_get_array(mesh, 'xEdge', xEdge)
      call mpas_pool_get_array(mesh, 'yEdge', yEdge)
      call mpas_pool_get_array(mesh, 'zEdge', zEdge)
      call mpas_pool_get_array(mesh, 'xVertex', xVertex)
      call mpas_pool_get_array(mesh, 'yVertex', yVertex)
      call mpas_pool_get_array(mesh, 'zVertex', zVertex)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_config(mesh, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(mesh, 'sphere_radius', sphere_radius)

      xCell(:) = xCell(:) * sphere_radius
      yCell(:) = yCell(:) * sphere_radius
      zCell(:) = zCell(:) * sphere_radius
      xVertex(:) = xVertex(:) * sphere_radius
      yVertex(:) = yVertex(:) * sphere_radius
      zVertex(:) = zVertex(:) * sphere_radius
      xEdge(:) = xEdge(:) * sphere_radius
      yEdge(:) = yEdge(:) * sphere_radius
      zEdge(:) = zEdge(:) * sphere_radius
      dvEdge(:) = dvEdge(:) * sphere_radius
      dcEdge(:) = dcEdge(:) * sphere_radius
      areaCell(:) = areaCell(:) * sphere_radius**2.0
      areaTriangle(:) = areaTriangle(:) * sphere_radius**2.0
      kiteAreasOnVertex(:,:) = kiteAreasOnVertex(:,:) * sphere_radius**2.0

      call mpas_pool_get_array(mesh, 'weightsOnEdge', weightsOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnEdge', nEdgesOnEdge)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnEdge', edgesOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'deriv_two', deriv_two)
      call mpas_pool_get_array(mesh, 'zb', zb)
      call mpas_pool_get_array(mesh, 'zb3', zb3)

      call mpas_pool_get_dimension(mesh, 'nVertLevels', nz1)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
      nz = nz1 + 1

      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'rdzu', rdzu)
      call mpas_pool_get_array(mesh, 'fzm', fzm)
      call mpas_pool_get_array(mesh, 'fzp', fzp)
      call mpas_pool_get_array(mesh, 'zxu', zxu)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(mesh, 'hx', hx)
      call mpas_pool_get_array(mesh, 'dss', dss)

      call mpas_pool_get_array(diag, 'exner_base', pb)
      call mpas_pool_get_array(diag, 'rho_base', rb)
      call mpas_pool_get_array(diag, 'rho_p', rr)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta_base', tb)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rtheta_base', rtb)
      call mpas_pool_get_array(diag, 'rtheta_p', rt)
      call mpas_pool_get_array(diag, 'exner', p)
      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'v', v)

      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
      call mpas_pool_get_array(state, 'theta_m', t, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(state, 'u', u, 1)
      call mpas_pool_get_array(state, 'w', w, 1)

      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'latVertex', latVertex)
      call mpas_pool_get_array(mesh, 'lonVertex', lonVertex)
      call mpas_pool_get_array(mesh, 'latEdge', latEdge)
      call mpas_pool_get_array(mesh, 'lonEdge', lonEdge)
      call mpas_pool_get_array(mesh, 'fEdge', fEdge)
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)

      call mpas_pool_get_array(mesh, 'cf1', cf1)
      call mpas_pool_get_array(mesh, 'cf2', cf2)
      call mpas_pool_get_array(mesh, 'cf3', cf3)

!.. initialization of moisture:
      scalars(:,:,:) = 0.0
      qsat(:,:)      = 0.0
      relhum(:,:)    = 0.0
      qv_2d(:,:)     = 0.0
!.. end initialization of moisture.

!  initialization of passive scalars for scaling tests
!  WCS 5 dec 2014
      call mpas_pool_get_dimension(state, 'index_p0', index_p0) 
      call mpas_pool_get_dimension(state, 'index_p1', index_p1) 
      call mpas_pool_get_dimension(state, 'index_p2', index_p2) 
      call mpas_pool_get_dimension(state, 'index_p3', index_p3) 
      call mpas_pool_get_dimension(state, 'index_p4', index_p4) 
      call mpas_pool_get_dimension(state, 'index_p5', index_p5) 
      call mpas_pool_get_dimension(state, 'index_p6', index_p6) 
      call mpas_pool_get_dimension(state, 'index_p7', index_p7) 
      call mpas_pool_get_dimension(state, 'index_p8', index_p8) 
      call mpas_pool_get_dimension(state, 'index_p9', index_p9) 
      do iCell=1,nCells
        stmp = max(0.,sin(9.*latCell(iCell))*sin(9.*lonCell(iCell)))
        do k=1,nz1
          scalars(index_p0,k,iCell) = stmp
          scalars(index_p1,k,iCell) = stmp
          scalars(index_p2,k,iCell) = stmp
          scalars(index_p3,k,iCell) = stmp
          scalars(index_p4,k,iCell) = stmp
          scalars(index_p5,k,iCell) = stmp
          scalars(index_p6,k,iCell) = stmp
          scalars(index_p7,k,iCell) = stmp
          scalars(index_p8,k,iCell) = stmp
          scalars(index_p9,k,iCell) = stmp
        end do
      end do

      surface_pressure(:) = 0.0

      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv) 

      xnutr = 0.0
      zd = 12000.0
      znut = eta_t

      etavs = (1.0 - 0.252) * pii/2.
      r_earth = sphere_radius
      omega_e = omega
      p0 = 1.0e+05

      write(0,*) ' point 1 in test case setup '

! We may pass in an hx(:,:) that has been precomputed elsewhere.
! For now it is independent of k

      do iCell=1,nCells
        do k=1,nz
          phi = latCell(iCell)
          hx(k,iCell) = u0/gravity*cos(etavs)**1.5                                   &
                      *((-2.*sin(phi)**6                                   &
                            *(cos(phi)**2+1./3.)+10./63.)                  &
                            *(u0)*cos(etavs)**1.5                          &
                       +(1.6*cos(phi)**3                                   &
                            *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)
        end do
      end do

      !     Metrics for hybrid coordinate and vertical stretching

      str = 1.5
      zt = 45000.
      dz = zt/float(nz1)

      write(0,*) ' hx computation complete '

      do k=1,nz
		
!           sh(k) is the stretching specified for height surfaces

            sh(k) = (real(k-1)*dz/zt)**str 
				
!           to specify specific heights zc(k) for coordinate surfaces,
!           input zc(k) and define sh(k) = zc(k)/zt
!           zw(k) is the hieght of zeta surfaces
!                zw(k) = (k-1)*dz yields constant dzeta
!                        and nonconstant dzeta/dz
!                zw(k) = sh(k)*zt yields nonconstant dzeta
!                        and nearly constant dzeta/dz 

            zw(k) = float(k-1)*dz
!            zw(k) = sh(k)*zt
!
!           ah(k) governs the transition between terrain-following 
!           and pureheight coordinates
!                ah(k) = 0 is a terrain-following coordinate
!                ah(k) = 1 is a height coordinate
 
            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
!            ah(k) = 0.
	    write(0,*) ' k, sh, zw, ah ',k,sh(k),zw(k),ah(k)
      end do
      do k=1,nz1
         dzw (k) = zw(k+1)-zw(k)
         rdzw(k) = 1./dzw(k)
         zu(k  ) = .5*(zw(k)+zw(k+1))
      end do
      do k=2,nz1
         dzu (k)  = .5*(dzw(k)+dzw(k-1))
         rdzu(k)  =  1./dzu(k)
         fzp (k)  = .5* dzw(k  )/dzu(k)
         fzm (k)  = .5* dzw(k-1)/dzu(k)
         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
      end do

!**********  how are we storing cf1, cf2 and cf3?

      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2) 
      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3) 
      CF1  = FZP(2) + COF1
      CF2  = FZM(2) - COF1 - COF2
      CF3  = COF2       

!      d1  = .5*dzw(1)
!      d2  = dzw(1)+.5*dzw(2)
!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))

      write(0,*) ' cf1, cf2, cf3 = ',cf1,cf2,cf3

      do iCell=1,nCells
        do k=1,nz	
          zgrid(k,iCell) = (1.-ah(k))*(sh(k)*(zt-hx(k,iCell))+hx(k,iCell))  &
                         + ah(k) * sh(k)* zt	
        end do
        do k=1,nz1
          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
        end do
      end do

      do i=1, nEdges
        iCell1 = cellsOnEdge(1,i)
        iCell2 = cellsOnEdge(2,i)
        do k=1,nz
          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
        end do
      end do
      do i=1, nCells
        do k=1,nz1
          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
          dss(k,i) = 0.
          ztemp = zgrid(k,i)
          if(ztemp.gt.zd+.1)  then
             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
          end if
        end do
      end do

      !do k=1,nz1
      !  write(0,*) ' k, zgrid(k,1),hx(k,1) ',k,zgrid(k,1),hx(k,1)
      !end do

      !do k=1,nz1
      !  write(0,*) ' k, zx(k,1) ',k,zx(k,1)
      !end do

      write(0,*) ' grid metrics setup complete '

!**************  section for 2d (z,lat) calc for zonal velocity

      dlat = 0.5*pii/float(nlat-1)
      do i = 1,nlat

        lat_2d(i) = float(i-1)*dlat
        phi = lat_2d(i)
        hx_1d    = u0/gravity*cos(etavs)**1.5                           &
                   *((-2.*sin(phi)**6                                   &
                         *(cos(phi)**2+1./3.)+10./63.)                  &
                         *(u0)*cos(etavs)**1.5                          &
                    +(1.6*cos(phi)**3                                   &
                         *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)

        do k=1,nz	
          zgrid_2d(k,i) = (1.-ah(k))*(sh(k)*(zt-hx_1d)+hx_1d)  &
                         + ah(k) * sh(k)* zt	
        end do
        do k=1,nz1
          zz_2d(k,i) = (zw(k+1)-zw(k))/(zgrid_2d(k+1,i)-zgrid_2d(k,i))
        end do

        do k=1,nz1
          ztemp    = .5*(zgrid_2d(k+1,i)+zgrid_2d(k,i))
          ppb_2d(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb_2d(k,i) = (ppb_2d(k,i)/p0)**(rgas/cp)
          rb_2d(k,i) = ppb_2d(k,i)/(rgas*t0b*zz_2d(k,i))
          tb_2d(k,i) = t0b/pb_2d(k,i)
          rtb_2d(k,i) = rb_2d(k,i)*tb_2d(k,i)
          p_2d(k,i) = pb_2d(k,i)
          pp_2d(k,i) = 0.0
          rr_2d(k,i) = 0.0
        end do


        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb_2d(k,i)+pp_2d(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k) >= znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do

          phi = lat_2d(i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*qv_2d(k,i))


            ztemp   = .5*(zgrid_2d(k,i)+zgrid_2d(k+1,i))
            ptemp   = ppb_2d(k,i) + pp_2d(k,i)

            !get moisture 
            if (moisture) then
              qv_2d(k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*qv_2d(k,i))
          end do

          do itrp = 1,25
            do k=1,nz1				
              rr_2d(k,i)  = (pp_2d(k,i)/(rgas*zz_2d(k,i)) - rb_2d(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                            &
                          *(1.25*(rr_2d(1,i)+rb_2d(1,i))*(1.+qv_2d(1,i))   &
                            -.25*(rr_2d(2,i)+rb_2d(2,i))*(1.+qv_2d(2,i)))

            ppi(1) = ppi(1)-ppb_2d(1,i)
            do k=1,nz1-1

              ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                       &
                            ( (rr_2d(k  ,i)+(rr_2d(k  ,i)+rb_2d(k  ,i))*qv_2d(k  ,i))*fzp(k+1)   &
                            + (rr_2d(k+1,i)+(rr_2d(k+1,i)+rb_2d(k+1,i))*qv_2d(k+1,i))*fzm(k+1) )
            end do

            do k=1,nz1
              pp_2d(k,i) = .2*ppi(k)+.8*pp_2d(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1
          rho_2d(k,i) = rr_2d(k,i)+rb_2d(k,i)
          etavs_2d(k,i) = ((ppb_2d(k,i)+pp_2d(k,i))/p0 - 0.252)*pii/2.
          u_2d(k,i) = u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(k,i))**1.5)
        end do

      end do  ! end loop over latitudes for 2D zonal wind field calc

      !SHP-balance:: in case of rebalacing for geostrophic wind component
      if (rebalance) then

        do i=1,nlat-1
          do k=1,nz1
            zx_2d(k,i) = (zgrid_2d(k,i+1)-zgrid_2d(k,i))/(dlat*r_earth)
          end do
        end do

        call init_atm_recompute_geostrophic_wind(u_2d, rho_2d, pp_2d, qv_2d, lat_2d, zz_2d, zx_2d,     &
                                        cf1, cf2, cf3, fzm, fzp, rdzw, nz1, nlat, dlat, sphere_radius)

      end if

!******************************************************************      

!
!---- baroclinc wave initialization ---------------------------------
!
!     reference sounding based on dry isothermal atmosphere
!
      do i=1, nCells
        do k=1,nz1
          ztemp    = .5*(zgrid(k+1,i)+zgrid(k,i))
          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
          rb (k,i) = ppb(k,i)/(rgas*t0b*zz(k,i))
          tb (k,i) = t0b/pb(k,i)
          rtb(k,i) = rb(k,i)*tb(k,i)
          p  (k,i) = pb(k,i)
          pp (k,i) = 0.
          rr (k,i) = 0.
        end do

!       if(i == 1) then
!         do k=1,nz1
!           write(0,*) ' k, ppb, pb, rb, tb (k,1) ',k,ppb(k,1),pb(k,1),rb(k,1)*zz(k,1),tb(k,1)
!         end do
!       end if

      200 format(4i6,8(1x,e15.8))
      201 format(3i6,8(1x,e15.8))
      202 format(2i6,10(1x,e15.8))
      203 format(i6,10(1x,e15.8))

!     iterations to converge temperature as a function of pressure
!
        do itr = 1,10

          do k=1,nz1
            eta (k) = (ppb(k,i)+pp(k,i))/p0
            etav(k) = (eta(k)-.252)*pii/2.
            if(eta(k).ge.znut)  then
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
            else
              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
            end if
          end do
          phi = latCell(i)
          do k=1,nz1
            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
                            *sqrt(cos(etav(k)))*                   &
                              ((-2.*sin(phi)**6                    &
                                   *(cos(phi)**2+1./3.)+10./63.)   &
                                   *2.*u0*cos(etav(k))**1.5        &
                              +(1.6*cos(phi)**3                    &
                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*scalars(index_qv,k,i))

            ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
            ptemp   = ppb(k,i) + pp(k,i)

            !get moisture 
            if (moisture) then
 
                !scalars(index_qv,k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )

               if(ptemp < 50000.) then
                  relhum(k,i) = 0.0
               elseif(ptemp > p0) then
                  relhum(k,i) = 1.0
               else
                  relhum(k,i) = (1.-((p0-ptemp)/50000.)**1.25)
               end if
               relhum(k,i) = min(rh_max,relhum(k,i))

               !.. calculation of water vapor mixing ratio:
               if (temperature_1d(k) > 273.15) then
                   es  = 1000.*0.6112*exp(17.67*(temperature_1d(k)-273.15)/(temperature_1d(k)-29.65))
               else
                   es  = 1000.*0.6112*exp(21.8745584*(temperature_1d(k)-273.15)/(temperature_1d(k)-7.66))
               end if
               qsat(k,i) = (287.04/461.6)*es/(ptemp-es)
               if(relhum(k,i) .eq. 0.0) qsat(k,i) = 0.0
               scalars(index_qv,k,i) = relhum(k,i)*qsat(k,i)
            end if

            tt(k) = temperature_1d(k)*(1.+1.61*scalars(index_qv,k,i))

          end do
		
          do itrp = 1,25
            do k=1,nz1				
              rr(k,i)  = (pp(k,i)/(rgas*zz(k,i)) - rb(k,i)*(tt(k)-t0b))/tt(k)
            end do

            ppi(1) = p0-.5*dzw(1)*gravity                         &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

            ppi(1) = ppi(1)-ppb(1,i)
            do k=1,nz1-1

!              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                     &
!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)   &
!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))

               ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                                  &
                             ( (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i))*fzp(k+1)   &
                             + (rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))*fzm(k+1) )

            end do

            do k=1,nz1
              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
            end do

          end do  ! end inner iteration loop itrp

        end do  ! end outer iteration loop itr

        do k=1,nz1	
          p (k,i) = ((ppb(k,i)+pp(k,i))/p0)**(rgas/cp)
          t (k,i) = tt(k)/p(k,i)
          rt (k,i) = t(k,i)*rr(k,i)+rb(k,i)*(t(k,i)-tb(k,i))
          rho_zz (k,i) = rb(k,i) + rr(k,i)
        end do

        !calculation of surface pressure:
        surface_pressure(i) = 0.5*dzw(1)*gravity                                    &
                        * (1.25*(rr(1,i) + rb(1,i)) * (1. + scalars(index_qv,1,i))  &
                        -  0.25*(rr(2,i) + rb(2,i)) * (1. + scalars(index_qv,2,i)))
        surface_pressure(i) = surface_pressure(i) + pp(1,i) + ppb(1,i)

      end do  ! end loop over cells

      !write(0,*)
      !write(0,*) '--- initialization of water vapor:'
      !do iCell = 1, nCells
      !   if(iCell == 1 .or. iCell == nCells) then
      !      do k = nz1, 1, -1
      !         write(0,202) iCell,k,t(k,iCell),relhum(k,iCell),qsat(k,iCell),scalars(index_qv,k,iCell)
      !      end do
      !      write(0,*)
      !   end if
      !end do

      lat_pert = latitude_pert*pii/180.
      lon_pert = longitude_pert*pii/180.

      do iEdge=1,nEdges

         vtx1 = verticesOnEdge(1,iEdge)
         vtx2 = verticesOnEdge(2,iEdge)
         lat1 = latVertex(vtx1)
         lat2 = latVertex(vtx2)
         iCell1 = cellsOnEdge(1,iEdge)
         iCell2 = cellsOnEdge(2,iEdge)
         flux = (0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)

         if (config_init_case == 2) then
            r_pert = sphere_distance( latEdge(iEdge), lonEdge(iEdge), &
                                      lat_pert, lon_pert, 1.0_RKIND)/(pert_radius)
            u_pert = u_perturbation*exp(-r_pert**2)*(lat2-lat1) * sphere_radius / dvEdge(iEdge)

         else if (config_init_case == 3) then
            lon_Edge = lonEdge(iEdge)
            u_pert = u_perturbation*cos(k_x*(lon_Edge - lon_pert)) &
                         *(0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)
         else
            u_pert = 0.0
         end if

         if (rebalance) then

           call init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,dvEdge(iEdge),sphere_radius,u0,nz1,nlat)
           do k=1,nVertLevels
             fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))
             u(k,iEdge) = fluxk + u_pert
           end do

         else 

           do k=1,nVertLevels
             etavs = (0.5*(ppb(k,iCell1)+ppb(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
             fluxk = u0*flux*(cos(etavs)**1.5)
             u(k,iEdge) = fluxk + u_pert
           end do

         end if

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nz1
            ru(k,iEdge)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,iEdge)
         end do

      !
      ! Generate rotated Coriolis field
      !

         fEdge(iEdge) = 2.0 * omega_e * &
                                       ( -cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha_grid) + &
                                         sin(latEdge(iEdge)) * cos(alpha_grid) &
                                       )
      end do

      do iVtx=1,nVertices
         fVertex(iVtx) = 2.0 * omega_e * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha_grid) + &
                                          sin(latVertex(iVtx)) * cos(alpha_grid) &
                                         )
      end do

      !
      !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
      !

      !
      !     pre-calculation z-metric terms in omega eqn.
      !
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 1, nVertLevels

            if (config_theta_adv_order == 2) then

               z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.

            else if (config_theta_adv_order == 3 .or. config_theta_adv_order ==4) then !theta_adv_order == 3 or 4 

               d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
               d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)

!  WCS fix 20120711

                  do i=1, nEdgesOnCell(cell1)
                     if ( cellsOnCell(i,cell1) > 0)       &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
                  end do
                  do i=1, nEdgesOnCell(cell2)
                     if ( cellsOnCell(i,cell2) > 0)       &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
                  end do             

               z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
                             - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

               if (config_theta_adv_order == 3) then
                  z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
               else
                  z_edge3 = 0.
               end if

            end if

               zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1)
               zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2)
               zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1)
               zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2)

         end do

      end do

      ! for including terrain
      rw = 0.0
      w = 0.0
      do iEdge = 1,nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         do k = 2, nVertLevels
            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
            rw(k,cell2) = rw(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
            rw(k,cell1) = rw(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux

            if (config_theta_adv_order ==3) then 
               rw(k,cell2) = rw(k,cell2)    &
                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
               rw(k,cell1) = rw(k,cell1)    &
                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
            end if

         end do

      end do

      ! Compute w from rho_zz and rw
      do iCell=1,nCells
         do k=2,nVertLevels
            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
         end do
      end do


      !
      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
      !
      v(:,:) = 0.0
      do iEdge = 1, nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1, nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
           end do
         end do
      end do

      do i=1,10
        psurf = (cf1*(ppb(1,i)+pp(1,i)) + cf2*(ppb(2,i)+pp(2,i)) + cf3*(ppb(3,i)+pp(3,i)))/100.

            psurf = (ppb(1,i)+pp(1,i)) + .5*dzw(1)*gravity        &
                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))

        write(0,*) ' i, psurf, lat ',i,psurf,latCell(i)*180./3.1415828
      end do

      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine init_atm_case_jw

!----------------------------------------------------------------------------------------------------------

   subroutine init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1_in,lat2_in,dvEdge,a,u0,nz1,nlat)

      implicit none
   
      integer, intent(in) :: nz1,nlat
      real (kind=RKIND), dimension(nz1,nlat), intent(in) :: u_2d,etavs_2d
      real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
      real (kind=RKIND), dimension(nz1), intent(out) :: flux_zonal
      real (kind=RKIND), intent(in) :: lat1_in, lat2_in, dvEdge, a, u0
   
      integer :: k,i
      real (kind=RKIND) :: lat1, lat2, w1, w2
      real (kind=RKIND) :: dlat,da,db
   
      lat1 = abs(lat1_in)
      lat2 = abs(lat2_in)
      if(lat2 <= lat1) then
        lat1 = abs(lat2_in)
        lat2 = abs(lat1_in)
      end if
   
      do k=1,nz1
        flux_zonal(k) = 0.
      end do
   
      do i=1,nlat-1
        if( (lat1 <= lat_2d(i+1)) .and. (lat2 >= lat_2d(i)) ) then
   
        dlat = lat_2d(i+1)-lat_2d(i)
        da = (max(lat1,lat_2d(i))-lat_2d(i))/dlat
        db = (min(lat2,lat_2d(i+1))-lat_2d(i))/dlat
        w1 = (db-da) -0.5*(db-da)**2
        w2 = 0.5*(db-da)**2
   
        do k=1,nz1
          flux_zonal(k) = flux_zonal(k) + w1*u_2d(k,i) + w2*u_2d(k,i+1)
        end do
   
        end if
   
      end do
   
   !  renormalize for setting cell-face fluxes
   
      do k=1,nz1
        flux_zonal(k) = sign(1.0_RKIND,lat2_in-lat1_in)*flux_zonal(k)*dlat*a/dvEdge/u0
      end do
        
   end subroutine init_atm_calc_flux_zonal


   !SHP-balance
   subroutine init_atm_recompute_geostrophic_wind(u_2d,rho_2d,pp_2d,qv_2d,lat_2d,zz_2d,zx_2d,     &
                                         cf1,cf2,cf3,fzm,fzp,rdzw,nz1,nlat,dlat,rad)

      implicit none
   
      integer, intent(in) :: nz1,nlat
      real (kind=RKIND), dimension(nz1,nlat), intent(inout) :: u_2d
      real (kind=RKIND), dimension(nz1,nlat), intent(in) :: rho_2d, pp_2d, qv_2d, zz_2d
      real (kind=RKIND), dimension(nz1,nlat-1), intent(in) :: zx_2d
      real (kind=RKIND), dimension(nlat), intent(in) :: lat_2d
      real (kind=RKIND), dimension(nz1), intent(in) :: fzm, fzp, rdzw
      real (kind=RKIND), intent(in) :: cf1, cf2, cf3, dlat, rad
   
      !local variable
      real (kind=RKIND), dimension(nz1,nlat-1) :: pgrad, ru, u
      real (kind=RKIND), dimension(nlat-1) :: f
      real (kind=RKIND), dimension(nz1+1)  :: dpzx
   
   !   real (kind=RKIND), parameter :: omega_e = 7.29212e-05
      real (kind=RKIND) :: omega_e
   
      real (kind=RKIND) :: rdx, qtot, r_earth, phi
      integer :: k,i, itr
   
      r_earth  = rad
      omega_e = omega
      rdx = 1./(dlat*r_earth)
   
      do i=1,nlat-1
        do k=1,nz1
          pgrad(k,i) = rdx*(pp_2d(k,i+1)/zz_2d(k,i+1)-pp_2d(k,i)/zz_2d(k,i))
        end do
   
        dpzx(:) = 0.
   
        k=1
        dpzx(k) = .5*zx_2d(k,i)*(cf1*(pp_2d(k  ,i+1)+pp_2d(k  ,i))        &
                                +cf2*(pp_2d(k+1,i+1)+pp_2d(k+1,i))        &
                                +cf3*(pp_2d(k+2,i+1)+pp_2d(k+2,i)))
        do k=2,nz1
           dpzx(k) = .5*zx_2d(k,i)*(fzm(k)*(pp_2d(k  ,i+1)+pp_2d(k  ,i))   &
                                   +fzp(k)*(pp_2d(k-1,i+1)+pp_2d(k-1,i)))
        end do
   
        do k=1,nz1
           pgrad(k,i) = pgrad(k,i) - rdzw(k)*(dpzx(k+1)-dpzx(k))
        end do
      end do
   
   
      !initial value of v and rv -> that is from analytic sln. 
      do i=1,nlat-1
         do k=1,nz1
            u(k,i) = .5*(u_2d(k,i)+u_2d(k,i+1))
            ru(k,i) = u(k,i)*(rho_2d(k,i)+rho_2d(k,i+1))*.5
         end do
      end do
   
      write(0,*) "MAX U wind before REBALANCING ---->", maxval(abs(u))
   
      !re-calculate geostrophic wind using iteration 
      do itr=1,50
      do i=1,nlat-1
         phi = (lat_2d(i)+lat_2d(i+1))/2.
         f(i) = 2.*omega_e*sin(phi)
         do k=1,nz1
            if (f(i).eq.0.) then
              ru(k,i) = 0.
            else
              qtot = .5*(qv_2d(k,i)+qv_2d(k,i+1))
              ru(k,i) = - ( 1./(1.+qtot)*pgrad(k,i) + tan(phi)/r_earth*u(k,i)*ru(k,i) )/f(i)
            end if
              u(k,i) = ru(k,i)*2./(rho_2d(k,i)+rho_2d(k,i+1))
         end do
      end do
      end do
   
      write(0,*) "MAX U wind after REBALANCING ---->", maxval(abs(u))
   
      !update 2d ru
      do i=2,nlat-1
        do k=1,nz1
          u_2d(k,i) = (ru(k,i-1)+ru(k,i))*.5
        end do
      end do
   
      i=1
      do k=1,nz1
         u_2d(k,i) = (3.*u_2d(k,i+1)-u_2d(k,i+2))*.5
      end do
      i=nlat
      do k=1,nz1
         u_2d(k,i) = (3.*u_2d(k,i-1)-u_2d(k,i-2))*.5
      end do

   end subroutine init_atm_recompute_geostrophic_wind



   real (kind=RKIND) function env_qv( z, temperature, pressure, rh_max )

      implicit none
      real (kind=RKIND) :: z, temperature, pressure, ztr, es, qvs, p0, rh_max

      p0 = 100000.

!      ztr = 5000.
!
!      if(z .gt. ztr) then
!         env_qv = 0.
!      else
!         if(z.lt.2000.) then
!            env_qv = .5
!         else
!            env_qv = .5*(1.-(z-2000.)/(ztr-2000.))
!         end if
!      end if

       if (pressure .lt. 50000. ) then
           env_qv = 0.0
       else
           env_qv = (1.-((p0-pressure)/50000.)**1.25)
       end if

       env_qv = min(rh_max,env_qv)

! env_qv is the relative humidity, turn it into mixing ratio
       if (temperature .gt. 273.15) then
           es  = 1000.*0.6112*exp(17.67*(temperature-273.15)/(temperature-29.65))
       else
           es  = 1000.*0.6112*exp(21.8745584*(temperature-273.16)/(temperature-7.66))
       end if
       qvs = (287.04/461.6)*es/(pressure-es)

       ! qvs =  380.*exp(17.27*(temperature-273.)/(temperature-36.))/pressure

        env_qv = env_qv*qvs

   end function env_qv


   subroutine decouple_variables(mesh, nCells, nVertLevels, state, diag)

      implicit none

      type (mpas_pool_type), intent(in) :: mesh
      integer, intent(in) :: nCells
      integer, intent(in) :: nVertLevels
      type (mpas_pool_type), intent(inout) :: state
      type (mpas_pool_type), intent(inout) :: diag

      integer :: iCell, iEdge, k

      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: rdzw
      real (kind=RKIND), dimension(:,:), pointer :: zz, pp, ppb, rho, rho_zz, theta, theta_m
      real (kind=RKIND), dimension(:), pointer :: surface_pressure
      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      integer, pointer :: index_qv

      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'rdzw', rdzw)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(diag, 'pressure_p', pp)
      call mpas_pool_get_array(diag, 'pressure_base', ppb)
      call mpas_pool_get_array(diag, 'surface_pressure', surface_pressure)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
      call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)

      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
     
      ! Compute surface pressure
      do iCell=1,nCells
         surface_pressure(iCell) = 0.5*gravity/rdzw(1)                                        &
                                   * (1.25* rho_zz(1,iCell) * (1. + scalars(index_qv, 1, iCell))  &
                                      -  0.25* rho_zz(2,iCell) * (1. + scalars(index_qv, 2, iCell)))
         surface_pressure(iCell) = surface_pressure(iCell) + pp(1,iCell) + ppb(1,iCell)
      end do


      ! Compute rho and theta from rho_zz and theta_m
      do iCell=1,nCells
         do k=1,nVertLevels
            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
            theta(k,iCell) = theta_m(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
         end do
      end do

   end subroutine decouple_variables


   real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance


end module init_atm_cases
