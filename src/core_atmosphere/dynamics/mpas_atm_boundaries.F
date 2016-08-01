module mpas_atm_boundaries

    use mpas_derived_types, only : mpas_pool_type, mpas_clock_type, block_type, mpas_time_type, MPAS_NOW
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_subpool
    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_geometry_utils, only : mpas_arc_length
    use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_time

    integer, parameter :: OUTSIDE = 0
    integer, parameter :: INSIDE = 1
    integer, parameter :: BOUNDARY = 2
    integer, parameter :: BOUNDARY2 = 3
    integer, parameter :: BOUNDARY3 = 4
    integer, parameter :: BOUNDARY4 = 5
    integer, parameter :: BOUNDARY5 = 6

    character(len=StrKIND), private :: LBC_intv_start


    contains


    subroutine mpas_atm_setup_bdy_masks(mesh)

        implicit none

        type (mpas_pool_type), intent(inout) :: mesh

        integer, pointer :: nCells, nEdges, nVertices
        integer :: iCell, iEdge, iVertex, j
        integer, dimension(:), pointer :: bdyMaskCell, bdyMaskEdge, bdyMaskVertex
        integer, dimension(:), pointer :: nEdgesOnCell
        integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, cellsOnVertex
        real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
        real (kind=RKIND) :: r


        call mpas_pool_get_array(mesh, 'bdyMaskCell', bdyMaskCell)
        call mpas_pool_get_array(mesh, 'bdyMaskEdge', bdyMaskEdge)
        call mpas_pool_get_array(mesh, 'bdyMaskVertex', bdyMaskVertex)
        call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
        call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
        call mpas_pool_get_array(mesh, 'xCell', xCell)
        call mpas_pool_get_array(mesh, 'yCell', yCell)
        call mpas_pool_get_array(mesh, 'zCell', zCell)

        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
        call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)

        !
        ! Zero-out boundary masks
        !
        bdyMaskCell(:) = OUTSIDE


        !
        ! Identify interior cells
        !
        do iCell=1,nCells
            if (mpas_arc_length(6371229.0_RKIND, 0.0_RKIND, 0.0_RKIND, &
                                xCell(iCell), yCell(iCell), zCell(iCell)) < 4000000.0_RKIND) then
                bdyMaskCell(iCell) = INSIDE
            end if
        end do


        !
        ! Build the first layer of boundary cells
        !
        do iCell=1,nCells
            if (bdyMaskCell(iCell) == OUTSIDE) then
                do j=1,nEdgesOnCell(iCell)
                    if (bdyMaskCell(cellsOnCell(j,iCell)) == INSIDE) then
                        bdyMaskCell(iCell) = BOUNDARY
                        exit
                    end if
                end do
            end if
        end do


        !
        ! And so on...
        !
        do iCell=1,nCells
            if (bdyMaskCell(iCell) == OUTSIDE) then
                do j=1,nEdgesOnCell(iCell)
                    if (bdyMaskCell(cellsOnCell(j,iCell)) == BOUNDARY) then
                        bdyMaskCell(iCell) = BOUNDARY2
                        exit
                    end if
                end do
            end if
        end do

        do iCell=1,nCells
            if (bdyMaskCell(iCell) == OUTSIDE) then
                do j=1,nEdgesOnCell(iCell)
                    if (bdyMaskCell(cellsOnCell(j,iCell)) == BOUNDARY2) then
                        bdyMaskCell(iCell) = BOUNDARY3
                        exit
                    end if
                end do
            end if
        end do

        do iCell=1,nCells
            if (bdyMaskCell(iCell) == OUTSIDE) then
                do j=1,nEdgesOnCell(iCell)
                    if (bdyMaskCell(cellsOnCell(j,iCell)) == BOUNDARY3) then
                        bdyMaskCell(iCell) = BOUNDARY4
                        exit
                    end if
                end do
            end if
        end do

        do iCell=1,nCells
            if (bdyMaskCell(iCell) == OUTSIDE) then
                do j=1,nEdgesOnCell(iCell)
                    if (bdyMaskCell(cellsOnCell(j,iCell)) == BOUNDARY4) then
                        bdyMaskCell(iCell) = BOUNDARY5
                        exit
                    end if
                end do
            end if
        end do


        !
        ! With cell masks set, compute edge masks. An edge receives the largest
        ! mask value of either of its neighboring cells.
        !
        do iEdge=1,nEdges
            bdyMaskEdge(iEdge) = max(bdyMaskCell(cellsOnEdge(1,iEdge)), &
                                     bdyMaskCell(cellsOnEdge(2,iEdge))) 
        end do


        !
        ! Compute vertex masks similar to edge masks: A vertex receives the largest
        ! mask value of any of its neighboring cells.
        !
        do iVertex=1,nVertices
            bdyMaskVertex(iVertex) = max(bdyMaskCell(cellsOnVertex(1,iVertex)), &
                                         bdyMaskCell(cellsOnVertex(2,iVertex)), &
                                         bdyMaskCell(cellsOnVertex(3,iVertex)))
        end do

    end subroutine mpas_atm_setup_bdy_masks


    subroutine mpas_atm_setup_bdy_tend(clock, block)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (block_type), intent(inout) :: block

        type (mpas_pool_type), pointer :: mesh
        type (mpas_pool_type), pointer :: state
        type (mpas_pool_type), pointer :: lbc

        real (kind=RKIND), dimension(:,:), pointer :: u
        real (kind=RKIND), dimension(:,:), pointer :: w
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

        integer :: ierr
        type (MPAS_Time_Type) :: currTime


write(0,*) 'Setting up next tendency'
        call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
        call mpas_get_time(curr_time=currTime, dateTimeString=LBC_intv_start, ierr=ierr)
write(0,*) '    will be valid starting '//trim(LBC_intv_start)

        call mpas_pool_get_array(state, 'u', u, 1)
        call mpas_pool_get_array(state, 'w', w, 1)
        call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
        call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
        call mpas_pool_get_array(state, 'scalars', scalars, 1)

        call mpas_pool_get_array(lbc, 'lbc_tend_u', lbc_tend_u)
        call mpas_pool_get_array(lbc, 'lbc_tend_w', lbc_tend_w)
        call mpas_pool_get_array(lbc, 'lbc_tend_theta_m', lbc_tend_theta_m)
        call mpas_pool_get_array(lbc, 'lbc_tend_rho_zz', lbc_tend_rho_zz)
        call mpas_pool_get_array(lbc, 'lbc_tend_scalars', lbc_tend_scalars)

        lbc_tend_u(:,:) = u(:,:)
        lbc_tend_w(:,:) = w(:,:)
        lbc_tend_theta_m(:,:) = theta_m(:,:)
        lbc_tend_rho_zz(:,:) = rho_zz(:,:)
        lbc_tend_scalars(:,:,:) = scalars(:,:,:)

    end subroutine mpas_atm_setup_bdy_tend


    subroutine mpas_atm_compute_bdy_tend(clock, block, lbc_start_time)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (block_type), intent(inout) :: block
        character (len=StrKIND), intent(out) :: lbc_start_time

        type (mpas_pool_type), pointer :: mesh
        type (mpas_pool_type), pointer :: state
        type (mpas_pool_type), pointer :: lbc
        real (kind=RKIND) :: dt

        real (kind=RKIND), dimension(:,:), pointer :: u
        real (kind=RKIND), dimension(:,:), pointer :: w
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

        integer :: ierr
        type (MPAS_Time_Type) :: currTime


        lbc_start_time = LBC_intv_start

dt = 10800.0

write(0,*) 'Computing tendency'
write(0,*) '    from '//trim(lbc_start_time)
        call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
        call mpas_get_time(curr_time=currTime, dateTimeString=LBC_intv_start, ierr=ierr)

        call mpas_pool_get_array(state, 'u', u, 1)
        call mpas_pool_get_array(state, 'w', w, 1)
        call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
        call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)
        call mpas_pool_get_array(state, 'scalars', scalars, 1)

        call mpas_pool_get_array(lbc, 'lbc_tend_u', lbc_tend_u)
        call mpas_pool_get_array(lbc, 'lbc_tend_w', lbc_tend_w)
        call mpas_pool_get_array(lbc, 'lbc_tend_theta_m', lbc_tend_theta_m)
        call mpas_pool_get_array(lbc, 'lbc_tend_rho_zz', lbc_tend_rho_zz)
        call mpas_pool_get_array(lbc, 'lbc_tend_scalars', lbc_tend_scalars)

        lbc_tend_u(:,:) = (u(:,:) - lbc_tend_u(:,:)) / dt
        lbc_tend_w(:,:) = (w(:,:) - lbc_tend_w(:,:)) / dt
        lbc_tend_theta_m(:,:) = (theta_m(:,:) - lbc_tend_theta_m(:,:)) / dt
        lbc_tend_rho_zz(:,:) = (rho_zz(:,:) - lbc_tend_rho_zz(:,:)) / dt
        lbc_tend_scalars(:,:,:) = (scalars(:,:,:) - lbc_tend_scalars(:,:,:)) / dt

    end subroutine mpas_atm_compute_bdy_tend

end module mpas_atm_boundaries
