module mpas_atm_boundaries

    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
    use mpas_kind_types, only : RKIND
    use mpas_geometry_utils, only : mpas_arc_length

    integer, parameter :: nSpecLayers = 2
    integer, parameter :: nRelaxLayers = 5
    integer, parameter :: OUTSIDE = -1
    integer, parameter :: INSIDE = 0
    integer, dimension(nSpecLayers + nRelaxLayers), parameter :: BOUNDARY = (/ 1, 2, 3, 4, 5, 6, 7 /)


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
        real (kind=RKIND) :: xc, yc, zc
        integer :: n, bmin, bmax


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

        call convert_lx(xc, yc, zc, 6371229.0, 40.0*3.1459265359/180.0, -95.0*3.14159265359/180.0)

        !
        ! Identify interior cells
        !
        do iCell=1,nCells
            if (mpas_arc_length(xc, yc, zc, &
                                xCell(iCell), yCell(iCell), zCell(iCell)) < 2500000.0_RKIND) then
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
                        bdyMaskCell(iCell) = BOUNDARY(1)
                        exit
                    end if
                end do
            end if
        end do


        !
        ! And so on...
        !
        do n = 2, nSpecLayers+nRelaxLayers
            do iCell=1,nCells
                if (bdyMaskCell(iCell) == OUTSIDE) then
                    do j=1,nEdgesOnCell(iCell)
                        if (bdyMaskCell(cellsOnCell(j,iCell)) == BOUNDARY(n-1)) then
                            bdyMaskCell(iCell) = BOUNDARY(n)
                            exit
                        end if
                    end do
                end if
            end do
        end do


        !
        ! With cell masks set, compute edge masks. An edge receives the largest
        ! mask value of either of its neighboring cells.
        !
        do iEdge=1,nEdges
            bmin = min(bdyMaskCell(cellsOnEdge(1,iEdge)), &
                       bdyMaskCell(cellsOnEdge(2,iEdge))) 

            bmax = max(bdyMaskCell(cellsOnEdge(1,iEdge)), &
                       bdyMaskCell(cellsOnEdge(2,iEdge))) 

            if (bmin /= bmax .and. bmin == OUTSIDE) then
                bdyMaskEdge(iEdge) = bmax
            else
                bdyMaskEdge(iEdge) = bmin
            end if
        end do


        !
        ! Compute vertex masks similar to edge masks: A vertex receives the largest
        ! mask value of any of its neighboring cells.
        !
        do iVertex=1,nVertices
            bmin = min(bdyMaskCell(cellsOnVertex(1,iVertex)), &
                       bdyMaskCell(cellsOnVertex(2,iVertex)), &
                       bdyMaskCell(cellsOnVertex(3,iVertex)))

            bmax = max(bdyMaskCell(cellsOnVertex(1,iVertex)), &
                       bdyMaskCell(cellsOnVertex(2,iVertex)), &
                       bdyMaskCell(cellsOnVertex(3,iVertex)))

            if (bmin /= bmax .and. bmin == OUTSIDE) then
                bdyMaskVertex(iVertex) = bmax
            else
                bdyMaskVertex(iVertex) = bmin
            end if
        end do

    end subroutine mpas_atm_setup_bdy_masks


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SUBROUTINE CONVERT_LX
    !
    ! Convert (lat,lon) to an (x, y, z) location on a sphere with specified radius.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine convert_lx(x, y, z, radius, lat, lon)

       implicit none

       real (kind=RKIND), intent(in) :: radius
       real (kind=RKIND), intent(out) :: x, y, z
       real (kind=RKIND), intent(in) :: lat, lon

       z = radius * sin(lat)
       x = radius * cos(lon) * cos(lat)
       y = radius * sin(lon) * cos(lat)

    end subroutine convert_lx

end module mpas_atm_boundaries
