module mpas_atm_boundaries

    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
    use mpas_kind_types, only : RKIND
    use mpas_geometry_utils, only : mpas_arc_length

    integer, parameter :: OUTSIDE = 0
    integer, parameter :: INSIDE = 1
    integer, parameter :: BOUNDARY = 2
    integer, parameter :: BOUNDARY2 = 3
    integer, parameter :: BOUNDARY3 = 4
    integer, parameter :: BOUNDARY4 = 5
    integer, parameter :: BOUNDARY5 = 6


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

end module mpas_atm_boundaries
