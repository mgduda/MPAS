module mpas_cell_ordering

   use mpas_dmpar_types
   use mpas_kind_types
   use mpas_sort


   contains


   subroutine mpas_orderCells(dminfo, nCells, maxEdges, indexToCellID, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, cellOrder)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nCells
      integer, intent(in) :: maxEdges
      integer, dimension(nCells), intent(in) :: indexToCellID
      real (kind=RKIND), dimension(nCells), intent(in) :: xCell
      real (kind=RKIND), dimension(nCells), intent(in) :: yCell
      real (kind=RKIND), dimension(nCells), intent(in) :: zCell
      integer, dimension(maxEdges, nCells), intent(in) :: cellsOnCell
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(nCells), intent(out) :: cellOrder

      integer :: i
      real (kind=RKIND), dimension(2,nCells) :: temp

      do i=1,nCells
         temp(1,i) = xCell(i)
         temp(2,i) = indexToCellID(i)
      end do

      call MPAS_quicksort(nCells, temp)

      do i=1,nCells
         cellOrder(i) = temp(2,i)
      end do

   end subroutine mpas_orderCells


   subroutine MPAS_permuteCellFields(nCells, maxEdges, indexToCellID, nEdgesOnCell, cellsOnCell, verticesOnCell, edgesOnCell, cellOrder)

      implicit none

      integer, intent(in) :: nCells
      integer, intent(in) :: maxEdges
      integer, dimension(nCells), intent(inout) :: indexToCellID
      integer, dimension(nCells), intent(inout) :: nEdgesOnCell
      integer, dimension(maxEdges, nCells), intent(inout) :: cellsOnCell
      integer, dimension(maxEdges, nCells), intent(inout) :: verticesOnCell
      integer, dimension(maxEdges, nCells), intent(inout) :: edgesOnCell
      integer, dimension(nCells), intent(inout) :: cellOrder

      integer :: i, j
      integer, dimension(nCells) :: src
      integer, dimension(:), allocatable :: temp1
      integer, dimension(:,:), allocatable :: temp2

      ! 
      ! Given the current list of cell global indices (indexToCellID) and the new list (cellOrder), 
      !   we first create a permutation
      ! 
      allocate(temp2(2,nCells))
      do i=1,nCells
         temp2(1,i) = indexToCellID(i)
         temp2(2,i) = i
      end do
      call MPAS_quicksort(nCells, temp2)

      do i=1,nCells
         j = MPAS_binary_search(temp2, 2, 1, nCells, cellOrder(i))
         if (j <= nCells) then
            src(i) = temp2(2,j) 
         else
            write(0,*) 'Strange errors in MPAS_permuteCellFields...'
         end if
      end do
      deallocate(temp2)

do i=1,10
   write(0,*) 'src = ', src(i)
end do

      !
      ! Knowing the source for all cells in the new arrays, we can permute into a temp array,
      !   then copy back to the original arrays
      !
      allocate(temp1(nCells))
      do i=1,nCells
         temp1(i) = indexToCellID(src(i))
      end do
      indexToCellID(:) = temp1(:)

      do i=1,nCells
         temp1(i) = nEdgesOnCell(src(i))
      end do
      nEdgesOnCell(:) = temp1(:)
      deallocate(temp1)

      allocate(temp2(maxEdges, nCells))
      do i=1,nCells
         temp2(:,i) = cellsOnCell(:,src(i))
      end do
      cellsOnCell(:,:) = temp2(:,:)

      do i=1,nCells
         temp2(:,i) = edgesOnCell(:,src(i))
      end do
      edgesOnCell(:,:) = temp2(:,:)

      do i=1,nCells
         temp2(:,i) = verticesOnCell(:,src(i))
      end do
      verticesOnCell(:,:) = temp2(:,:)
      deallocate(temp2)

   end subroutine MPAS_permuteCellFields   

end module mpas_cell_ordering
