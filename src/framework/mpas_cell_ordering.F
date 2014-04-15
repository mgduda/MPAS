module mpas_cell_ordering

   use mpas_dmpar_types
   use mpas_kind_types
   use mpas_sort

   public :: MPAS_orderCells, MPAS_permuteCellFields

   private

   integer :: max_queue_length
   integer :: max_dictionary_size

   integer :: queue_head = 0
   integer :: queue_tail = 0
   integer :: queue_size = 0
   integer, dimension(:), allocatable :: queue_array

   integer :: int_size = 32
   integer, dimension(:), allocatable :: dictionary_array


   contains


   subroutine MPAS_orderCells(dminfo, nCells, maxEdges, indexToCellID, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, cellOrder)

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

      integer :: iCell, i, j, k, firstIdx, bdyCell, localCell
      integer, dimension(2,nCells) :: global_local
      logical :: found_non_cell


      !
      ! Begin by finding a cell that is on the boundary of the block; if no valid boundary cell
      !   is returned, we will assume that there are no boundaries and we can select any cell
      !
      bdyCell = find_bdy_cell(nCells, maxEdges, indexToCellID, nEdgesOnCell, cellsOnCell)
      if (bdyCell == -1) bdyCell = indexToCellID(1)

      max_dictionary_size = nCells
      allocate(dictionary_array(max_dictionary_size))
      call dictionary_reset()

      max_queue_length = nCells
      allocate(queue_array(0:max_queue_length-1))
      call queue_reset()


      !
      ! Create a table to help us translate a global index to a local index efficiently
      !
      do i=1,nCells
         global_local(1,i) = indexToCellID(i)
         global_local(2,i) = i
      end do
      call MPAS_quicksort(nCells, global_local)

      iCell = nCells
      
      k = MPAS_binary_search(global_local, 2, 1, nCells, bdyCell)
      localCell = global_local(2,k)
      cellOrder(iCell) = bdyCell     ! add global ID to cellOrder array
      iCell = iCell - 1
      call queue_insert(localCell)
      call dictionary_insert(localCell)

      !
      ! Process cells
      !
      do while (queue_size > 0)
         i = queue_remove()


         ! Find first index that we can add to the queue after finding an index that we cannot
         found_non_cell = .false.
         firstIdx = -1
         do j=1,nEdgesOnCell(i)
            k = MPAS_binary_search(global_local, 2, 1, nCells, cellsOnCell(j,i))
            if (k <= nCells) then
               if (.not. dictionary_search(localCell)) then
                  if (found_non_cell) then
                     firstIdx = j
                     exit
                  end if
               else
                  found_non_cell = .true.
               end if
            else
               found_non_cell = .true.
            end if
         end do
         if (firstIdx == -1) then
            firstIdx = 1
         end if
         

         do j=firstIdx,nEdgesOnCell(i)

            k = MPAS_binary_search(global_local, 2, 1, nCells, cellsOnCell(j,i))
            if (k <= nCells) then
               localCell = global_local(2,k)
               if (.not. dictionary_search(localCell)) then
                  cellOrder(iCell) = global_local(1,k)     ! add global ID to cellOrder array
                  iCell = iCell - 1
                  call queue_insert(localCell)
                  call dictionary_insert(localCell)
               else
               end if
            end if
         end do

         do j=1,firstIdx-1

            k = MPAS_binary_search(global_local, 2, 1, nCells, cellsOnCell(j,i))
            if (k <= nCells) then
               localCell = global_local(2,k)
               if (.not. dictionary_search(localCell)) then
                  cellOrder(iCell) = global_local(1,k)     ! add global ID to cellOrder array
                  iCell = iCell - 1
                  call queue_insert(localCell)
                  call dictionary_insert(localCell)
               else
               end if
            end if
         end do
    
      end do


      deallocate(dictionary_array)
      deallocate(queue_array)

   end subroutine MPAS_orderCells


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


      !
      ! Knowing the source for all cells in the new arrays, we can permute into a temp array,
      !   then copy back to the original arrays
      !
      allocate(temp1(nCells))

      ! indexToCellID
      do i=1,nCells
         temp1(i) = indexToCellID(src(i))
      end do
      indexToCellID(:) = temp1(:)

      ! nEdgesOnCell
      do i=1,nCells
         temp1(i) = nEdgesOnCell(src(i))
      end do
      nEdgesOnCell(:) = temp1(:)

      deallocate(temp1)

      allocate(temp2(maxEdges, nCells))

      ! cellsOnCell
      do i=1,nCells
         temp2(:,i) = cellsOnCell(:,src(i))
      end do
      cellsOnCell(:,:) = temp2(:,:)

      ! edgesOnCell
      do i=1,nCells
         temp2(:,i) = edgesOnCell(:,src(i))
      end do
      edgesOnCell(:,:) = temp2(:,:)

      ! verticesOnCell
      do i=1,nCells
         temp2(:,i) = verticesOnCell(:,src(i))
      end do
      verticesOnCell(:,:) = temp2(:,:)

      deallocate(temp2)

   end subroutine MPAS_permuteCellFields   


   integer function find_bdy_cell(nCells, maxEdges, indexToCellID, nEdgesOnCell, cellsOnCell)

      implicit none

      integer, intent(in) :: nCells
      integer, intent(in) :: maxEdges
      integer, dimension(nCells), intent(in) :: indexToCellID
      integer, dimension(nCells), intent(in) :: nEdgesOnCell
      integer, dimension(maxEdges, nCells), intent(in) :: cellsOnCell

      integer :: i, j
      integer, dimension(2,nCells) :: lookup


      !
      ! Add all cell IDs to a reverse look-up array
      !
      do i=1,nCells
         lookup(1,i) = indexToCellID(i)
         lookup(2,i) = i
      end do
      call MPAS_quicksort(nCells, lookup)

      !
      ! Loop through cells until we find one with a neighbor not in the lookup table
      !
      do i=1,nCells
         do j=1,nEdgesOnCell(i)
            find_bdy_cell = MPAS_binary_search(lookup, 2, 1, nCells, cellsOnCell(j,i))
            if (find_bdy_cell > nCells) then
               find_bdy_cell = indexToCellID(i)
               return
            end if 
         end do
      end do
    
      find_bdy_cell = -1

   end function find_bdy_cell


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Insert a new integer into the queue
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine queue_insert(i)
   
       implicit none
   
       integer, intent(in) :: i
   
       if (queue_size == max_queue_length) then
           write(0,*) 'Error: queue overrun'
           return
       end if
       queue_size = queue_size + 1
       queue_array(queue_head) = i
       queue_head = mod(queue_head + 1, max_queue_length)
   
   end subroutine queue_insert
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Remove the oldest integer from the queue
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function queue_remove()
   
       implicit none
   
       if (queue_size <= 0) then
           write(0,*) 'Error: queue underrun'
           queue_remove = -1
           return
       end if
       queue_size = queue_size - 1
       queue_remove = queue_array(queue_tail)
       queue_tail = mod(queue_tail + 1, max_queue_length)
   
   end function queue_remove
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Reset the queue to an empty state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine queue_reset()
   
       implicit none
   
       queue_head = 0
       queue_tail = 0
       queue_size = 0
   
   end subroutine queue_reset
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Insert an integer into the dictionary
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine dictionary_insert(i)
   
       implicit none
   
       integer, intent(in) :: i
   
       integer :: n_integer
       integer :: n_bit
   
       n_integer = ((i-1) / int_size) + 1
       n_bit = mod((i-1), int_size)
   
       if (n_integer > max_dictionary_size) then
           write(0,*) 'Error: dictionary insert out of bounds'
           return
       end if
   
       dictionary_array(n_integer) = ibset(dictionary_array(n_integer), n_bit)
   
   end subroutine dictionary_insert
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Search for an integer in the dictionary
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   logical function dictionary_search(i)
   
       implicit none
   
       integer, intent(in) :: i
   
       integer :: n_integer
       integer :: n_bit
   
       n_integer = ((i-1) / int_size) + 1
       n_bit = mod((i-1), int_size)
   
       if (n_integer > max_dictionary_size) then
           write(0,*) 'Error: dictionary search out of bounds'
           dictionary_search = .false.
           return
       end if
   
       dictionary_search = btest(dictionary_array(n_integer), n_bit)
   
   end function dictionary_search
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Reset the dictionary to an empty state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine dictionary_reset()
   
       implicit none
   
       dictionary_array(:) = 0
   
   end subroutine dictionary_reset

end module mpas_cell_ordering
