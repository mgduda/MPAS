module mpas_debugging

   use mpas_grid_types
   use mpas_dmpar_types

   logical :: write_fields
   type (mpas_pool_type), pointer :: all_fields 
   type (dm_info), pointer :: dminfo 

   public :: mpas_debugging_init, mpas_debugging_check, mpas_debugging_finalize

   private

#ifdef _MPI
include 'mpif.h'
   integer, parameter :: MPI_INTEGERKIND = MPI_INTEGER

#ifdef SINGLE_PRECISION
   integer, parameter :: MPI_REALKIND = MPI_REAL
#else
   integer, parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
#endif
#endif

   type exchange_list
      integer :: procID
      integer :: nlist
      integer, dimension(:), pointer :: list
      type (exchange_list), pointer :: next
      real (kind=RKIND), dimension(:), pointer :: rbuffer
      integer, dimension(:), pointer           :: ibuffer
      integer :: reqID
   end type exchange_list


   interface dmpar_alltoall_field
      module procedure dmpar_alltoall_field1dInteger
      module procedure dmpar_alltoall_field2dInteger
      module procedure dmpar_alltoall_field1dReal
      module procedure dmpar_alltoall_field2dReal
      module procedure dmpar_alltoall_field3dReal
   end interface

   integer, pointer :: nCells, nCellsSolve, &
                       nEdges, nEdgesSolve, &
                       nVertices, nVerticesSolve 
   integer :: nCellsGlobal, nEdgesGlobal, nVerticesGlobal

   type (exchange_list), pointer :: cellSendList, cellRecvList
   type (exchange_list), pointer :: edgeSendList, edgeRecvList
   type (exchange_list), pointer :: vertexSendList, vertexRecvList


   contains


   subroutine mpas_debugging_init(arg_dminfo, arg_write_fields, field_pool, mesh_pool)

      implicit none

      type (dm_info), pointer :: arg_dminfo 
      logical, intent(in) :: arg_write_fields
      type (mpas_pool_type), pointer :: field_pool
      type (mpas_pool_type), pointer :: mesh_pool

      integer :: i
      integer :: mpi_ierr
      integer, dimension(:), pointer :: indexToCellID, indexToEdgeID, indexToVertexID
      integer, dimension(:), pointer :: neededList

      dminfo => arg_dminfo
      write_fields = arg_write_fields
      all_fields => field_pool

      call mpas_pool_get_dimension(mesh_pool, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh_pool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh_pool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh_pool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(mesh_pool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh_pool, 'nVerticesSolve', nVerticesSolve)

      call MPI_Allreduce(nCellsSolve, nCellsGlobal, 1, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
      call MPI_Allreduce(nEdgesSolve, nEdgesGlobal, 1, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)
      call MPI_Allreduce(nVerticesSolve, nVerticesGlobal, 1, MPI_INTEGERKIND, MPI_SUM, dminfo % comm, mpi_ierr)

      call mpas_pool_get_array(mesh_pool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(mesh_pool, 'indexToEdgeID', indexToEdgeID)
      call mpas_pool_get_array(mesh_pool, 'indexToVertexID', indexToVertexID)

      if (dminfo % my_proc_id == 0) then
         allocate(neededList(nCellsGlobal))
         do i=1,nCellsGlobal
            neededList(i) = i
         end do
         call dmpar_get_owner_list(dminfo, &
                                   nCellsSolve, nCellsGlobal, &
                                   indexToCellID(1:nCellsSolve), neededList, &
                                   cellSendList, cellRecvList)
         deallocate(neededList)

         allocate(neededList(nEdgesGlobal))
         do i=1,nEdgesGlobal
            neededList(i) = i
         end do
         call dmpar_get_owner_list(dminfo, &
                                   nEdgesSolve, nEdgesGlobal, &
                                   indexToEdgeID(1:nEdgesSolve), neededList, &
                                   edgeSendList, edgeRecvList)
         deallocate(neededList)

         allocate(neededList(nVerticesGlobal))
         do i=1,nVerticesGlobal
            neededList(i) = i
         end do
         call dmpar_get_owner_list(dminfo, &
                                   nVerticesSolve, nVerticesGlobal, &
                                   indexToVertexID(1:nVerticesSolve), neededList, &
                                   vertexSendList, vertexRecvList)
         deallocate(neededList)
      else
         allocate(neededList(0))
         call dmpar_get_owner_list(dminfo, &
                                   nCellsSolve, 0, &
                                   indexToCellID(1:nCellsSolve), neededList, &
                                   cellSendList, cellRecvList)
         call dmpar_get_owner_list(dminfo, &
                                   nEdgesSolve, 0, &
                                   indexToEdgeID(1:nEdgesSolve), neededList, &
                                   edgeSendList, edgeRecvList)
         call dmpar_get_owner_list(dminfo, &
                                   nVerticesSolve, 0, &
                                   indexToVertexID(1:nVerticesSolve), neededList, &
                                   vertexSendList, vertexRecvList)
         deallocate(neededList)
      end if

      if (write_fields) then
         open(91,file='checkpoints.dat',form='unformatted',status='replace',action='write')
      else
         open(91,file='checkpoints.dat',form='unformatted',status='old',action='read')
      end if

   end subroutine mpas_debugging_init


   subroutine mpas_debugging_check(checkpoint, fieldname)

        implicit none

        character(len=*), intent(in) :: checkpoint

        integer :: i, j, k, outdim, outdimglobal
        type (mpas_pool_iterator_type) :: poolItr
        type (field1DInteger), pointer :: int1DField
        type (field2DInteger), pointer :: int2DField
        type (field1DReal), pointer :: real1DField
        type (field2DReal), pointer :: real2DField
        type (field3DReal), pointer :: real3DField
        character(len=*), intent(in), optional :: fieldname

        integer, dimension(:), pointer   :: arrayOut_i1, arrayComp_i1
        integer, dimension(:,:), pointer :: arrayOut_i2, arrayComp_i2
        real (kind=RKIND), dimension(:), pointer     :: arrayOut_r1, arrayComp_r1
        real (kind=RKIND), dimension(:,:), pointer   :: arrayOut_r2, arrayComp_r2
        real (kind=RKIND), dimension(:,:,:), pointer :: arrayOut_r3, arrayComp_r3
        type (exchange_list), pointer :: sendList, recvList
        type (mpas_pool_field_info_type) :: field_info


        if (.not. present(fieldname)) then

           call mpas_pool_begin_iteration(all_fields)
           do while (mpas_pool_get_next_member(all_fields, poolItr))
               if (poolItr % memberType == MPAS_POOL_FIELD) then
                   if (poolItr % dataType == MPAS_POOL_REAL) then
   
                       if (poolItr % nDims == 1) then
                           call mpas_pool_get_field(all_fields, poolItr % memberName, real1DField)
                           call compare_r1(real1DField, checkpoint)
   
                       else if (poolItr % nDims == 2) then
   
                           call mpas_pool_get_field(all_fields, poolItr % memberName, real2DField)
                           call compare_r2(real2DField, checkpoint)
   
                       else if (poolItr % nDims == 3) then
   
                           call mpas_pool_get_field(all_fields, poolItr % memberName, real3DField)
                           call compare_r3(real3DField, checkpoint)
   
                       end if
   
                   else if (poolItr % dataType == MPAS_POOL_INTEGER) then
   
                       if (poolItr % nDims == 1) then
   
                           call mpas_pool_get_field(all_fields, poolItr % memberName, int1DField)
                           call compare_i1(int1DField, checkpoint)
   
                       else if (poolItr % nDims == 2) then
   
                           call mpas_pool_get_field(all_fields, poolItr % memberName, int2DField)
                           call compare_i2(int2DField, checkpoint)
   
                        end if
   
                   end if
               end if
           end do

        else

           call mpas_pool_get_field_info(all_fields, fieldname, field_info)
   
           select case(field_info % fieldType)
   
              case (MPAS_POOL_REAL)
   
                 if (field_info % nDims == 1) then
   
                    call mpas_pool_get_field(all_fields, fieldname, real1DField)
                    call compare_r1(real1DField, checkpoint)
   
                 else if (field_info % nDims == 2) then
   
                    call mpas_pool_get_field(all_fields, fieldname, real2DField)
                    call compare_r2(real2DField, checkpoint)
   
                 else if (field_info % nDims == 3) then
   
                    call mpas_pool_get_field(all_fields, fieldname, real3DField)
                    call compare_r3(real3DField, checkpoint)

                 end if
   
              case (MPAS_POOL_INTEGER)
   
                 if (field_info % nDims == 1) then
   
                    call mpas_pool_get_field(all_fields, fieldname, int1DField)
                    call compare_i1(int1DField, checkpoint)
   
                 else if (field_info % nDims == 2) then
   
                    call mpas_pool_get_field(all_fields, fieldname, int2DField)
                    call compare_i2(int2DField, checkpoint)
   
                 end if
   
           end select

        end if

   end subroutine mpas_debugging_check


   subroutine mpas_debugging_finalize()

      implicit none

      close(91)

   end subroutine mpas_debugging_finalize


   subroutine dmpar_get_owner_list(dminfo, &
                                   nOwnedList, nNeededList, &
                                   ownedList, neededList, &
                                   sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nOwnedList, nNeededList
      integer, dimension(nOwnedList), intent(in) :: ownedList
      integer, dimension(nNeededList), intent(in) :: neededList
      type (exchange_list), pointer :: sendList
      type (exchange_list), pointer :: recvList

      integer :: i, j, k, kk
      integer :: totalSize, nMesgRecv, nMesgSend, recvNeighbor, sendNeighbor, currentProc
      integer :: numToSend, numToRecv
      integer, dimension(nOwnedList) :: recipientList
      integer, dimension(2,nOwnedList) :: ownedListSorted
      integer, allocatable, dimension(:) :: ownerListIn, ownerListOut
      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: mpi_ierr, mpi_rreq, mpi_sreq

#ifdef _MPI
      allocate(sendList)
      allocate(recvList)
      nullify(sendList % next)
      nullify(recvList % next)
      sendListPtr => sendList
      recvListPtr => recvList

      do i=1,nOwnedList
         ownedListSorted(1,i) = ownedList(i)
         ownedListSorted(2,i) = i
      end do
      call quicksort(nOwnedList, ownedListSorted)

      call MPI_Allreduce(nNeededList, totalSize, 1, MPI_INTEGER, MPI_MAX, dminfo % comm, mpi_ierr)

      allocate(ownerListIn(totalSize))
      allocate(ownerListOut(totalSize))

      nMesgRecv = nNeededList
      ownerListIn(1:nNeededList) = neededList(1:nNeededList)

      recvNeighbor = mod(dminfo % my_proc_id + dminfo % nprocs - 1, dminfo % nprocs)
      sendNeighbor = mod(dminfo % my_proc_id + 1, dminfo % nprocs)

      do i=1, dminfo % nprocs

         recipientList(:) = -1
         numToSend = 0

         currentProc = mod(dminfo % my_proc_id + dminfo % nprocs - i + 1, dminfo % nprocs)
         do j=1,nMesgRecv
            if (ownerListIn(j) > 0) then
               k = binary_search(ownedListSorted, 2, 1, nOwnedList, ownerListIn(j))
               if (k <= nOwnedList) then
                  ownerListOut(j) = -1 * dminfo % my_proc_id
                  numToSend = numToSend + 1
                  recipientList(ownedListSorted(2,k)) = numToSend
               else
                  ownerListOut(j) = ownerListIn(j)
               end if
            else
               ownerListOut(j) = ownerListIn(j)
            end if
         end do

         if (numToSend > 0) then
            allocate(sendListPtr % next)
            sendListPtr => sendListPtr % next
            sendListPtr % procID = currentProc
            sendListPtr % nlist = numToSend
            allocate(sendListPtr % list(numToSend))
            nullify(sendListPtr % next)
            kk = 1
            do j=1,nOwnedList
               if (recipientList(j) /= -1) then
                  sendListPtr % list(recipientList(j)) = j
                  kk = kk + 1
               end if
            end do
         end if

         nMesgSend = nMesgRecv
         call MPI_Irecv(nMesgRecv, 1, MPI_INTEGER, recvNeighbor, i, dminfo % comm, mpi_rreq, mpi_ierr)
         call MPI_Isend(nMesgSend, 1, MPI_INTEGER, sendNeighbor, i, dminfo % comm, mpi_sreq, mpi_ierr)
         call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
         call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
         call MPI_Irecv(ownerListIn, nMesgRecv, MPI_INTEGER, recvNeighbor, i, dminfo % comm, mpi_rreq, mpi_ierr)
         call MPI_Isend(ownerListOut, nMesgSend, MPI_INTEGER, sendNeighbor, i, dminfo % comm, mpi_sreq, mpi_ierr)
         call MPI_Wait(mpi_rreq, MPI_STATUS_IGNORE, mpi_ierr)
         call MPI_Wait(mpi_sreq, MPI_STATUS_IGNORE, mpi_ierr)
      end do

      do i=0, dminfo % nprocs - 1

         numToRecv = 0
         do j=1,nNeededList
            if (ownerListIn(j) == -i) numToRecv = numToRecv + 1
         end do
         if (numToRecv > 0) then
            allocate(recvListPtr % next)
            recvListPtr => recvListPtr % next
            recvListPtr % procID = i
            recvListPtr % nlist = numToRecv
            allocate(recvListPtr % list(numToRecv))
            nullify(recvListPtr % next)
            kk = 1
            do j=1,nNeededList
               if (ownerListIn(j) == -i) then
                  recvListPtr % list(kk) = j
                  kk = kk + 1
               end if
            end do
         end if

      end do

      deallocate(ownerListIn)
      deallocate(ownerListOut)

      sendListPtr => sendList
      sendList => sendList % next
      deallocate(sendListPtr)

      recvListPtr => recvList
      recvList => recvList % next
      deallocate(recvListPtr)

#else
      allocate(recvList)
      recvList % procID = dminfo % my_proc_id
      recvList % nlist = nNeededList
      allocate(recvList % list(nNeededList))
      nullify(recvList % next)
      do j=1,nNeededList
         recvList % list(j) = j
      end do

      allocate(sendList)
      sendList % procID = dminfo % my_proc_id
      sendList % nlist = nOwnedList
      allocate(sendList % list(nOwnedList))
      nullify(sendList % next)
      do j=1,nOwnedList
         sendList % list(j) = j
      end do
#endif

   end subroutine dmpar_get_owner_list


   subroutine dmpar_alltoall_field1dInteger(dminfo, arrayIn, arrayOut, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, dimension(nOwnedList), intent(in) :: arrayIn
      integer, dimension(nNeededList), intent(inout) :: arrayOut
      integer, intent(in) :: nOwnedList, nNeededList
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(recvListPtr % list(i)) = arrayIn(sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            allocate(recvListPtr % ibuffer(recvListPtr % nlist))
            call MPI_Irecv(recvListPtr % ibuffer, recvListPtr % nlist, MPI_INTEGER, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            allocate(sendListPtr % ibuffer(sendListPtr % nlist))
            call packSendBuf1dInteger(nOwnedList, arrayIn, sendListPtr, 1, sendListPtr % nlist, &
                                   sendListPtr % ibuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % ibuffer, sendListPtr % nlist, MPI_INTEGER, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            call unpackRecvBuf1dInteger(nNeededList, arrayOut, recvListPtr, 1, recvListPtr % nlist, &
                                     recvListPtr % ibuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % ibuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % ibuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, arrayIn and arrayOut dims must match.'
         call dmpar_abort(dminfo)
      else
         arrayOut(:) = arrayIn(:)
      end if
#endif

   end subroutine dmpar_alltoall_field1dInteger


   subroutine dmpar_alltoall_field2dInteger(dminfo, arrayIn, arrayOut, dim1, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, nOwnedList, nNeededList
      integer, dimension(dim1,nOwnedList), intent(in) :: arrayIn
      integer, dimension(dim1,nNeededList), intent(inout) :: arrayOut
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i, d2

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(:,recvListPtr % list(i)) = arrayIn(:,sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * recvListPtr % nlist
            allocate(recvListPtr % ibuffer(d2))
            call MPI_Irecv(recvListPtr % ibuffer, d2, MPI_INTEGER, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * sendListPtr % nlist
            allocate(sendListPtr % ibuffer(d2))
            call packSendBuf2dInteger(1, dim1, nOwnedList, arrayIn, sendListPtr, 1, d2, &
                                   sendListPtr % ibuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % ibuffer, d2, MPI_INTEGER, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d2 = dim1 * recvListPtr % nlist
            call unpackRecvBuf2dInteger(1, dim1, nNeededList, arrayOut, recvListPtr, 1, d2, &
                                     recvListPtr % ibuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % ibuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % ibuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, arrayIn and arrayOut dims must match.'
         call dmpar_abort(dminfo)
      else
         arrayOut(:,:) = arrayIn(:,:)
      end if
#endif

   end subroutine dmpar_alltoall_field2dInteger


   subroutine dmpar_alltoall_field1dReal(dminfo, arrayIn, arrayOut, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      real (kind=RKIND), dimension(nOwnedList), intent(in) :: arrayIn
      real (kind=RKIND), dimension(nNeededList), intent(inout) :: arrayOut
      integer, intent(in) :: nOwnedList, nNeededList
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(recvListPtr % list(i)) = arrayIn(sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            allocate(recvListPtr % rbuffer(recvListPtr % nlist))
            call MPI_Irecv(recvListPtr % rbuffer, recvListPtr % nlist, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            allocate(sendListPtr % rbuffer(sendListPtr % nlist))
            call packSendBuf1dReal(nOwnedList, arrayIn, sendListPtr, 1, sendListPtr % nlist, &
                                   sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, sendListPtr % nlist, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            call unpackRecvBuf1dReal(nNeededList, arrayOut, recvListPtr, 1, recvListPtr % nlist, &
                                     recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, arrayIn and arrayOut dims must match.'
         call dmpar_abort(dminfo)
      else
         arrayOut(:) = arrayIn(:)
      end if
#endif

   end subroutine dmpar_alltoall_field1dReal


   subroutine dmpar_alltoall_field2dReal(dminfo, arrayIn, arrayOut, dim1, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, nOwnedList, nNeededList
      real (kind=RKIND), dimension(dim1,nOwnedList), intent(in) :: arrayIn
      real (kind=RKIND), dimension(dim1,nNeededList), intent(inout) :: arrayOut
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i, d2

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(:,recvListPtr % list(i)) = arrayIn(:,sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * recvListPtr % nlist
            allocate(recvListPtr % rbuffer(d2))
            call MPI_Irecv(recvListPtr % rbuffer, d2, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d2 = dim1 * sendListPtr % nlist
            allocate(sendListPtr % rbuffer(d2))
            call packSendBuf2dReal(1, dim1, nOwnedList, arrayIn, sendListPtr, 1, d2, &
                                   sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, d2, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d2 = dim1 * recvListPtr % nlist
            call unpackRecvBuf2dReal(1, dim1, nNeededList, arrayOut, recvListPtr, 1, d2, &
                                     recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, arrayIn and arrayOut dims must match.'
         call dmpar_abort(dminfo)
      else
         arrayOut(:,:) = arrayIn(:,:)
      end if
#endif

   end subroutine dmpar_alltoall_field2dReal

  
   subroutine dmpar_alltoall_field3dReal(dminfo, arrayIn, arrayOut, dim1, dim2, nOwnedList, nNeededList, sendList, recvList)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: dim1, dim2, nOwnedList, nNeededList
      real (kind=RKIND), dimension(dim1,dim2,nOwnedList), intent(in) :: arrayIn
      real (kind=RKIND), dimension(dim1,dim2,nNeededList), intent(inout) :: arrayOut
      type (exchange_list), pointer :: sendList, recvList

      type (exchange_list), pointer :: sendListPtr, recvListPtr
      integer :: lastPackedIdx, lastUnpackedIdx, nPacked, nUnpacked
      integer :: mpi_ierr
      integer :: i, d3

#ifdef _MPI

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID == dminfo % my_proc_id) exit
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID == dminfo % my_proc_id) exit
         recvListPtr => recvListPtr % next
      end do

      if (associated(recvListPtr) .and. associated(sendListPtr)) then
         do i=1,recvListPtr % nlist
            arrayOut(:,:,recvListPtr % list(i)) = arrayIn(:,:,sendListPtr % list(i))
         end do
      end if

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * recvListPtr % nlist
            allocate(recvListPtr % rbuffer(d3))
            call MPI_Irecv(recvListPtr % rbuffer, d3, MPI_REALKIND, &
                           recvListPtr % procID, recvListPtr % procID, dminfo % comm, recvListPtr % reqID, mpi_ierr)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            d3 = dim1 * dim2 * sendListPtr % nlist
            allocate(sendListPtr % rbuffer(d3))
            call packSendBuf3dReal(1, dim1, 1, dim2, nOwnedList, arrayIn, sendListPtr, 1, d3, &
                                   sendListPtr % rbuffer, nPacked, lastPackedIdx)
            call MPI_Isend(sendListPtr % rbuffer, d3, MPI_REALKIND, &
                           sendListPtr % procID, dminfo % my_proc_id, dminfo % comm, sendListPtr % reqID, mpi_ierr)
         end if
         sendListPtr => sendListPtr % next
      end do

      recvListPtr => recvList
      do while (associated(recvListPtr))
         if (recvListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(recvListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            d3 = dim1 * dim2 * recvListPtr % nlist
            call unpackRecvBuf3dReal(1, dim1, 1, dim2, nNeededList, arrayOut, recvListPtr, 1, d3, &
                                     recvListPtr % rbuffer, nUnpacked, lastUnpackedIdx)
            deallocate(recvListPtr % rbuffer)
         end if
         recvListPtr => recvListPtr % next
      end do

      sendListPtr => sendList
      do while (associated(sendListPtr))
         if (sendListPtr % procID /= dminfo % my_proc_id) then
            call MPI_Wait(sendListPtr % reqID, MPI_STATUS_IGNORE, mpi_ierr)
            deallocate(sendListPtr % rbuffer)
         end if
         sendListPtr => sendListPtr % next
      end do

#else
      if (nOwnedList /= nNeededList) then
         write(0,*) 'Error in dmpar_alltoall_field: For non-dmpar, arrayIn and arrayOut dims must match.'
         call dmpar_abort(dminfo)
      else
         arrayOut(:,:,:) = arrayIn(:,:,:)
      end if
#endif

   end subroutine dmpar_alltoall_field3dReal

  
   subroutine packSendBuf1dInteger(nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startPackIdx
      integer, dimension(nField), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      integer, dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + 1
         if (nPacked > nBuffer) then
            nPacked = nPacked - 1
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked) = field(sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine packSendBuf1dInteger


   subroutine packSendBuf2dInteger(ds, de, nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startPackIdx
      integer, dimension(ds:de,1:nField), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      integer, dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i, n

      n = de-ds+1

      if (n > nBuffer) then
         write(0,*) 'packSendBuf2dInteger: Not enough space in buffer to fit a single slice.'
         return
      end if

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + n
         if (nPacked > nBuffer) then
            nPacked = nPacked - n
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked-n+1:nPacked) = field(ds:de,sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine packSendBuf2dInteger


   subroutine packSendBuf1dReal(nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startPackIdx
      real (kind=RKIND), dimension(nField), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      real (kind=RKIND), dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + 1
         if (nPacked > nBuffer) then
            nPacked = nPacked - 1
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked) = field(sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine packSendBuf1dReal


   subroutine packSendBuf2dReal(ds, de, nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startPackIdx
      real (kind=RKIND), dimension(ds:de,1:nField), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      real (kind=RKIND), dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i, n

      n = de-ds+1

      if (n > nBuffer) then
         write(0,*) 'packSendBuf2dReal: Not enough space in buffer to fit a single slice.'
         return
      end if

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + n
         if (nPacked > nBuffer) then
            nPacked = nPacked - n
            lastPackedIdx = i - 1
            return
         end if
         buffer(nPacked-n+1:nPacked) = field(ds:de,sendList % list(i))
      end do
      lastPackedIdx = sendList % nlist

   end subroutine packSendBuf2dReal


   subroutine packSendBuf3dReal(d1s, d1e, d2s, d2e, nField, field, sendList, startPackIdx, nBuffer, buffer, nPacked, lastPackedIdx)

      implicit none

      integer, intent(in) :: d1s, d1e, d2s, d2e, nField, nBuffer, startPackIdx
      real (kind=RKIND), dimension(d1s:d1e,d2s:d2e,1:nField), intent(in) :: field
      type (exchange_list), intent(in) :: sendList
      real (kind=RKIND), dimension(nBuffer), intent(out) :: buffer
      integer, intent(inout) :: nPacked, lastPackedIdx

      integer :: i, j, k, n

      n = (d1e-d1s+1) * (d2e-d2s+1)

      if (n > nBuffer) then
         write(0,*) 'packSendBuf2dReal: Not enough space in buffer to fit a single slice.'
         return
      end if

      nPacked = 0
      do i=startPackIdx, sendList % nlist
         nPacked = nPacked + n
         if (nPacked > nBuffer) then
            nPacked = nPacked - n
            lastPackedIdx = i - 1
            return
         end if
         k = nPacked-n+1
         do j=d2s,d2e
            buffer(k:k+d1e-d1s) = field(d1s:d1e,j,sendList % list(i))
            k = k + d1e-d1s+1
         end do
      end do
      lastPackedIdx = sendList % nlist

   end subroutine packSendBuf3dReal


   subroutine unpackRecvBuf1dInteger(nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startUnpackIdx
      integer, dimension(nField), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      integer, dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + 1
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - 1
            lastUnpackedIdx = i - 1
            return
         end if
         field(recvList % list(i)) = buffer(nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine unpackRecvBuf1dInteger


   subroutine unpackRecvBuf2dInteger(ds, de, nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startUnpackIdx
      integer, dimension(ds:de,1:nField), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      integer, dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i, n

      n = de-ds+1

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + n
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - n
            lastUnpackedIdx = i - 1
            return
         end if
         field(ds:de,recvList % list(i)) = buffer(nUnpacked-n+1:nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine unpackRecvBuf2dInteger


   subroutine unpackRecvBuf1dReal(nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: nField, nBuffer, startUnpackIdx
      real (kind=RKIND), dimension(nField), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      real (kind=RKIND), dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + 1
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - 1
            lastUnpackedIdx = i - 1
            return
         end if
         field(recvList % list(i)) = buffer(nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine unpackRecvBuf1dReal


   subroutine unpackRecvBuf2dReal(ds, de, nField, field, recvList, startUnpackIdx, nBuffer, buffer, nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: ds, de, nField, nBuffer, startUnpackIdx
      real (kind=RKIND), dimension(ds:de,1:nField), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      real (kind=RKIND), dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i, n

      n = de-ds+1

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + n
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - n
            lastUnpackedIdx = i - 1
            return
         end if
         field(ds:de,recvList % list(i)) = buffer(nUnpacked-n+1:nUnpacked)
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine unpackRecvBuf2dReal


   subroutine unpackRecvBuf3dReal(d1s, d1e, d2s, d2e, nField, field, recvList, startUnpackIdx, nBuffer, buffer, &
                                  nUnpacked, lastUnpackedIdx)

      implicit none

      integer, intent(in) :: d1s, d1e, d2s, d2e, nField, nBuffer, startUnpackIdx
      real (kind=RKIND), dimension(d1s:d1e,d2s:d2e,1:nField), intent(inout) :: field
      type (exchange_list), intent(in) :: recvList
      real (kind=RKIND), dimension(nBuffer), intent(in) :: buffer
      integer, intent(inout) :: nUnpacked, lastUnpackedIdx

      integer :: i, j, k, n

      n = (d1e-d1s+1) * (d2e-d2s+1)

      nUnpacked = 0
      do i=startUnpackIdx, recvList % nlist
         nUnpacked = nUnpacked + n
         if (nUnpacked > nBuffer) then
            nUnpacked = nUnpacked - n
            lastUnpackedIdx = i - 1
            return
         end if
         k = nUnpacked-n+1
         do j=d2s,d2e
            field(d1s:d1e,j,recvList % list(i)) = buffer(k:k+d1e-d1s)
            k = k + d1e-d1s+1
         end do
      end do
      lastUnpackedIdx = recvList % nlist

   end subroutine unpackRecvBuf3dReal


   integer function binary_search(array, d1, n1, n2, key)

      implicit none

      integer, intent(in) :: d1, n1, n2, key
      integer, dimension(d1,n1:n2), intent(in) :: array

      integer :: l, u, k

      binary_search = n2+1

      l = n1
      u = n2
      k = (l+u)/2
      do while (u >= l)
         if (array(1,k) == key) then
            binary_search = k
            exit
         else if (array(1,k) < key) then
            l = k + 1
            k = (l+u)/2
         else
            u = k - 1
            k = (l+u)/2
         end if
      end do

   end function binary_search


   subroutine quicksort(nArray, array)

      implicit none

      integer, intent(in) :: nArray
      integer, dimension(2,nArray), intent(inout) :: array

      integer :: i, j, top, l, r, pivot, s
      integer :: pivot_value
      integer, dimension(2) :: temp
      integer, dimension(1000) :: lstack, rstack

      if (nArray < 1) return

      top = 1
      lstack(top) = 1
      rstack(top) = nArray

      do while (top > 0)

         l = lstack(top)
         r = rstack(top)
         top = top - 1

         pivot = (l+r)/2

         pivot_value = array(1,pivot)
         temp(:) = array(:,pivot)
         array(:,pivot) = array(:,r)
         array(:,r) = temp(:)

         s = l
         do i=l,r-1
            if (array(1,i) <= pivot_value) then
               temp(:) = array(:,s)
               array(:,s) = array(:,i)
               array(:,i) = temp(:)
               s = s + 1
            end if
         end do

         temp(:) = array(:,s)
         array(:,s) = array(:,r)
         array(:,r) = temp(:)

         if (s-1 > l) then
            top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
            lstack(top) = l
            rstack(top) = s-1
         end if

         if (r > s+1) then
            top = top + 1
if (top > 1000) write(0,*) 'Error: Quicksort exhausted its stack.'
            lstack(top) = s+1
            rstack(top) = r
         end if
      end do

   end subroutine quicksort


   subroutine compare_r1(real1DField, checkpoint)

      implicit none

      type (field1DReal), pointer :: real1DField
      character(len=*), intent(in) :: checkpoint

      real (kind=RKIND), dimension(:), pointer :: arrayOut_r1, arrayComp_r1
      integer :: i, outdim, outdimglobal
      type (exchange_list), pointer :: sendList, recvList


      if (trim(real1DField % dimnames(1)) == 'nCells') then
         outdim = nCellsSolve
         outdimglobal = nCellsGlobal
         sendList => cellSendList
         recvList => cellRecvList
      else if (trim(real1DField % dimnames(1)) == 'nEdges') then
         outdim = nEdgesSolve
         outdimglobal = nEdgesGlobal
         sendList => edgeSendList
         recvList => edgeRecvList
      else if (trim(real1DField % dimnames(1)) == 'nVertices') then
         outdim = nVerticesSolve
         outdimglobal = nVerticesGlobal
         sendList => vertexSendList
         recvList => vertexRecvList
      else
         nullify(sendList)
         nullify(recvList)
      end if

      if (dminfo % my_proc_id == 0) then
         if (associated(sendList)) then

         allocate(arrayOut_r1(outdimglobal))
         allocate(arrayComp_r1(outdimglobal))
         call dmpar_alltoall_field1dReal(dminfo, real1DField % array(1:outdim), arrayOut_r1, &
                             outdim, outdimglobal, sendList, recvList)
         if (write_fields) then
            write(91) arrayOut_r1
         else
            read(91) arrayComp_r1

            OUT1D: do i=1,size(arrayComp_r1)
               if (arrayComp_r1(i) /= arrayOut_r1(i)) then
                  write(stderrUnit,*) trim(checkpoint)//' :: Field '//trim(real1DField % fieldName)//' differs'
                  exit OUT1D
               end if
            end do OUT1D
         end if
         deallocate(arrayOut_r1)
         deallocate(arrayComp_r1)

         end if
      else
         if (associated(sendList)) then

         allocate(arrayOut_r1(0))
         call dmpar_alltoall_field1dReal(dminfo, real1DField % array(1:outdim), arrayOut_r1, &
                             outdim, 0, sendList, recvList)
         deallocate(arrayOut_r1)

         end if
      end if

   end subroutine compare_r1


   subroutine compare_r2(real2DField, checkpoint)

      implicit none

      type (field2DReal), pointer :: real2DField
      character(len=*), intent(in) :: checkpoint

      real (kind=RKIND), dimension(:,:), pointer :: arrayOut_r2, arrayComp_r2
      integer :: i, j, outdim, outdimglobal
      type (exchange_list), pointer :: sendList, recvList

      if (trim(real2DField % dimnames(2)) == 'nCells') then
         outdim = nCellsSolve
         outdimglobal = nCellsGlobal
         sendList => cellSendList
         recvList => cellRecvList
      else if (trim(real2DField % dimnames(2)) == 'nEdges') then
         outdim = nEdgesSolve
         outdimglobal = nEdgesGlobal
         sendList => edgeSendList
         recvList => edgeRecvList
      else if (trim(real2DField % dimnames(2)) == 'nVertices') then
         outdim = nVerticesSolve
         outdimglobal = nVerticesGlobal
         sendList => vertexSendList
         recvList => vertexRecvList
      else
         nullify(sendList)
         nullify(recvList)
      end if

      if (dminfo % my_proc_id == 0) then
         if (associated(sendList)) then

         allocate(arrayOut_r2(real2DField % dimsizes(1),outdimglobal))
         allocate(arrayComp_r2(real2DField % dimsizes(1),outdimglobal))
         call dmpar_alltoall_field2dReal(dminfo, real2DField % array(:,1:outdim), arrayOut_r2, &
                             real2DField % dimsizes(1), outdim, outdimglobal, sendList, recvList)
         if (write_fields) then
            write(91) arrayOut_r2
         else
            read(91) arrayComp_r2

            OUT2D: do j=1,size(arrayComp_r2,2)
            do i=1,size(arrayComp_r2,1)
               if (arrayComp_r2(i,j) /= arrayOut_r2(i,j)) then
                  write(stderrUnit,*) trim(checkpoint)//' :: Field '//trim(real2DField % fieldName)//' differs'
                  exit OUT2D
               end if
            end do
            end do OUT2D
         end if
         deallocate(arrayOut_r2)
         deallocate(arrayComp_r2)

         end if
      else
         if (associated(sendList)) then

         allocate(arrayOut_r2(0,0))
         call dmpar_alltoall_field2dReal(dminfo, real2DField % array(:,1:outdim), arrayOut_r2, &
                             real2DField % dimsizes(1), outdim, 0, sendList, recvList)
         deallocate(arrayOut_r2)

         end if
      end if

   end subroutine compare_r2


   subroutine compare_r3(real3DField, checkpoint)

      implicit none

      type (field3DReal), pointer :: real3DField
      character(len=*), intent(in) :: checkpoint

      real (kind=RKIND), dimension(:,:,:), pointer :: arrayOut_r3, arrayComp_r3
      integer :: i, j, k, outdim, outdimglobal
      type (exchange_list), pointer :: sendList, recvList


      if (trim(real3DField % dimnames(3)) == 'nCells') then
         outdim = nCellsSolve
         outdimglobal = nCellsGlobal
         sendList => cellSendList
         recvList => cellRecvList
      else if (trim(real3DField % dimnames(3)) == 'nEdges') then
         outdim = nEdgesSolve
         outdimglobal = nEdgesGlobal
         sendList => edgeSendList
         recvList => edgeRecvList
      else if (trim(real3DField % dimnames(3)) == 'nVertices') then
         outdim = nVerticesSolve
         outdimglobal = nVerticesGlobal
         sendList => vertexSendList
         recvList => vertexRecvList
      else
         nullify(sendList)
         nullify(recvList)
      end if

      if (dminfo % my_proc_id == 0) then
         if (associated(sendList)) then

         allocate(arrayOut_r3(real3DField % dimsizes(1),real3DField % dimsizes(2),outdimglobal))
         allocate(arrayComp_r3(real3DField % dimsizes(1),real3DField % dimsizes(2),outdimglobal))
         call dmpar_alltoall_field3dReal(dminfo, real3DField % array(:,:,1:outdim), arrayOut_r3, &
                             real3DField % dimsizes(1), real3DField % dimsizes(2), outdim, outdimglobal, sendList, recvList)
         if (write_fields) then
            write(91) arrayOut_r3
         else
            read(91) arrayComp_r3

            OUT3D: do k=1,size(arrayComp_r3,3)
            do j=1,size(arrayComp_r3,2)
            do i=1,size(arrayComp_r3,1)
               if (arrayComp_r3(i,j,k) /= arrayOut_r3(i,j,k)) then
                  write(stderrUnit,*) trim(checkpoint)//' :: Field '//trim(real3DField % fieldName)//' differs'
                  exit OUT3D
               end if
            end do
            end do
            end do OUT3D
         end if
         deallocate(arrayOut_r3)
         deallocate(arrayComp_r3)

         end if
      else
         if (associated(sendList)) then

         allocate(arrayOut_r3(0,0,0))
         call dmpar_alltoall_field3dReal(dminfo, real3DField % array(:,:,1:outdim), arrayOut_r3, &
                             real3DField % dimsizes(1), real3DField % dimsizes(2), outdim, 0, sendList, recvList)
         deallocate(arrayOut_r3)

         end if
      end if

   end subroutine compare_r3


   subroutine compare_i1(int1DField, checkpoint)

      implicit none

      type (field1DInteger), pointer :: int1DField
      character(len=*), intent(in) :: checkpoint

      integer, dimension(:), pointer :: arrayOut_i1, arrayComp_i1
      integer :: i, outdim, outdimglobal
      type (exchange_list), pointer :: sendList, recvList


      if (trim(int1DField % dimnames(1)) == 'nCells') then
         outdim = nCellsSolve
         outdimglobal = nCellsGlobal
         sendList => cellSendList
         recvList => cellRecvList
      else if (trim(int1DField % dimnames(1)) == 'nEdges') then
         outdim = nEdgesSolve
         outdimglobal = nEdgesGlobal
         sendList => edgeSendList
         recvList => edgeRecvList
      else if (trim(int1DField % dimnames(1)) == 'nVertices') then
         outdim = nVerticesSolve
         outdimglobal = nVerticesGlobal
         sendList => vertexSendList
         recvList => vertexRecvList
      else
         nullify(sendList)
         nullify(recvList)
      end if

      if (dminfo % my_proc_id == 0) then
         if (associated(sendList)) then

         allocate(arrayOut_i1(outdimglobal))
         allocate(arrayComp_i1(outdimglobal))
         call dmpar_alltoall_field1dInteger(dminfo, int1DField % array(1:outdim), arrayOut_i1, &
                             outdim, outdimglobal, sendList, recvList)
         if (write_fields) then
            write(91) arrayOut_i1
         else
            read(91) arrayComp_i1

            OUT1DI: do i=1,size(arrayComp_i1)
               if (arrayComp_i1(i) /= arrayOut_i1(i)) then
                  write(stderrUnit,*) trim(checkpoint)//' :: Field '//trim(int1DField % fieldName)//' differs'
                  exit OUT1DI
               end if
            end do OUT1DI
         end if
         deallocate(arrayOut_i1)
         deallocate(arrayComp_i1)

         end if
      else
         if (associated(sendList)) then

         allocate(arrayOut_i1(0))
         call dmpar_alltoall_field1dInteger(dminfo, int1DField % array(1:outdim), arrayOut_i1, &
                             outdim, 0, sendList, recvList)
         deallocate(arrayOut_i1)

         end if
      end if

   end subroutine compare_i1


   subroutine compare_i2(int2DField, checkpoint)

      implicit none

      type (field2DInteger), pointer :: int2DField
      character(len=*), intent(in) :: checkpoint

      integer, dimension(:,:), pointer :: arrayOut_i2, arrayComp_i2
      integer :: i, j, outdim, outdimglobal
      type (exchange_list), pointer :: sendList, recvList


      if (trim(int2DField % dimnames(2)) == 'nCells') then
         outdim = nCellsSolve
         outdimglobal = nCellsGlobal
         sendList => cellSendList
         recvList => cellRecvList
      else if (trim(int2DField % dimnames(2)) == 'nEdges') then
         outdim = nEdgesSolve
         outdimglobal = nEdgesGlobal
         sendList => edgeSendList
         recvList => edgeRecvList
      else if (trim(int2DField % dimnames(2)) == 'nVertices') then
         outdim = nVerticesSolve
         outdimglobal = nVerticesGlobal
         sendList => vertexSendList
         recvList => vertexRecvList
      else
         nullify(sendList)
         nullify(recvList)
      end if

      if (dminfo % my_proc_id == 0) then
         if (associated(sendList)) then

         allocate(arrayOut_i2(int2DField % dimsizes(1),outdimglobal))
         allocate(arrayComp_i2(int2DField % dimsizes(1),outdimglobal))
         call dmpar_alltoall_field2dInteger(dminfo, int2DField % array(:,1:outdim), arrayOut_i2, &
                             int2DField % dimsizes(1), outdim, outdimglobal, sendList, recvList)
         if (write_fields) then
            write(91) arrayOut_i2
         else
            read(91) arrayComp_i2

            OUT2DI: do j=1,size(arrayComp_i2,2)
            do i=1,size(arrayComp_i2,1)
               if (arrayComp_i2(i,j) /= arrayOut_i2(i,j)) then
                  write(stderrUnit,*) trim(checkpoint)//' :: Field '//trim(int2DField % fieldName)//' differs'
                  exit OUT2DI
               end if
            end do
            end do OUT2DI
         end if
         deallocate(arrayOut_i2)
         deallocate(arrayComp_i2)

         end if
      else
         if (associated(sendList)) then

         allocate(arrayOut_i2(0,0))
         call dmpar_alltoall_field2dInteger(dminfo, int2DField % array(:,1:outdim), arrayOut_i2, &
                             int2DField % dimsizes(1), outdim, 0, sendList, recvList)
         deallocate(arrayOut_i2)

         end if
      end if

   end subroutine compare_i2

end module mpas_debugging
