module mpas_debugging

   use mpas_grid_types

   logical :: write_fields
   type (mpas_pool_type), pointer :: all_fields 


   contains


   subroutine mpas_debugging_init(arg_write_fields, field_pool)

      implicit none

      logical, intent(in) :: arg_write_fields
      type (mpas_pool_type), pointer :: field_pool

return
      write_fields = arg_write_fields
      all_fields => field_pool

      if (write_fields) then
         open(91,file='checkpoints.dat',form='unformatted',status='replace',action='write')
      else
         open(91,file='checkpoints.dat',form='unformatted',status='old',action='read')
      end if

   end subroutine mpas_debugging_init


   subroutine mpas_debugging_check(checkpoint)

        implicit none

        character(len=*), intent(in) :: checkpoint

        integer :: i, j, k, outdim
        type (mpas_pool_iterator_type) :: poolItr
        type (field2DReal), pointer :: real2DField, real2DField_comp
        type (field3DReal), pointer :: real3DField, real3DField_comp

return
        call mpas_pool_begin_iteration(all_fields)
        do while (mpas_pool_get_next_member(all_fields, poolItr))
            if (poolItr % memberType == MPAS_POOL_FIELD) then
                if (poolItr % dataType == MPAS_POOL_REAL) then
                    if (poolItr % nDims == 2) then
                        call mpas_pool_get_field(all_fields, poolItr % memberName, real2DField)
                        if (write_fields) then
                           if (trim(real2DField % dimnames(2)) == 'nCells' .or. &
                               trim(real2DField % dimnames(2)) == 'nEdges' .or. &
                               trim(real2DField % dimnames(2)) == 'nVertices') then
                              outdim = real2DField % dimsizes(2)-1
                              write(91) real2DField % array(:,1:real2DField % dimsizes(2)-1)
                           else
                              outdim = real2DField % dimsizes(2)
                              write(91) real2DField % array(:,:)
                           end if
                        else
                           call mpas_duplicate_field(real2DField, real2DField_comp)
                           if (trim(real2DField % dimnames(2)) == 'nCells' .or. &
                               trim(real2DField % dimnames(2)) == 'nEdges' .or. &
                               trim(real2DField % dimnames(2)) == 'nVertices') then
                              outdim = real2DField % dimsizes(2)-1
                              read(91) real2DField_comp % array(:,1:real2DField % dimsizes(2)-1)
                           else
                              outdim = real2DField % dimsizes(2)
                              read(91) real2DField_comp % array(:,:)
                           end if

                           OUT2D: do j=1,outdim
                           do i=1,real2DField % dimsizes(1)
                              if (real2DField % array(i,j) /= real2DField_comp % array(i,j)) then
                                 write(stderrUnit,*) trim(checkpoint)//' :: Field '//trim(real2DField % fieldName)//' differs'
                                 exit OUT2D
                              end if
                           end do
                           end do OUT2D

                           call mpas_deallocate_field(real2DField_comp)
                        end if
                    else if (poolItr % nDims == 3) then
                        call mpas_pool_get_field(all_fields, poolItr % memberName, real3DField)
                        if (write_fields) then
                           if (trim(real3DField % dimnames(3)) == 'nCells' .or. &
                               trim(real3DField % dimnames(3)) == 'nEdges' .or. &
                               trim(real3DField % dimnames(3)) == 'nVertices') then
                              outdim = real3DField % dimsizes(3)-1
                              write(91) real3DField % array(:,:,1:real3DField % dimsizes(3)-1)
                           else
                              outdim = real3DField % dimsizes(3)
                              write(91) real3DField % array(:,:,:)
                           end if
                        else
                           call mpas_duplicate_field(real3DField, real3DField_comp)
                           if (trim(real3DField % dimnames(3)) == 'nCells' .or. &
                               trim(real3DField % dimnames(3)) == 'nEdges' .or. &
                               trim(real3DField % dimnames(3)) == 'nVertices') then
                              outdim = real3DField % dimsizes(3)-1
                              read(91) real3DField_comp % array(:,:,1:real3DField % dimsizes(3)-1)
                           else
                              outdim = real3DField % dimsizes(3)
                              read(91) real3DField_comp % array(:,:,:)
                           end if

                           OUT3D: do k=1,outdim
                           do j=1,real3DField % dimsizes(2)
                           do i=1,real3DField % dimsizes(1)
                              if (real3DField % array(i,j,k) /= real3DField_comp % array(i,j,k)) then
                                 write(stderrUnit,*) trim(checkpoint)//' :: Field '//trim(real3DField % fieldName)//' differs'
                                 exit OUT3D
                              end if
                           end do
                           end do
                           end do OUT3D

                           call mpas_deallocate_field(real3DField_comp)
                        end if
                    end if
                end if
            end if
        end do

   end subroutine mpas_debugging_check


   subroutine mpas_debugging_finalize()

      implicit none

return
      close(91)

   end subroutine mpas_debugging_finalize

end module mpas_debugging
