module mpas_debugging

   use mpas_grid_types

   contains

   subroutine get_indices(ioffset, dims, indices)

      implicit none

      integer, intent(in) :: ioffset
      integer, dimension(:), intent(in) :: dims
      integer, dimension(:), intent(inout) :: indices

      integer :: i, j, m

      ! Convert offset to 0-based indexing
      i = ioffset - 1

      ! Get indices in 0-based indexing
      indices(1) = mod(i,dims(1))
      i = (i - indices(1))/dims(1)

      do j=2,size(dims)
         indices(j) = mod(i,dims(j))
         i = (i - indices(j))/dims(j)
      end do

      ! Convert to 1-based indexing
      indices(:) = indices(:) + 1

   end subroutine get_indices


#include "nan_checks.inc"


end module mpas_debugging
