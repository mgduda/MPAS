      integer function num_parthds()
      use omp_lib
!!$OMP PARALLEL
!     num_parthds = omp_get_num_threads()
      num_parthds = 6
!     num_parthds = 4
!!$OMP END PARALLEL
      return
      end

 
