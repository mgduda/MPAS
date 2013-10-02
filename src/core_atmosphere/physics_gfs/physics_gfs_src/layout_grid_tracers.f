      module layout_grid_tracers
      use machine, only : kind_evod
      implicit none
      save
cc
      integer xhalo
      integer yhalo
cc
      real(kind=kind_evod) ,allocatable :: rg1_h(:,:,:)
      real(kind=kind_evod) ,allocatable :: rg2_h(:,:,:)
      real(kind=kind_evod) ,allocatable :: rg3_h(:,:,:)
cc
      real(kind=kind_evod) ,allocatable :: rg1_a(:,:,:)
      real(kind=kind_evod) ,allocatable :: rg2_a(:,:,:)
      real(kind=kind_evod) ,allocatable :: rg3_a(:,:,:)
cc
      end module layout_grid_tracers
