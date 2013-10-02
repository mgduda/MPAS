      module layout1
!     use resol_def
      implicit none
      save
cc
      integer           nodes, nodes_comp,nodes_io,
     x                  me,
     x                  ls_dim,
     x                  ls_max_node,
     x                  lats_dim_a,
     x                  lats_dim_r,
!    x                  lats_dim_ext,
     x                  lats_node_a,
     x                  lats_node_a_max,
!    x                  lats_node_ext,
     x                  lats_node_r,
     x                  lats_node_r_max,
     x                  ipt_lats_node_a,
!    x                  ipt_lats_node_ext,
     x                  ipt_lats_node_r,
     x                  len_trie_ls,
     x                  len_trio_ls,
     x                  len_trie_ls_max,
     x                  len_trio_ls_max,
     x                  lon_dim_a,
     x                  lon_dim_r,
     x                  me_l_0
cc
      INTEGER ,ALLOCATABLE :: lat1s_a(:),lat1s_r(:)
!mjr .  lon_dims_a(:),lon_dims_ext(:),lon_dims_r(:)
!mjr .  lon_dims_a(:),                lon_dims_r(:)

!hmhj ndsl
      integer   lonfull,lonhalf,lonpart,lonlenmax,mylonlen
      integer   latfull,lathalf,latpart,latlenmax,mylatlen

      integer, allocatable :: lonstr(:),lonlen(:)
      integer, allocatable :: latstr(:),latlen(:)
      real, allocatable :: cosglat(:)
      real, allocatable :: gglat(:),gglati(:),gslati(:),ggfact(:,:)
      real, allocatable :: gglon(:),ggloni(:)

      end module layout1
