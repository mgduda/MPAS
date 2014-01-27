      module gfsmisc_def
      use machine

      implicit none
      save

      logical, parameter    :: fix_ncld_hr=.true.
      integer, parameter    :: nrcmax=30 ! Maximum number of random clouds per 1200s
      real,    parameter    :: rlapse=0.65e-2, omz1=10.0

      real (kind=kind_ior)  :: CONS0, CONS0P5, CONS1200, CONS3600
      real (kind=kind_phys) :: pdryini
      real (kind=kind_evod) :: phour,zhour,deltim
      real (kind=kind_evod) :: slag,sdec,cdec,batah
      integer               :: kdt,iret,iprint,MAXSTP
      logical               :: lsout

#if !defined(mpas)
      integer,  allocatable :: lonsperlar(:),global_lats_r(:)
      integer,  allocatable :: jindx1(:),jindx2(:)
      real,     allocatable :: ozplin(:,:,:,:)     !OZONE PL Coeff
      real (kind=kind_phys), allocatable ::  xlat(:,:), xlon(:,:),
     &     coszdg(:,:), hprime(:,:,:), fluxr(:,:,:), sfalb(:,:), 
     &     swh(:,:,:), hlw(:,:,:),
     &     sinlat_r2(:,:),coslat_r2(:,:)
      real (kind=kind_phys), allocatable :: phy_f3d(:,:,:,:), 
     &     phy_f2d(:,:,:), ddy(:), fscav(:)
#endif
      real (kind=kind_phys), allocatable :: fscav(:)
      end module gfsmisc_def
