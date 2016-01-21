! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_tbd_ddt
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_vars
 use machine
 use nuopc_physics

 implicit none
 private
 public:: mpas_tbd_ddt

 type,extends(tbd_ddt):: mpas_tbd_ddt

    real(kind=kind_phys),dimension(:),pointer    :: dpshc_l   => null()   
    real(kind=kind_phys),dimension(:,:),pointer  :: slc_l     => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: smc_l     => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: stc_l     => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: updmf_l   => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: dwnmf_l   => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: detmf_l   => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: rann_l    => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: phyf2d_l  => null()
    real(kind=kind_phys),dimension(:,:,:),pointer:: phyf3d_l  => null()
    real(kind=kind_phys),dimension(:,:,:),pointer:: ozmix_l   => null()

    !... new pointers added in nuopc_physics.F90 version 2:
    real(kind=kind_phys),dimension(:),pointer    :: fscav_l   => null()
    real(kind=kind_phys),dimension(:),pointer    :: fswtr_l   => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: phyfctd_l => null()

    contains
       procedure:: mpas_set       => mpas_tbd_set
       procedure:: mpas_write     => mpas_tbd_write
       procedure:: allocate_mem   => mpas_tbd_ddt_allocate
       procedure:: deallocate_mem => mpas_tbd_ddt_deallocate 

 end type

!modifications:
!--------------
! * added the local arrays fscav_l and fswtr_l for tracer scavenging, and phyfctd_l for the Chikira-Sugiyama
!   convection scheme with the release of nuopc_physics.F90 version 2. 
!   Laura D. Fowler (laura@ucar.edu) / 2015-09-22.
! * temporarily fixed the array rann_l to 0.6 for debugging.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-19.
! * added deallocation of local array ozmix_l.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-27.
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.
! * in subroutine mpas_tbd_set, now fills the local arrays phyf2d_l and phyf3d_l with the one-dimensional arrays
!   gscond_psp, gscond_psp1, gscond_tp, gscond_qp, gscond_tp1, and gscond_qp1 so that they can be saved in the
!   restart stream for high-resolution meshes. in subroutine mpas_tbd_write, now saves the arrays phy_f2d and
!   phy_f3d between time-steps.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-21.


 contains


!=================================================================================================================
 subroutine mpas_tbd_ddt_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_tbd_ddt),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nSoilLevels,nOznLevels,nOznPlCoeff
 integer:: ncld,nrcm,ntrac,nctp,num_p2d,num_p3d

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_tbd_ddt_allocate:'

 ncld    = mdl_parameters%ncld
 nrcm    = mdl_parameters%nrcm
 ntrac   = mdl_parameters%ntrac
 nctp    = mdl_parameters%nctp
 num_p2d = mdl_parameters%num_p2d
 num_p3d = mdl_parameters%num_p3d

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)

 if(.not.associated(this%dpshc_l)  ) allocate(this%dpshc_l(1:nCellsSolve)              )
 if(.not.associated(this%slc_l)    ) allocate(this%slc_l(1:nCellsSolve,1:nSoilLevels)  )
 if(.not.associated(this%smc_l)    ) allocate(this%smc_l(1:nCellsSolve,1:nSoilLevels)  )
 if(.not.associated(this%stc_l)    ) allocate(this%stc_l(1:nCellsSolve,1:nSoilLevels)  )
 if(.not.associated(this%updmf_l)  ) allocate(this%updmf_l(1:nCellsSolve,1:nVertLevels))
 if(.not.associated(this%dwnmf_l)  ) allocate(this%dwnmf_l(1:nCellsSolve,1:nVertLevels))
 if(.not.associated(this%detmf_l)  ) allocate(this%detmf_l(1:nCellsSolve,1:nVertLevels))

 if(.not.associated(this%rann_l)   ) allocate(this%rann_l(1:nCellsSolve,nrcm)                   )
 if(.not.associated(this%phyf2d_l) ) allocate(this%phyf2d_l(1:nCellsSolve,num_p2d)              )
 if(.not.associated(this%phyf3d_l) ) allocate(this%phyf3d_l(1:nCellsSolve,1:nVertLevels,num_p3d))

 if(.not.associated(this%ozmix_l)  ) allocate(this%ozmix_l(1:nCellsSolve,1:nOznLevels,1:nOznPlCoeff))

!... new pointers added in nuopc_physics.F90 version 2:
 if(.not.associated(this%fscav_l)  ) allocate(this%fscav_l(ntrac-ncld+2)        )
 if(.not.associated(this%fswtr_l)  ) allocate(this%fswtr_l(ntrac-ncld+2)        )
 if(.not.associated(this%phyfctd_l)) allocate(this%phyfctd_l(1:nCellsSolve,nctp))

!write(0,*) '--- end subroutine mpas_tbd_ddt_allocate'

 end subroutine mpas_tbd_ddt_allocate

!=================================================================================================================
 subroutine mpas_tbd_ddt_deallocate(this)
!=================================================================================================================

 class(mpas_tbd_ddt),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_tbd_ddt_deallocate:'

 if(associated(this%dpshc_l)  ) deallocate(this%dpshc_l  )
 if(associated(this%slc_l)    ) deallocate(this%slc_l    )
 if(associated(this%smc_l)    ) deallocate(this%smc_l    )
 if(associated(this%stc_l)    ) deallocate(this%stc_l    )
 if(associated(this%updmf_l)  ) deallocate(this%updmf_l  )
 if(associated(this%dwnmf_l)  ) deallocate(this%dwnmf_l  )
 if(associated(this%detmf_l)  ) deallocate(this%detmf_l  )
 if(associated(this%rann_l)   ) deallocate(this%rann_l   )
 if(associated(this%phyf2d_l) ) deallocate(this%phyf2d_l )
 if(associated(this%phyf3d_l) ) deallocate(this%phyf3d_l )

 if(associated(this%ozmix_l)  ) deallocate(this%ozmix_l  )

!... new pointers added in nuopc_physics.F90 version 2:
 if(associated(this%fscav_l)  ) deallocate(this%fscav_l  )
 if(associated(this%fswtr_l)  ) deallocate(this%fswtr_l  )
 if(associated(this%phyfctd_l)) deallocate(this%phyfctd_l)

!write(0,*) '--- end subroutine mpas_tbd_ddt_deallocate'

 end subroutine mpas_tbd_ddt_deallocate

!=================================================================================================================
 subroutine mpas_tbd_set(this,mesh,nggps_atm_inputs,nggps_tbd_ddt,statein)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in) :: mesh
 type(mpas_pool_type),intent(in) :: nggps_atm_inputs
 type(mpas_pool_type),intent(in) :: nggps_tbd_ddt
 type(state_fields_in),intent(in):: statein

!inout arguments:
 class(mpas_tbd_ddt),intent(inout):: this

!local pointers:
 integer,pointer:: numP2d,numP3d
 integer,pointer:: nCellsSolve,nVertLevels,nSoilLevels,nOznLevels,nOznPlCoeff
 real(kind=RKIND),dimension(:),pointer    :: poz
 real(kind=RKIND),dimension(:),pointer    :: acv,acvb,acvt,tprcp,srflag,tref,z_c,c_0,c_d,w_0,w_d
 real(kind=RKIND),dimension(:,:),pointer  :: slc,smc,stc,upd_mf,dwn_mf,det_mf,phy_f2d
 real(kind=RKIND),dimension(:,:,:),pointer:: ozmix
 real(kind=RKIND),dimension(:,:,:),pointer:: phy_f3d,phy_fctd

 real(kind=kind_phys),dimension(:),pointer    :: pgr,dpshc_l
 real(kind=kind_phys),dimension(:,:),pointer  :: slc_l,smc_l,stc_l,updmf_l,dwnmf_l, detmf_l,rann_l,phyf2d_l
 real(kind=kind_phys),dimension(:,:,:),pointer:: phyf3d_l,ozmix_l

 real(kind=kind_phys),dimension(:),pointer    :: gscond_psp,gscond_psp1
 real(kind=kind_phys),dimension(:,:),pointer  :: gscond_tp,gscond_tp1,gscond_qp,gscond_qp1

!... new pointers added in nuopc_physics.F90 version 2:
 real(kind=kind_phys),dimension(:),pointer    :: fscav_l,fswtr_l
 real(kind=kind_phys),dimension(:,:),pointer  :: phyfctd_l

!local variables:
 logical:: test

 integer:: iCell,k,n
 integer:: ncld,ntrac,nctp,nrcm

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_tbd_set:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)
 call mpas_pool_get_dimension(mesh,'numP2d'     ,numP2d     )
 call mpas_pool_get_dimension(mesh,'numP3d'     ,numP3d     )

 call mpas_pool_get_array(nggps_tbd_ddt,'acv'     ,acv     )
 call mpas_pool_get_array(nggps_tbd_ddt,'acvt'    ,acvt    )
 call mpas_pool_get_array(nggps_tbd_ddt,'acvb'    ,acvb    )
 call mpas_pool_get_array(nggps_tbd_ddt,'tprcp'   ,tprcp   )
 call mpas_pool_get_array(nggps_tbd_ddt,'srflag'  ,srflag  )
 call mpas_pool_get_array(nggps_tbd_ddt,'tref'    ,tref    )
 call mpas_pool_get_array(nggps_tbd_ddt,'z_c'     ,z_c     )
 call mpas_pool_get_array(nggps_tbd_ddt,'c_0'     ,c_0     )
 call mpas_pool_get_array(nggps_tbd_ddt,'c_d'     ,c_d     )
 call mpas_pool_get_array(nggps_tbd_ddt,'w_0'     ,w_0     )
 call mpas_pool_get_array(nggps_tbd_ddt,'w_d'     ,w_d     )
 call mpas_pool_get_array(nggps_tbd_ddt,'slc'     ,slc     )
 call mpas_pool_get_array(nggps_tbd_ddt,'smc'     ,smc     )
 call mpas_pool_get_array(nggps_tbd_ddt,'stc'     ,stc     )
 call mpas_pool_get_array(nggps_tbd_ddt,'upd_mf'  ,upd_mf  )
 call mpas_pool_get_array(nggps_tbd_ddt,'dwn_mf'  ,dwn_mf  )
 call mpas_pool_get_array(nggps_tbd_ddt,'det_mf'  ,det_mf  )

 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_psp' ,gscond_psp )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_psp1',gscond_psp1)
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_tp'  ,gscond_tp  )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_tp1' ,gscond_tp1 )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_qp'  ,gscond_qp  )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_qp1' ,gscond_qp1 )

 call mpas_pool_get_array(nggps_tbd_ddt,'phy_f2d' ,phy_f2d )
 call mpas_pool_get_array(nggps_tbd_ddt,'phy_f3d' ,phy_f3d )

 call mpas_pool_get_array(nggps_atm_inputs,'pl_pres',poz  )
 call mpas_pool_get_array(nggps_atm_inputs,'ozmix'  ,ozmix)

 ncld  = mdl_parameters%ncld
 ntrac = mdl_parameters%ntrac
 nctp  = mdl_parameters%nctp
 nrcm  = mdl_parameters%nrcm

 dpshc_l   => this%dpshc_l
 slc_l     => this%slc_l
 smc_l     => this%smc_l
 stc_l     => this%stc_l
 updmf_l   => this%updmf_l
 dwnmf_l   => this%dwnmf_l
 detmf_l   => this%detmf_l
 rann_l    => this%rann_l
 phyf2d_l  => this%phyf2d_l
 phyf3d_l  => this%phyf3d_l
 ozmix_l   => this%ozmix_l
!... new pointers added in nuopc_physics.F90 version 2:
 fscav_l   => this%fscav_l
 fswtr_l   => this%fswtr_l
 phyfctd_l => this%phyfctd_l

 pgr       => statein%pgr

 do iCell = 1, nCellsSolve
    dpshc_l(iCell) = 0.3_kind_phys * pgr(iCell)
    do k = 1, nSoilLevels
       slc_l(iCell,k) = slc(k,iCell)
       smc_l(iCell,k) = smc(k,iCell)
       stc_l(iCell,k) = stc(k,iCell)
    enddo
    do k = 1, nVertLevels
       updmf_l(iCell,k) = upd_mf(k,iCell) 
       dwnmf_l(iCell,k) = dwn_mf(k,iCell)
       detmf_l(iCell,k) = det_mf(k,iCell)
    enddo
 enddo
 do iCell = 1, nCellsSolve
    phyf2d_l(iCell,1) = gscond_psp(iCell)
    phyf2d_l(iCell,2) = gscond_psp1(iCell)
    phyf2d_l(iCell,3) = 0._kind_phys
 enddo
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       phyf3d_l(iCell,k,1) = gscond_tp(k,iCell)
       phyf3d_l(iCell,k,2) = gscond_qp(k,iCell)
       phyf3d_l(iCell,k,3) = gscond_tp1(k,iCell)
       phyf3d_l(iCell,k,4) = gscond_qp1(k,iCell)
    enddo
 enddo
 do iCell = 1, nCellsSolve
    do n = 1, nrcm
       rann_l(iCell,n) = 0.6_RKIND
    enddo
 enddo
 do iCell = 1, nCellsSolve
    do k = 1, nOznLevels
       do n = 1, nOznPlCoeff
          ozmix_l(iCell,k,n) = ozmix(k,n,iCell)
       enddo
    enddo
 enddo

!... new pointers added in nuopc_physics.F90 version 2: the options to use fscav_l and fwstr_l (tracer scavenging)
!... phyfctd_l (Chikira-Sugiyama) are not used in the current testing of the nggps driver and not linked to the
!... MPAS registry:
 do n = 1, ntrac-ncld+2
    fscav_l(n) = 0._kind_phys
    fswtr_l(n) = 0._kind_phys
 enddo
 do iCell = 1, nCellsSolve
    do n = 1, nctp
       phyfctd_l(iCell,n) = 0._kind_phys
    enddo 
 enddo

 call this%set( &
           dpshc_l , ozmix_l , poz       , rann_l   , xkzm_m , xkzm_h  , &
           xkzm_s  , psautco , prautco   , evpco    , wminco , acv     , &
           acvb    , acvt    , slc_l     , smc_l    , stc_l  , updmf_l , &
           dwnmf_l , detmf_l , phyf3d_l  , phyf2d_l , tprcp  , srflag  , &
           tref    , z_c     , c_0       , c_d      , w_0    , w_d     , &
           fscav_l , fswtr_l , phyfctd_l) 

!write(0,*) '--- end subroutine mpas_tbd_set'

 end subroutine mpas_tbd_set

!=================================================================================================================
 subroutine mpas_tbd_write(this,mesh,nggps_tbd_ddt)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 class(mpas_tbd_ddt),intent(in) :: this

!inout arguments:
 type(mpas_pool_type),intent(inout) :: nggps_tbd_ddt

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nSoilLevels
 real(kind=RKIND),dimension(:),pointer  :: acv,acvb,acvt,tprcp,srflag,tref,z_c,c_0,c_d,w_0,w_d
 real(kind=RKIND),dimension(:,:),pointer:: slc,smc,stc,upd_mf,dwn_mf,det_mf

 real(kind=RKIND),dimension(:),pointer  :: gscond_psp,gscond_psp1
 real(kind=RKIND),dimension(:,:),pointer:: gscond_tp,gscond_tp1,gscond_qp,gscond_qp1  

 real(kind=kind_phys),dimension(:,:),pointer  :: slc_l,smc_l,stc_l,updmf_l,dwnmf_l, detmf_l
 real(kind=kind_phys),dimension(:,:),pointer  :: phyf2d_l
 real(kind=kind_phys),dimension(:,:,:),pointer:: phyf3d_l

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_tbd_write:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)

 call mpas_pool_get_array(nggps_tbd_ddt,'acv'   ,acv   )
 call mpas_pool_get_array(nggps_tbd_ddt,'acvt'  ,acvt  )
 call mpas_pool_get_array(nggps_tbd_ddt,'acvb'  ,acvb  )
 call mpas_pool_get_array(nggps_tbd_ddt,'tprcp' ,tprcp )
 call mpas_pool_get_array(nggps_tbd_ddt,'srflag',srflag)
 call mpas_pool_get_array(nggps_tbd_ddt,'tref'  ,tref  )
 call mpas_pool_get_array(nggps_tbd_ddt,'z_c'   ,z_c   )
 call mpas_pool_get_array(nggps_tbd_ddt,'c_0'   ,c_0   )
 call mpas_pool_get_array(nggps_tbd_ddt,'c_d'   ,c_d   )
 call mpas_pool_get_array(nggps_tbd_ddt,'w_0'   ,w_0   )
 call mpas_pool_get_array(nggps_tbd_ddt,'w_d'   ,w_d   )
 call mpas_pool_get_array(nggps_tbd_ddt,'slc'   ,slc   )
 call mpas_pool_get_array(nggps_tbd_ddt,'smc'   ,smc   )
 call mpas_pool_get_array(nggps_tbd_ddt,'stc'   ,stc   )
 call mpas_pool_get_array(nggps_tbd_ddt,'upd_mf',upd_mf)
 call mpas_pool_get_array(nggps_tbd_ddt,'dwn_mf',dwn_mf)
 call mpas_pool_get_array(nggps_tbd_ddt,'det_mf',det_mf)

 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_psp' ,gscond_psp )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_psp1',gscond_psp1)
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_tp'  ,gscond_tp  )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_tp1' ,gscond_tp1 )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_qp'  ,gscond_qp  )
 call mpas_pool_get_array(nggps_tbd_ddt,'gscond_qp1' ,gscond_qp1 )

!... back to Registry:
 acv     => this%acv
 acvb    => this%acvb
 acvt    => this%acvt
 tprcp   => this%tprcp
 srflag  => this%srflag
 tref    => this%tref
 z_c     => this%z_c
 c_0     => this%c_0
 c_d     => this%c_d
 w_0     => this%w_0
 w_d     => this%w_d

 slc_l   => this%slc
 smc_l   => this%smc
 stc_l   => this%stc
 updmf_l => this%upd_mf
 dwnmf_l => this%dwn_mf
 detmf_l => this%det_mf

 phyf2d_l => this%phyf2d_l
 phyf3d_l => this%phyf3d_l

 do iCell = 1, nCellsSolve
    do k = 1, nSoilLevels
       slc(k,iCell) = slc_l(iCell,k)
       smc(k,iCell) = smc_l(iCell,k)
       stc(k,iCell) = stc_l(iCell,k)
    enddo
    do k = 1, nVertLevels
       upd_mf(k,iCell) = updmf_l(iCell,k)
       dwn_mf(k,iCell) = dwnmf_l(iCell,k)
       det_mf(k,iCell) = detmf_l(iCell,k)
    enddo
 enddo

 do iCell = 1, nCellsSolve
    gscond_psp(iCell)  = phyf2d_l(iCell,1)
    gscond_psp1(iCell) = phyf2d_l(iCell,2)
 enddo
 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       gscond_tp(k,iCell)  = phyf3d_l(iCell,k,1)
       gscond_qp(k,iCell)  = phyf3d_l(iCell,k,2)
       gscond_tp1(k,iCell) = phyf3d_l(iCell,k,3)
       gscond_qp1(k,iCell) = phyf3d_l(iCell,k,4)
    enddo
 enddo

!write(0,*) '--- end subroutine mpas_tbd_write'

 end subroutine mpas_tbd_write
 
!=================================================================================================================
 end module mpas_atmphys_nggps_tbd_ddt
!=================================================================================================================
