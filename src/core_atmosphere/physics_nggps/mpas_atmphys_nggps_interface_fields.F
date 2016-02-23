! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_interface_fields
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use machine
 use nuopc_physics
 use module_radlw_parameters, only: sfcflw_type
 use module_radsw_parameters, only: sfcfsw_type

 implicit none
 private
 public:: mpas_interface_fields

 type,extends(interface_fields):: mpas_interface_fields
       type(sfcflw_type),dimension(:),pointer:: sfcflw_l => null()
       type(sfcfsw_type),dimension(:),pointer:: sfcfsw_l => null()

       real(kind=kind_phys),dimension(:),pointer    :: sfcnsw_l => null()
       real(kind=kind_phys),dimension(:,:),pointer  :: htrlw0_l => null()
       real(kind=kind_phys),dimension(:,:),pointer  :: htrsw0_l => null()
       real(kind=kind_phys),dimension(:,:,:),pointer:: htrlwb_l => null()
       real(kind=kind_phys),dimension(:,:,:),pointer:: htrswb_l => null()

    contains
       procedure:: mpas_setrad    => mpas_interface_fld_setrad
       procedure:: mpas_setphys   => mpas_interface_fld_setphys
       procedure:: mpas_write     => mpas_interface_fld_write
       procedure:: allocate_mem   => mpas_interface_fld_allocate
       procedure:: deallocate_mem => mpas_interface_fld_deallocate 

 end type

!modifications to sourcecode:
!----------------------------
! * in subroutine mpas_interface_fld_setphys, removed the variables sfcdsw, sfcnsw, and sfcdlw, which are not
!   needed and not defined in Registry.physics_nggps.interface_fields.xml.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-19.
! * in mpas_interface_fld_setrad, removed arrays htrlwb and htrswb declared in Registry.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-09.
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.

 contains


!=================================================================================================================
 subroutine mpas_interface_fld_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_interface_fields),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nBands

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_interface_fld_allocate:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nBands'     ,nBands     )

 if(.not.associated(this%sfcnsw_l) ) allocate(this%sfcnsw_l(1:nCellsSolve))
 if(.not.associated(this%sfcflw_l) ) allocate(this%sfcflw_l(1:nCellsSolve))
 if(.not.associated(this%sfcfsw_l) ) allocate(this%sfcfsw_l(1:nCellsSolve))
 if(.not.associated(this%htrlw0_l) ) allocate(this%htrlw0_l(1:nCellsSolve,1:nVertLevels))
 if(.not.associated(this%htrsw0_l) ) allocate(this%htrsw0_l(1:nCellsSolve,1:nVertLevels))
 if(.not.associated(this%htrlwb_l) ) allocate(this%htrlwb_l(1:nCellsSolve,1:nVertLevels,1:nBands))
 if(.not.associated(this%htrswb_l) ) allocate(this%htrswb_l(1:nCellsSolve,1:nVertLevels,1:nBands))

!write(0,*) '--- end subroutine mpas_interface_fld_allocate'

 end subroutine mpas_interface_fld_allocate

!=================================================================================================================
 subroutine mpas_interface_fld_deallocate(this)
!=================================================================================================================

 class(mpas_interface_fields),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_interface_fld_deallocate'

 if(associated(this%sfcnsw_l) ) deallocate(this%sfcnsw_l)
 if(associated(this%sfcflw_l) ) deallocate(this%sfcflw_l)
 if(associated(this%sfcfsw_l) ) deallocate(this%sfcfsw_l)
 if(associated(this%htrlw0_l) ) deallocate(this%htrlw0_l)
 if(associated(this%htrsw0_l) ) deallocate(this%htrsw0_l)
 if(associated(this%htrlwb_l) ) deallocate(this%htrlwb_l)
 if(associated(this%htrswb_l) ) deallocate(this%htrswb_l)

!write(0,*) '--- end subroutine mpas_interface_fld_deallocate'

 end subroutine mpas_interface_fld_deallocate

!=================================================================================================================
 subroutine mpas_interface_fld_setphys(this,mesh,nggps_interface_fields)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_interface_fields

!inout arguments:
 class(mpas_interface_fields),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve

 real(kind=R8KIND),dimension(:),pointer:: sfclwdn,sfcswdn,sfcswup
 real(kind=R8KIND),dimension(:),pointer:: sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu, &
                                         sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd
 real(kind=R8KIND),dimension(:),pointer:: dusfc_cpl,dvsfc_cpl,dtsfc_cpl,dqsfc_cpl,dlwsfc_cpl,dswsfc_cpl,        &
                                         dnirbm_cpl,dnirdf_cpl,dvisbm_cpl,dvisdf_cpl,rain_cpl,nlwsfc_cpl,      &
                                         nswsfc_cpl,nnirbm_cpl,nnirdf_cpl,nvisbm_cpl,nvisdf_cpl,dusfci_cpl,    &
                                         dvsfci_cpl,dtsfci_cpl,dqsfci_cpl,dlwsfci_cpl,dswsfci_cpl,dnirbmi_cpl, & 
                                         dnirdfi_cpl,dvisbmi_cpl,dvisdfi_cpl,nlwsfci_cpl,nswsfci_cpl,          &
                                         nnirbmi_cpl,nnirdfi_cpl,nvisbmi_cpl,nvisdfi_cpl,t2mi_cpl,q2mi_cpl,    &
                                         u10mi_cpl,v10mi_cpl,tseai_cpl,psurfi_cpl,oro_cpl,slmsk_cpl
 real(kind=kind_phys),dimension(:),pointer:: xt,xs,xu,xv,xz,zm,xtts,xzts,d_conv,ifd,dt_cool,qrain

 real(kind=kind_phys),dimension(:),pointer:: sfcnsw_l

!local variables:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_interface_fld_setphys:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(nggps_interface_fields,'sfclwdn'    ,sfclwdn  )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswdn'    ,sfcswdn  )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswup'    ,sfcswup  )
 call mpas_pool_get_array(nggps_interface_fields,'sfcnirbmu'  ,sfcnirbmu)
 call mpas_pool_get_array(nggps_interface_fields,'sfcnirdfu'  ,sfcnirdfu)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisbmu'  ,sfcvisbmu)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisdfu'  ,sfcvisdfu)
 call mpas_pool_get_array(nggps_interface_fields,'sfcnirbmd'  ,sfcnirbmd)
 call mpas_pool_get_array(nggps_interface_fields,'sfcnirdfd'  ,sfcnirdfd)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisbmd'  ,sfcvisbmd)
 call mpas_pool_get_array(nggps_interface_fields,'sfcvisdfd'  ,sfcvisdfd)

 call mpas_pool_get_array(nggps_interface_fields,'dusfc_cpl'  ,dusfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvsfc_cpl'  ,dvsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dtsfc_cpl'  ,dtsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dqsfc_cpl'  ,dqsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dlwsfc_cpl' ,dlwsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dswsfc_cpl' ,dswsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirbm_cpl' ,dnirbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirdf_cpl' ,dnirdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisbm_cpl' ,dvisbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisdf_cpl' ,dvisdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'rain_cpl'   ,rain_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nlwsfc_cpl' ,nlwsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nswsfc_cpl' ,nswsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirbm_cpl' ,nnirbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirdf_cpl' ,nnirdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisbm_cpl' ,nvisbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisdf_cpl' ,nvisdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'xt'         ,xt)
 call mpas_pool_get_array(nggps_interface_fields,'xs'         ,xs)
 call mpas_pool_get_array(nggps_interface_fields,'xu'         ,xu)
 call mpas_pool_get_array(nggps_interface_fields,'xv'         ,xv)
 call mpas_pool_get_array(nggps_interface_fields,'xz'         ,xz)
 call mpas_pool_get_array(nggps_interface_fields,'zm'         ,zm)
 call mpas_pool_get_array(nggps_interface_fields,'xtts'       ,xtts)
 call mpas_pool_get_array(nggps_interface_fields,'xzts'       ,xzts)
 call mpas_pool_get_array(nggps_interface_fields,'d_conv'     ,d_conv)
 call mpas_pool_get_array(nggps_interface_fields,'ifd'        ,ifd)
 call mpas_pool_get_array(nggps_interface_fields,'dt_cool'    ,dt_cool)
 call mpas_pool_get_array(nggps_interface_fields,'qrain'      ,qrain)
 call mpas_pool_get_array(nggps_interface_fields,'dusfci_cpl' ,dusfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvsfci_cpl' ,dvsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dtsfci_cpl' ,dtsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dqsfci_cpl' ,dqsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dlwsfci_cpl',dlwsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dswsfci_cpl',dswsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirbmi_cpl',dnirbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirdfi_cpl',dnirdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisbmi_cpl',dvisbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisdfi_cpl',dvisdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nlwsfci_cpl',nlwsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nswsfci_cpl',nswsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirbmi_cpl',nnirbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirdfi_cpl',nnirdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisbmi_cpl',nvisbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisdfi_cpl',nvisdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'t2mi_cpl'   ,t2mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'q2mi_cpl'   ,q2mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'u10mi_cpl'  ,u10mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'v10mi_cpl'  ,v10mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'tseai_cpl'  ,tseai_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'psurfi_cpl' ,psurfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'oro_cpl'    ,oro_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'slmsk_cpl'  ,slmsk_cpl)

 sfcnsw_l => this%sfcnsw_l

 do iCell = 1, nCellsSolve
    sfcnsw_l(iCell) = sfcswdn(iCell) - sfcswup(iCell)
 enddo

 call this%setphys( &
           sfcswdn     , sfcnsw_l    , sfclwdn     , sfcnirbmu   , sfcnirdfu   , sfcvisbmu   , &
           sfcvisdfu   , sfcnirbmd   , sfcnirdfd   , sfcvisbmd   , sfcvisdfd   , dusfc_cpl   , &
           dvsfc_cpl   , dtsfc_cpl   , dqsfc_cpl   , dlwsfc_cpl  , dswsfc_cpl  , dnirbm_cpl  , &
           dnirdf_cpl  , dvisbm_cpl  , dvisdf_cpl  , rain_cpl    , nlwsfc_cpl  , nswsfc_cpl  , &
           nnirbm_cpl  , nnirdf_cpl  , nvisbm_cpl  , nvisdf_cpl  , xt          , xs          , &
           xu          , xv          , xz          , zm          , xtts        , xzts        , &
           d_conv      , ifd         , dt_cool     , qrain       , dusfci_cpl  , dvsfci_cpl  , &
           dtsfci_cpl  , dqsfci_cpl  , dlwsfci_cpl , dswsfci_cpl , dnirbmi_cpl , dnirdfi_cpl , &
           dvisbmi_cpl , dvisdfi_cpl , nlwsfci_cpl , nswsfci_cpl , nnirbmi_cpl , nnirdfi_cpl , &
           nvisbmi_cpl , nvisdfi_cpl , t2mi_cpl    , q2mi_cpl    , u10mi_cpl   , v10mi_cpl   , &
           tseai_cpl   , psurfi_cpl  , oro_cpl     , slmsk_cpl)

!write(0,*) '--- end subroutine mpas_interface_fld_setphys'

 end subroutine mpas_interface_fld_setphys

!=================================================================================================================
 subroutine mpas_interface_fld_setrad(this,mesh,nggps_interface_fields)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_interface_fields

!inout arguments:
 class(mpas_interface_fields),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nBands
 real(kind=R8KIND),dimension(:),pointer    :: sfclwup,sfclwdn,sfclwupc,sfclwdnc
 real(kind=R8KIND),dimension(:),pointer    :: sfcswup,sfcswdn,sfcswupc,sfcswdnc
 real(kind=R8KIND),dimension(:,:),pointer  :: htrlw0,htrsw0

!local nggps pointers:
 type(sfcflw_type),dimension(:),pointer:: sfcflw_l
 type(sfcfsw_type),dimension(:),pointer:: sfcfsw_l

 real(kind=kind_phys),dimension(:,:),pointer  :: htrlw0_l,htrsw0_l
 real(kind=kind_phys),dimension(:,:,:),pointer:: htrlwb_l,htrswb_l

!local variables:
 integer:: iCell,k,n

!-----------------------------------------------------------------------------------------------------------------
!write(0,*) '--- enter subroutine mpas_interface_fld_setrad'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nBands'     ,nBands     )

 call mpas_pool_get_array(nggps_interface_fields,'sfclwup' ,sfclwup )
 call mpas_pool_get_array(nggps_interface_fields,'sfclwdn' ,sfclwdn )
 call mpas_pool_get_array(nggps_interface_fields,'sfclwupc',sfclwupc)
 call mpas_pool_get_array(nggps_interface_fields,'sfclwdnc',sfclwdnc)
 call mpas_pool_get_array(nggps_interface_fields,'sfcswup' ,sfcswup )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswdn' ,sfcswdn )
 call mpas_pool_get_array(nggps_interface_fields,'sfcswupc',sfcswupc)
 call mpas_pool_get_array(nggps_interface_fields,'sfcswdnc',sfcswdnc)
 call mpas_pool_get_array(nggps_interface_fields,'htrlw0'  ,htrlw0  )
 call mpas_pool_get_array(nggps_interface_fields,'htrsw0'  ,htrsw0  )

 sfcflw_l => this%sfcflw_l
 sfcfsw_l => this%sfcfsw_l
 htrlw0_l => this%htrlw0_l
 htrsw0_l => this%htrsw0_l
 htrlwb_l => this%htrlwb_l
 htrswb_l => this%htrswb_l

 do iCell = 1, nCellsSolve
    sfcflw_l(iCell)%upfxc = sfclwup(iCell)
    sfcflw_l(iCell)%dnfxc = sfclwdn(iCell)
    sfcflw_l(iCell)%upfx0 = sfclwupc(iCell)
    sfcflw_l(iCell)%dnfx0 = sfclwdnc(iCell)
    sfcfsw_l(iCell)%upfxc = sfcswup(iCell)
    sfcfsw_l(iCell)%dnfxc = sfcswdn(iCell)
    sfcfsw_l(iCell)%upfx0 = sfcswupc(iCell)
    sfcfsw_l(iCell)%dnfx0 = sfcswdnc(iCell)
    do k = 1, nVertLevels
       htrlw0_l(iCell,k) = htrlw0(k,iCell)
       htrsw0_l(iCell,k) = htrsw0(k,iCell)
       do n = 1, nBands
          htrlwb_l(iCell,k,n) = 0._kind_phys
          htrswb_l(iCell,k,n) = 0._kind_phys
       enddo
    enddo
 enddo

 call this%setrad(sfcfsw_l,sfcflw_l,htrlw0_l,htrsw0_l,htrswb_l,htrlwb_l)

!write(0,*) '--- end subroutine mpas_interface_fld_setrad'

 end subroutine mpas_interface_fld_setrad

!=================================================================================================================
 subroutine mpas_interface_fld_write(this,mesh,nggps_interface_fields)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)        :: mesh
 class(mpas_interface_fields),intent(in):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_interface_fields

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels,nBands
 real(kind=R8KIND),dimension(:),pointer:: sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu, &
                                         sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd
 real(kind=R8KIND),dimension(:),pointer:: dusfc_cpl,dvsfc_cpl,dtsfc_cpl,dqsfc_cpl,dlwsfc_cpl,dswsfc_cpl,        &
                                         dnirbm_cpl,dnirdf_cpl,dvisbm_cpl,dvisdf_cpl,rain_cpl,nlwsfc_cpl,      &
                                         nswsfc_cpl,nnirbm_cpl,nnirdf_cpl,nvisbm_cpl,nvisdf_cpl,dusfci_cpl,    &
                                         dvsfci_cpl,dtsfci_cpl,dqsfci_cpl,dlwsfci_cpl,dswsfci_cpl,dnirbmi_cpl, &
                                         dnirdfi_cpl,dvisbmi_cpl,dvisdfi_cpl,nlwsfci_cpl,nswsfci_cpl,          &
                                         nnirbmi_cpl,nnirdfi_cpl,nvisbmi_cpl,nvisdfi_cpl,t2mi_cpl,q2mi_cpl,    &
                                         u10mi_cpl,v10mi_cpl,tseai_cpl,psurfi_cpl,oro_cpl,slmsk_cpl
 real(kind=R8KIND),dimension(:),pointer:: xt,xs,xu,xv,xz,zm,xtts,xzts,d_conv,ifd,dt_cool,qrain

 real(kind=R8KIND),dimension(:,:),pointer  :: htrlw0,htrsw0

!local nggps pointers:
 real(kind=kind_phys),dimension(:,:),pointer  :: htrlw0_l,htrsw0_l

!local variables:
 integer:: iCell,k,n

!-----------------------------------------------------------------------------------------------------------------
!write(0,*) '--- enter subroutine mpas_interface_fld_write:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nBands'     ,nBands     )

!... from physics:
 call mpas_pool_get_array(nggps_interface_fields,'dusfc_cpl'  ,dusfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvsfc_cpl'  ,dvsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvsfc_cpl'  ,dvsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dqsfc_cpl'  ,dqsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dlwsfc_cpl' ,dlwsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dswsfc_cpl' ,dswsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirbm_cpl' ,dnirbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirdf_cpl' ,dnirdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisbm_cpl' ,dvisbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisdf_cpl' ,dvisdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'rain_cpl'   ,rain_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nlwsfc_cpl' ,nlwsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nswsfc_cpl' ,nswsfc_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirbm_cpl' ,nnirbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirdf_cpl' ,nnirdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisbm_cpl' ,nvisbm_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisdf_cpl' ,nvisdf_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'xt'         ,xt)
 call mpas_pool_get_array(nggps_interface_fields,'xs'         ,xs)
 call mpas_pool_get_array(nggps_interface_fields,'xu'         ,xu)
 call mpas_pool_get_array(nggps_interface_fields,'xv'         ,xv)
 call mpas_pool_get_array(nggps_interface_fields,'xz'         ,xz)
 call mpas_pool_get_array(nggps_interface_fields,'zm'         ,zm)
 call mpas_pool_get_array(nggps_interface_fields,'xtts'       ,xtts)
 call mpas_pool_get_array(nggps_interface_fields,'xzts'       ,xzts)
 call mpas_pool_get_array(nggps_interface_fields,'d_conv'     ,d_conv)
 call mpas_pool_get_array(nggps_interface_fields,'ifd'        ,ifd)
 call mpas_pool_get_array(nggps_interface_fields,'dt_cool'    ,dt_cool)
 call mpas_pool_get_array(nggps_interface_fields,'qrain'      ,qrain)
 call mpas_pool_get_array(nggps_interface_fields,'dusfci_cpl' ,dusfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvsfci_cpl' ,dvsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dtsfci_cpl' ,dtsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dqsfci_cpl' ,dqsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dlwsfci_cpl',dlwsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dswsfci_cpl',dswsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirbmi_cpl',dnirbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dnirdfi_cpl',dnirdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisbmi_cpl',dvisbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'dvisdfi_cpl',dvisdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nlwsfci_cpl',nlwsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nswsfci_cpl',nswsfci_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirbmi_cpl',nnirbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nnirdfi_cpl',nnirdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisbmi_cpl',nvisbmi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'nvisdfi_cpl',nvisdfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'t2mi_cpl'   ,t2mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'q2mi_cpl'   ,q2mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'u10mi_cpl'  ,u10mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'v10mi_cpl'  ,v10mi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'tseai_cpl'  ,tseai_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'psurfi_cpl' ,psurfi_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'oro_cpl'    ,oro_cpl)
 call mpas_pool_get_array(nggps_interface_fields,'slmsk_cpl'  ,slmsk_cpl)

!... from radiation:
 call mpas_pool_get_array(nggps_interface_fields,'htrlw0'  ,htrlw0  )
 call mpas_pool_get_array(nggps_interface_fields,'htrsw0'  ,htrsw0  )

!... back to Registry:
 dusfc_cpl  => this%dusfc_cpl
 dvsfc_cpl  => this%dvsfc_cpl
 dvsfc_cpl  => this%dvsfc_cpl
 dqsfc_cpl  => this%dqsfc_cpl
 dlwsfc_cpl => this%dlwsfc_cpl
 dswsfc_cpl => this%dswsfc_cpl
 dnirbm_cpl => this%dnirbm_cpl
 dnirdf_cpl => this%dnirdf_cpl
 dvisbm_cpl => this%dvisbm_cpl
 dvisdf_cpl => this%dvisdf_cpl
 rain_cpl   => this%rain_cpl
 nlwsfc_cpl => this%nlwsfc_cpl
 nswsfc_cpl => this%nswsfc_cpl
 nnirbm_cpl => this%nnirbm_cpl
 nnirdf_cpl => this%nnirdf_cpl
 nvisbm_cpl => this%nvisbm_cpl
 nvisdf_cpl => this%nvisdf_cpl
 xt         => this%xt
 xs         => this%xs
 xu         => this%xu
 xv         => this%xv
 xz         => this%xz
 zm         => this%zm
 xtts       => this%xtts
 xzts       => this%xzts
 d_conv     => this%d_conv
 ifd        => this%ifd
 dt_cool    => this%dt_cool
 qrain      => this%qrain
 dusfci_cpl => this%dusfci_cpl
 dvsfci_cpl => this%dvsfci_cpl
 dtsfci_cpl => this%dtsfci_cpl
 dqsfci_cpl => this%dqsfci_cpl
 dlwsfci_cpl=> this%dlwsfci_cpl
 dswsfci_cpl=> this%dswsfci_cpl
 dnirbmi_cpl=> this%dnirbmi_cpl
 dnirdfi_cpl=> this%dnirdfi_cpl
 dvisbmi_cpl=> this%dvisbmi_cpl
 dvisdfi_cpl=> this%dvisdfi_cpl
 nlwsfci_cpl=> this%nlwsfci_cpl
 nswsfci_cpl=> this%nswsfci_cpl
 nnirbmi_cpl=> this%nnirbmi_cpl
 nnirdfi_cpl=> this%nnirdfi_cpl
 nvisbmi_cpl=> this%nvisbmi_cpl
 nvisdfi_cpl=> this%nvisdfi_cpl
 t2mi_cpl   => this%t2mi_cpl
 q2mi_cpl   => this%q2mi_cpl
 u10mi_cpl  => this%u10mi_cpl
 v10mi_cpl  => this%v10mi_cpl
 tseai_cpl  => this%tseai_cpl
 psurfi_cpl => this%psurfi_cpl
 oro_cpl    => this%oro_cpl
 slmsk_cpl  => this%slmsk_cpl

 htrlw0_l   => this%htrlw0
 htrsw0_l   => this%htrsw0

 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       htrlw0(k,iCell) = htrlw0_l(iCell,k) 
       htrsw0(k,iCell) = htrsw0_l(iCell,k)
    enddo
 enddo

!write(0,*) '--- end subroutine mpas_interface_fld_write'

 end subroutine mpas_interface_fld_write

!=================================================================================================================
 end module mpas_atmphys_nggps_interface_fields
!=================================================================================================================
