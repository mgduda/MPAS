! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_cloud_properties
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_vars
 use machine
 use nuopc_physics

 implicit none
 private
 public:: mpas_cloud_properties

 type,extends(cloud_properties):: mpas_cloud_properties

    real(kind=kind_phys),dimension(:,:),pointer:: cnvqc_l  => null()
    real(kind=kind_phys),dimension(:,:),pointer:: fcice_l  => null()
    real(kind=kind_phys),dimension(:,:),pointer:: frain_l  => null()
    real(kind=kind_phys),dimension(:,:),pointer:: rrime_l  => null()
    real(kind=kind_phys),dimension(:,:),pointer:: cldcov_l => null()
    real(kind=kind_phys),dimension(:,:),pointer:: deltaq_l => null()
    real(kind=kind_phys),dimension(:,:),pointer:: cnvw_l   => null()
    real(kind=kind_phys),dimension(:,:),pointer:: cnvc_l   => null()

    contains
       procedure:: mpas_setphys   => mpas_cld_prop_setphys
       procedure:: mpas_setrad    => mpas_cld_prop_setrad
       procedure:: mpas_write     => mpas_cld_prop_write
       procedure:: allocate_mem   => mpas_cld_prop_allocate
       procedure:: deallocate_mem => mpas_cld_prop_deallocate 

 end type

!modifications to sourcecode:
!----------------------------
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.

 contains


!=================================================================================================================
 subroutine mpas_cld_prop_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_cloud_properties),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_cld_prop_allocate:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 if(.not.associated(this%cnvqc_l) ) allocate(this%cnvqc_l(1:nCellsSolve,1:nVertLevels) )
 if(.not.associated(this%fcice_l) ) allocate(this%fcice_l(1:nCellsSolve,1:nVertLevels) )
 if(.not.associated(this%frain_l) ) allocate(this%frain_l(1:nCellsSolve,1:nVertLevels) )
 if(.not.associated(this%rrime_l) ) allocate(this%rrime_l(1:nCellsSolve,1:nVertLevels) )
 if(.not.associated(this%cldcov_l)) allocate(this%cldcov_l(1:nCellsSolve,1:nVertLevels))
 if(.not.associated(this%deltaq_l)) allocate(this%deltaq_l(1:nCellsSolve,1:nVertLevels))
 if(.not.associated(this%cnvw_l)  ) allocate(this%cnvw_l(1:nCellsSolve,1:nVertLevels)  )
 if(.not.associated(this%cnvc_l)  ) allocate(this%cnvc_l(1:nCellsSolve,1:nVertLevels)  )

!write(0,*) '--- end subroutine mpas_cld_prop_allocate'

 end subroutine mpas_cld_prop_allocate

!=================================================================================================================
 subroutine mpas_cld_prop_deallocate(this)
!=================================================================================================================

 class(mpas_cloud_properties),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_cld_prop_deallocate:'

 if(associated(this%cnvqc_l) ) deallocate(this%cnvqc_l )
 if(associated(this%fcice_l) ) deallocate(this%fcice_l )
 if(associated(this%frain_l) ) deallocate(this%frain_l )
 if(associated(this%rrime_l) ) deallocate(this%rrime_l )
 if(associated(this%cldcov_l)) deallocate(this%cldcov_l)
 if(associated(this%deltaq_l)) deallocate(this%deltaq_l)
 if(associated(this%cnvw_l)  ) deallocate(this%cnvw_l  )
 if(associated(this%cnvc_l)  ) deallocate(this%cnvc_l  )

!write(0,*) '--- end subroutine mpas_cld_prop_deallocate'

 end subroutine mpas_cld_prop_deallocate

!=================================================================================================================
 subroutine mpas_cld_prop_setphys(this,mesh,nggps_cloud_properties)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_cloud_properties

!inout arguments:
 class(mpas_cloud_properties),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

 real(kind=R8KIND),dimension(:),pointer  :: cv,cvt,cvb
 real(kind=R8KIND),dimension(:,:),pointer:: cnvqc_v

 real(kind=kind_phys),dimension(:,:),pointer:: cnvqc_l

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_cld_prop_setphys:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_cloud_properties,'cv'      ,cv     )
 call mpas_pool_get_array(nggps_cloud_properties,'cvt'     ,cvt    )
 call mpas_pool_get_array(nggps_cloud_properties,'cvb'     ,cvb    )
 call mpas_pool_get_array(nggps_cloud_properties,'cnvqc_v' ,cnvqc_v)

 cnvqc_l => this%cnvqc_l

 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       cnvqc_l(iCell,k) = cnvqc_v(k,iCell)
    enddo
 enddo

 call this%setphys(flgmin, cv , cvt , cvb , cnvqc_l , sup)

!write(0,*) '--- end subroutine mpas_cld_prop_setphys'

 end subroutine mpas_cld_prop_setphys

!=================================================================================================================
 subroutine mpas_cld_prop_setrad(this,mesh,nggps_cloud_properties)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_cloud_properties

!inout arguments:
 class(mpas_cloud_properties),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

 real(kind=R8KIND),dimension(:),pointer  :: cv,cvt,cvb
 real(kind=R8KIND),dimension(:,:),pointer:: fcice,frain,rrime,cldcov,deltaq,cnvw,cnvc

 real(kind=kind_phys),dimension(:,:),pointer:: fcice_l,frain_l,rrime_l,deltaq_l,cnvw_l,cnvc_l, &
                                               cldcov_l

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------
!write(0,*) '--- enter subroutine mpas_cld_prop_setrad:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_cloud_properties,'cv'    ,cv    )
 call mpas_pool_get_array(nggps_cloud_properties,'cvt'   ,cvt   )
 call mpas_pool_get_array(nggps_cloud_properties,'cvb'   ,cvb   )
 call mpas_pool_get_array(nggps_cloud_properties,'fcice' ,fcice )
 call mpas_pool_get_array(nggps_cloud_properties,'frain' ,frain )
 call mpas_pool_get_array(nggps_cloud_properties,'rrime' ,rrime )
 call mpas_pool_get_array(nggps_cloud_properties,'deltaq',deltaq)
 call mpas_pool_get_array(nggps_cloud_properties,'cnvw'  ,cnvw  )
 call mpas_pool_get_array(nggps_cloud_properties,'cnvc'  ,cnvc  )
 call mpas_pool_get_array(nggps_cloud_properties,'cldcov',cldcov)

 fcice_l  => this%fcice_l
 frain_l  => this%frain_l
 rrime_l  => this%rrime_l
 deltaq_l => this%deltaq_l
 cnvw_l   => this%cnvw_l
 cnvc_l   => this%cnvc_l
 cldcov_l => this%cldcov_l

 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       fcice_l(iCell,k)  = fcice(k,iCell)
       frain_l(iCell,k)  = frain(k,iCell)
       rrime_l(iCell,k)  = rrime(k,iCell)
       deltaq_l(iCell,k) = deltaq(k,iCell)
       cnvw_l(iCell,k)   = cnvw(k,iCell)
       cnvc_l(iCell,k)   = cnvc(k,iCell)
       cldcov_l(iCell,k) = cldcov(k,iCell)
    enddo
 enddo

 call this%setrad(&
           cv       , cvt , cvb    , fcice_l , frain_l , rrime_l , flgmin, cldcov_l, &
           deltaq_l , sup , cnvw_l , cnvc_l)

!write(0,*) '--- end subroutine mpas_cld_prop_setrad'

 end subroutine mpas_cld_prop_setrad

!=================================================================================================================
 subroutine mpas_cld_prop_write(this,mesh,nggps_cloud_properties)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)        :: mesh
 class(mpas_cloud_properties),intent(in):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_cloud_properties

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

 real(kind=R8KIND),dimension(:),pointer  :: cv,cvt,cvb
 real(kind=R8KIND),dimension(:,:),pointer:: fcice,frain,rrime,cldcov,deltaq,cnvw,cnvc,cnvqc_v

 real(kind=kind_phys),dimension(:,:),pointer:: fcice_l,frain_l,rrime_l,deltaq_l,cnvw_l,cnvc_l,cnvqc_l

!local variables:
 integer:: iCell,k

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_cld_prop_write:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_cloud_properties,'cv'     ,cv     )
 call mpas_pool_get_array(nggps_cloud_properties,'cvt'    ,cvt    )
 call mpas_pool_get_array(nggps_cloud_properties,'cvb'    ,cvb    )
 call mpas_pool_get_array(nggps_cloud_properties,'cnvqc_v',cnvqc_v)
 call mpas_pool_get_array(nggps_cloud_properties,'fcice'  ,fcice  )
 call mpas_pool_get_array(nggps_cloud_properties,'frain'  ,frain  )
 call mpas_pool_get_array(nggps_cloud_properties,'rrime'  ,rrime  )
 call mpas_pool_get_array(nggps_cloud_properties,'deltaq' ,deltaq )
 call mpas_pool_get_array(nggps_cloud_properties,'cnvw'   ,cnvw   )
 call mpas_pool_get_array(nggps_cloud_properties,'cnvc'   ,cnvc   )

!... back to Registry:
 cv       => this%cv
 cvt      => this%cvt
 cvb      => this%cvb

 cnvqc_l  => this%cnvqc_v
 fcice_l  => this%fcice
 frain_l  => this%frain
 rrime_l  => this%rrime
 deltaq_l => this%deltaq
 cnvw_l   => this%cnvw
 cnvc_l   => this%cnvc

 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       cnvqc_v(k,iCell) = cnvqc_l(iCell,k)
       fcice(k,iCell)   = fcice_l(iCell,k)
       frain(k,iCell)   = frain_l(iCell,k)
       rrime(k,iCell)   = rrime_l(iCell,k)
       deltaq(k,iCell)  = deltaq_l(iCell,k)
       cnvw(k,iCell)    = cnvw_l(iCell,k)
       cnvc(k,iCell)    = cnvc_l(iCell,k)
    enddo
 enddo

!write(0,*) '--- end subroutine mpas_cld_prop_write'

 end subroutine mpas_cld_prop_write

!=================================================================================================================
 end module mpas_atmphys_nggps_cloud_properties
!=================================================================================================================
