! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_radiation_tendencies
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use machine
 use nuopc_physics

 implicit none
 private
 public:: mpas_radiation_tendencies

 type,extends(radiation_tendencies):: mpas_radiation_tendencies
    real(kind=kind_phys),dimension(:,:),pointer  :: htrsw_l => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: htrlw_l => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: dtdtr_l => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: swhc_l  => null()
    real(kind=kind_phys),dimension(:,:),pointer  :: hlwc_l  => null()
    real(kind=kind_phys),dimension(:,:,:),pointer:: hlwd_l  => null()

    contains
       procedure:: mpas_setrad    => mpas_rad_tend_set
       procedure:: mpas_write     => mpas_rad_tend_write
       procedure:: allocate_mem   => mpas_rad_tend_set_allocate
       procedure:: deallocate_mem => mpas_rad_tend_set_deallocate 

 end type

!modifications to sourcecode:
!----------------------------
! * in subroutine mpas_rad_tend_set, initialized the variables rqtk and dtdtr_l to zero.
!   Laura D. Fowler (laura@ucar.edu) / 2015-11-15. 
! * commented out some print statements.
!   laura D. Fowler (laura@ucar.edu) / 2016-01.14.

 contains


!=================================================================================================================
 subroutine mpas_rad_tend_set_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_radiation_tendencies),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve
 integer,pointer:: nVertLevels

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_rad_tend_set_allocate:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 if(.not.associated(this%htrsw_l) ) allocate(this%htrsw_l(1:nCellsSolve,1:nVertLevels) )
 if(.not.associated(this%htrlw_l) ) allocate(this%htrlw_l(1:nCellsSolve,1:nVertLevels) )
 if(.not.associated(this%dtdtr_l) ) allocate(this%dtdtr_l(1:nCellsSolve,1:nVertLevels) )
 if(.not.associated(this%swhc_l)  ) allocate(this%swhc_l(1:nCellsSolve,1:nVertLevels)  )
 if(.not.associated(this%hlwc_l)  ) allocate(this%hlwc_l(1:nCellsSolve,1:nVertLevels)  )
 if(.not.associated(this%hlwd_l)  ) allocate(this%hlwd_l(1:nCellsSolve,1:nVertLevels,6))

!write(0,*) '--- end subroutine mpas_rad_tend_set_allocate'

 end subroutine mpas_rad_tend_set_allocate

!=================================================================================================================
 subroutine mpas_rad_tend_set_deallocate(this)
!=================================================================================================================

 class(mpas_radiation_tendencies),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_rad_tend_set_deallocate:'

 if(associated(this%htrsw_l)) deallocate(this%htrsw_l)
 if(associated(this%htrlw_l)) deallocate(this%htrlw_l)
 if(associated(this%dtdtr_l)) deallocate(this%dtdtr_l)
 if(associated(this%swhc_l) ) deallocate(this%swhc_l )
 if(associated(this%hlwc_l) ) deallocate(this%hlwc_l )
 if(associated(this%hlwd_l) ) deallocate(this%hlwd_l )

!write(0,*) '--- end subroutine mpas_rad_tend_set_deallocate'
!write(0,*)

 end subroutine mpas_rad_tend_set_deallocate

!=================================================================================================================
 subroutine mpas_rad_tend_set(this,mesh,nggps_radiation_tendencies)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: nggps_radiation_tendencies

!inout arguments:
 class(mpas_radiation_tendencies),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 real(kind=RKIND),dimension(:),pointer    :: sfalb,coszen,tsflw,semis,coszdg,rqtk
 real(kind=RKIND),dimension(:,:),pointer  :: htrsw,htrlw,dtdtr,swhc,hlwc
 real(kind=RKIND),dimension(:,:,:),pointer:: hlwd

 real(kind=kind_phys),dimension(:,:),pointer  :: htrsw_l,htrlw_l,dtdtr_l,swhc_l,hlwc_l
 real(kind=kind_phys),dimension(:,:,:),pointer:: hlwd_l

!local variables:
 integer:: iCell,k,n

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_rad_tend_set:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_radiation_tendencies,'sfalb' ,sfalb )
 call mpas_pool_get_array(nggps_radiation_tendencies,'coszen',coszen)
 call mpas_pool_get_array(nggps_radiation_tendencies,'tsflw' ,tsflw )
 call mpas_pool_get_array(nggps_radiation_tendencies,'semis' ,semis )
 call mpas_pool_get_array(nggps_radiation_tendencies,'coszdg',coszdg)
 call mpas_pool_get_array(nggps_radiation_tendencies,'rqtk'  ,rqtk  )
 call mpas_pool_get_array(nggps_radiation_tendencies,'htrsw' ,htrsw )
 call mpas_pool_get_array(nggps_radiation_tendencies,'htrlw' ,htrlw )
 call mpas_pool_get_array(nggps_radiation_tendencies,'dtdtr' ,dtdtr )
 call mpas_pool_get_array(nggps_radiation_tendencies,'swhc'  ,swhc  )
 call mpas_pool_get_array(nggps_radiation_tendencies,'hlwc'  ,hlwc  )
 call mpas_pool_get_array(nggps_radiation_tendencies,'hlwd'  ,hlwd  )

 htrsw_l => this%htrsw_l
 htrlw_l => this%htrlw_l
 dtdtr_l => this%dtdtr_l
 swhc_l  => this%swhc_l
 hlwc_l  => this%hlwc_l
 hlwd_l  => this%hlwd_l

 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       htrsw_l(iCell,k) = htrsw(k,iCell)
       htrlw_l(iCell,k) = htrlw(k,iCell)
!      dtdtr_l(iCell,k) = dtdtr(k,iCell)
       swhc_l(iCell,k)  = swhc(k,iCell)
       hlwc_l(iCell,k)  = hlwc(k,iCell)
       do n = 1, 6
          hlwd_l(iCell,k,n) = hlwd(n,k,iCell)
       enddo   
    enddo
 enddo

!In gbphys, the variables dtdtr and rqtk are accumulated, but I am not sure at what point and where they are
!initialized. To avoid the occurrence of an initialized variables, I reset them to zero prior to the set call:
 do iCell = 1, nCellsSolve
    rqtk(iCell) = 0._kind_phys
    do k = 1, nVertLevels
       dtdtr_l(iCell,k) = 0._kind_phys
    enddo
 enddo

 call this%set(htrsw_l,sfalb,coszen,htrlw_l,tsflw,semis,coszdg,rqtk,hlwd_l,dtdtr_l,swhc_l,hlwc_l)

!write(0,*) '--- end subroutine mpas_rad_tend_set'

 end subroutine mpas_rad_tend_set

!=================================================================================================================
 subroutine mpas_rad_tend_write(this,mesh,nggps_radiation_tendencies)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in)            :: mesh
 class(mpas_radiation_tendencies),intent(in):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_radiation_tendencies

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 real(kind=RKIND),dimension(:),pointer    :: rqtk
 real(kind=RKIND),dimension(:,:),pointer  :: dtdtr,swhc,hlwc
 real(kind=RKIND),dimension(:,:,:),pointer:: hlwd

 real(kind=kind_phys),dimension(:,:),pointer  :: dtdtr_l,swhc_l,hlwc_l
 real(kind=kind_phys),dimension(:,:,:),pointer:: hlwd_l

!local variables:
 integer:: iCell,k,n

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_rad_tend_write:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(nggps_radiation_tendencies,'rqtk'  ,rqtk  )
 call mpas_pool_get_array(nggps_radiation_tendencies,'dtdtr' ,dtdtr )
 call mpas_pool_get_array(nggps_radiation_tendencies,'swhc'  ,swhc  )
 call mpas_pool_get_array(nggps_radiation_tendencies,'hlwc'  ,hlwc  )
 call mpas_pool_get_array(nggps_radiation_tendencies,'hlwd'  ,hlwd  )

!... back to Registry:
 rqtk   => this%rqtk

 dtdtr_l => this%dtdtr
 swhc_l  => this%swhc
 hlwc_l  => this%hlwc
 hlwd_l  => this%hlwd

 do iCell = 1, nCellsSolve
    do k = 1, nVertLevels
       dtdtr(k,iCell) = dtdtr_l(iCell,k) 
       swhc(k,iCell)  = swhc_l(iCell,k) 
       hlwc(k,iCell)  = hlwc_l(iCell,k) 
       do n = 1, 6
          hlwd(n,k,iCell) = hlwd_l(iCell,k,n)
       enddo   
    enddo
 enddo

!write(0,*) '--- end subroutine mpas_rad_tend_write'

 end subroutine mpas_rad_tend_write

!=================================================================================================================
 end module mpas_atmphys_nggps_radiation_tendencies
!=================================================================================================================
