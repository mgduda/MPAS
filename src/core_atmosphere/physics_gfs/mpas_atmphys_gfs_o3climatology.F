! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_gfs_o3climatology
 use mpas_dmpar
 use mpas_kind_types
 use mpas_grid_types

 use mpas_atmphys_gfs_constants, only: degrad
 use mpas_atmphys_gfs_vars, only: blatc,dphiozc,latsozp,levozp,ntoz,kozpl,latsozc,pl_coeff,timeoz

 implicit none
 private
 public:: init_o3climatology, &
          o3climatology_from_MPAS

!>\brief
!> mpas_atmphys_gfs_o3climatology contains the subroutines needed to initialize, interpolate, and
!> update the input coefficients needed in the parameterization of prognostic ozone. 
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-01-27.


 contains
 
!==================================================================================================
 subroutine init_o3climatology(dminfo,mesh,atm_input)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(dm_info),intent(in):: dminfo

!inout arguments:
 type(atm_input_type),intent(inout):: atm_input

!local variables:
 integer,parameter:: open_ok  = 0

 integer:: iCell
 integer:: it,k,nc
 integer:: ilat,ilev,itime,istat
 integer:: jo,j1,j2

 real(kind=RKIND):: dy2y1,ry2y1,ry2y,ryy1,ylat
 real(kind=4),dimension(:),allocatable:: pl_lat4,pl_pres4,pl_time4,tempin
 real(kind=RKIND),dimension(:),allocatable:: pl_lat
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: ozplin

 real(kind=RKIND),dimension(:),pointer:: pl_pres,pl_time
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ozmixm

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine init_o3climatology: ntoz =', ntoz

 pl_pres => atm_input % pl_pres % array
 pl_time => atm_input % pl_time % array
 ozmixm  => atm_input % ozmixm  % array

 if(ntoz .le. 0) then
    write(0,*) '--- ozone parameterization option no longer supported. Reset ntoz to 2.'
    call mpas_dmpar_abort(dminfo)
 else                       
!open input file needed in the parameterization of prognostic ozone:
    write(0,*) '--- opening file for prognostic ozone:'
    open(kozpl,file='global_o3prdlos.f77',action='READ',form='UNFORMATTED', &
         status='old',iostat=istat)
    if(istat /= open_ok) then
       write(0,*) '--- subroutine init_o3climatology: failure opening global_o3prdlos.f77'
       call mpas_dmpar_abort(dminfo)
    endif
    rewind (kozpl)
    read (kozpl) pl_coeff,latsozp,levozp,timeoz
   
    if(levozp .ne. mesh%nOznLevels .or. pl_coeff .ne. mesh%nOznPlCoeff .or. &
       timeoz .ne. mesh%nMonths) then
 
       write(0,*) 'nOznLevels  = ',mesh%nOznLevels
       write(0,*) 'nMonths     = ',mesh%nMonths
       write(0,*) 'nOznPlCoeff = ',mesh%nOznPlCoeff
       write(0,*)
       write(0,*) 'latsozp  = ',latsozp
       write(0,*) 'levozp   = ',levozp
       write(0,*) 'timeoz   = ',timeoz
       write(0,*) 'pl_coeff = ',pl_coeff
       write(0,*) '--- subroutine init_o3climatology: ozone dimensions do not match Registry'
       call mpas_dmpar_abort(dminfo)
    endif
    if(.not. allocated(pl_lat)  ) allocate(pl_lat(latsozp)   )
    if(.not. allocated(pl_lat4) ) allocate(pl_lat4(latsozp)  )
    if(.not. allocated(pl_pres4)) allocate(pl_pres4(levozp)  )
    if(.not. allocated(pl_time4)) allocate(pl_time4(timeoz+1))

    rewind (kozpl)
    read (kozpl) pl_coeff,latsozp,levozp,timeoz,pl_lat4,pl_pres4,pl_time4
    pl_pres(:) = pl_pres4(:)
    pl_lat(:)  = pl_lat4(:)
    pl_time(:) = pl_time4(:)
    latsozc = 2
    blatc   = 0.0
 endif
 dphiozc    = -(blatc+blatc)/(latsozc-1)
 pl_pres(:) = log(100._RKIND*pl_pres(:))

!write(0,*) 'latsozp  = ',latsozp
!write(0,*) 'levozp   = ',levozp
!write(0,*) 'timeoz   = ',timeoz
!write(0,*) 'pl_coeff = ',pl_coeff
!write(0,*) 'dphiozc  =', dphiozc
!write(0,*) '... p_lat:'
!do ilat = 1,latsozp
!   write(0,102) ilat, pl_lat(ilat)
!enddo
!write(0,*) '... pl_pres:'
!do ilev = 1,levozp
!   write(0,102) ilev, pl_pres(ilev)
!enddo
!write(0,*) '... pl_time:'
!do itime = 1,timeoz+1
!   write(0,102) itime,pl_time(itime)
!enddo

!read input ozone data:
 if(.not. allocated(tempin)) allocate(tempin(latsozp))
 if(.not. allocated(ozplin)) allocate(ozplin(latsozp,levozp,pl_coeff,timeoz))
 do it = 1, timeoz
    do nc = 1, pl_coeff
       do k = 1,levozp
          read(kozpl) tempin
!         write(0,*) tempin
          ozplin(:,k,nc,it) = tempin(:)
       enddo
    enddo
 enddo

!interpolation of the input ozone data to the MPAS mesh:
 do iCell = 1, mesh%nCells

    ylat = mesh%latCell%array(iCell)/degrad
    do jo = 1, latsozp-1
       if(ylat.ge.pl_lat(jo) .and. ylat.lt.pl_lat(jo+1)) then
          j1 = jo
          j2 = jo+1
          dy2y1 = pl_lat(j2)-pl_lat(j1)
          ry2y  = (pl_lat(j2)-ylat)/dy2y1
          ryy1  = (ylat-pl_lat(j1))/dy2y1
       elseif(ylat .lt. pl_lat(1)) then
          j1 = 1
          j2 = 1
          ry2y = 1.0
          ryy1 = 0.0
       elseif(ylat .ge. pl_lat(latsozp)) then
          j1 = latsozp
          j2 = latsozp
          ry2y = 0.0
          ryy1 = 1.0
       endif
    enddo

    do it = 1, timeoz
       do k = 1, levozp
          do nc = 1, pl_coeff
             ozmixm(it,k,nc,iCell) = ry2y*ozplin(j1,k,nc,it) + ryy1*ozplin(j2,k,nc,it)
          enddo
       enddo
    enddo
 enddo

!deallocation of local arrays:
 if(allocated(pl_lat)  ) deallocate(pl_lat  )
 if(allocated(pl_lat4) ) deallocate(pl_lat4 )
 if(allocated(pl_pres4)) deallocate(pl_pres4)
 if(allocated(pl_time4)) deallocate(pl_time4)
 if(allocated(tempin)  ) deallocate(tempin  )
 if(allocated(ozplin)  ) deallocate(ozplin  )

 write(0,*) '--- end subroutine init_o3climatology'

!formats:
 102 format(i3,10(1x,e15.8))

 end subroutine init_o3climatology

!==================================================================================================
 subroutine o3climatology_from_MPAS(mesh,atm_input,diag_physics)
!==================================================================================================

!input arguments:
 type(mesh_type),intent(in):: mesh
 type(atm_input_type),intent(in):: atm_input

!inout arguments:
 type(diag_physics_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell,k,nc
 integer:: nOznLevels,nOznPlCoeff,nMonths,nCellsSolve

 real(kind=RKIND),dimension(:,:,:),pointer:: o3clim
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ozmixm

!--------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine o3climatology_from_MPAS:'

 nOznLevels  = mesh % nOznLevels
 nOznPlCoeff = mesh % nOznPlCoeff
 nCellsSolve = mesh % nCellsSolve
 nMonths     = mesh % nMonths

 ozmixm => atm_input    % ozmixm % array
 o3clim => diag_physics % o3clim % array

!for now, hard wire o3clim to ozmixm for October:
 do k = 1, nOznLevels
    do nc = 1, nOznPlCoeff
       do iCell = 1, nCellsSolve
          o3clim(k,nc,iCell) = ozmixm(10,k,nc,iCell)
       enddo
    enddo
 enddo

 do iCell = 1, 10
    do k = 1, nOZnLevels
       write(0,201) iCell,k,(o3clim(k,nc,iCell),nc=1,nOznPlCoeff)
    enddo
 enddo

 write(0,*) '--- end subroutine o3climatology_from_MPAS'

!format:
 201 format(2i9,10(1x,e15.8))

 end subroutine o3climatology_from_MPAS

!==================================================================================================
 end module mpas_atmphys_gfs_o3climatology
!==================================================================================================
