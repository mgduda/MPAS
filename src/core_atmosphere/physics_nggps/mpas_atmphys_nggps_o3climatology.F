! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_o3climatology
 use mpas_dmpar
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_date_time
 use mpas_atmphys_nggps_constants, only: degrad
 use mpas_atmphys_nggps_vars, only: blatc,dphiozc,latsozp,levozp,kozpl,latsozc,pl_coeff,timeoz

 implicit none
 private
 public:: init_o3climatology, &
          o3climatology_from_MPAS

!modifications to sourcecode:
!----------------------------
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.
! * relaxed hard-wiring of the array ozmix in subroutine o3climatology_from_MPAS. ozmix is now interpolated to the
!   current date with a call to subroutine monthly_interp_to_date.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-22. 

 contains

 
!=================================================================================================================
 subroutine init_o3climatology(dminfo,mesh,atm_input,ntoz)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(dm_info),intent(in):: dminfo

 integer,intent(in):: ntoz

!inout arguments:
 type(mpas_pool_type),intent(inout):: atm_input

!local pointers:
 integer,pointer:: nCells
 integer,pointer:: nOznLevels,nOznPlCoeff,nMonths

 real(kind=RKIND),dimension(:),pointer:: latCell
 real(kind=PHYSKIND),dimension(:),pointer:: pl_pres,pl_time
 real(kind=PHYSKIND),dimension(:,:,:,:),pointer:: ozmixm

!local variables:
 integer,parameter:: open_ok  = 0

 integer:: iCell
 integer:: it,k,nc
 integer:: ilat,ilev,itime,istat
 integer:: jo,j1,j2

 real(kind=PHYSKIND):: dy2y1,ry2y1,ry2y,ryy1,ylat
 real(kind=4),dimension(:),allocatable:: pl_lat4,pl_pres4,pl_time4,tempin
 real(kind=PHYSKIND),dimension(:),allocatable:: pl_lat
 real(kind=PHYSKIND),dimension(:,:,:,:),allocatable:: ozplin

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine init_o3climatology: ntoz =', ntoz

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nMonths'    ,nMonths    )
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)

 call mpas_pool_get_array(mesh,'latCell',latCell)

 call mpas_pool_get_array(atm_input,'pl_pres',pl_pres)
 call mpas_pool_get_array(atm_input,'pl_time',pl_time)
 call mpas_pool_get_array(atm_input,'ozmixm' ,ozmixm )

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
   
    if(levozp .ne. nOznLevels .or. pl_coeff .ne. nOznPlCoeff .or. &
       timeoz .ne. nMonths) then
 
       write(0,*) 'nOznLevels  = ',nOznLevels
       write(0,*) 'nMonths     = ',nMonths
       write(0,*) 'nOznPlCoeff = ',nOznPlCoeff
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
 pl_pres(:) = log(100._PHYSKIND*pl_pres(:))

 write(0,*) '... latsozp  = ',latsozp
 write(0,*) '... levozp   = ',levozp
 write(0,*) '... timeoz   = ',timeoz
 write(0,*) '... pl_coeff = ',pl_coeff
 write(0,*) '... dphiozc  =', dphiozc
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
          do ilat = 1, latsozp
             ozplin(ilat,k,nc,it) = tempin(ilat)
          enddo
       enddo
    enddo
 enddo
!do it = 1,timeoz
!   do ilat = 1, latsozp
!      do k = 1, levozp
!         write(0,203) it,ilat,k,(ozplin(ilat,k,nc,it),nc=1,pl_coeff)
!      enddo
!   enddo
!   write(0,*)
!enddo
!203 format(3i9,10(1x,e15.8))

!interpolation of the input ozone data to the MPAS mesh:
 do iCell = 1, nCells

    ylat = latCell(iCell)/degrad
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

!write(0,*) '--- end subroutine init_o3climatology'

!formats:
 102 format(i3,10(1x,e15.8))

 end subroutine init_o3climatology

!=================================================================================================================
 subroutine o3climatology_from_MPAS(current_date,mesh,nggps_atm_inputs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=*),intent(in):: current_date

!inout arguments:
 type(mpas_pool_type),intent(inout):: nggps_atm_inputs

!local pointers:
 integer,pointer:: nOznLevels,nOznPlCoeff,nMonths
 integer,pointer:: nCellsSolve

 real(kind=PHYSKIND),dimension(:,:,:),pointer:: ozmix
 real(kind=PHYSKIND),dimension(:,:,:,:),pointer:: ozmixm

!local variables:
 integer:: iCell,iM,k,nc
 real(kind=PHYSKIND),dimension(:),allocatable  :: datout
 real(kind=PHYSKIND),dimension(:,:),allocatable:: datin

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine o3climatology_from_MPAS:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nOznLevels' ,nOznLevels )
 call mpas_pool_get_dimension(mesh,'nOznPlCoeff',nOznPlCoeff)
 call mpas_pool_get_dimension(mesh,'nMonths'    ,nMonths    )

 call mpas_pool_get_array(nggps_atm_inputs,'ozmixm',ozmixm)
 call mpas_pool_get_array(nggps_atm_inputs,'ozmix' ,ozmix )

 write(0,*) '--- nOznLevels  = ', nOznLevels
 write(0,*) '--- nOznPlCoeff = ', nOznPlCoeff

 if(.not.allocated(datout)) allocate(datout(1:nCellsSolve))
 if(.not.allocated(datin) ) allocate(datin(1:nMonths,1:nCellsSolve))

 do k = 1, nOznLevels
    do nc = 1, nOznPlCoeff
       datin(1:nMonths,1:nCellsSolve) = ozmixm(1:nMonths,k,nc,1:nCellsSolve)
       call monthly_interp_to_date(nCellsSolve,current_date,datin,datout)
       ozmix(k,nc,1:nCellsSolve) = datout(1:nCellsSolve)
    enddo
 enddo

 if(allocated(datin) ) deallocate(datin)
 if(allocated(datout)) deallocate(datout)

!for now, hard wire ozmix to ozmixm for August:
!do k = 1, nOznLevels
!   do nc = 1, nOznPlCoeff
!      do iCell = 1, nCellsSolve
!         ozmix(k,nc,iCell) = ozmixm(10,k,nc,iCell)
!      enddo
!   enddo
!enddo

 write(0,*) '--- end subroutine o3climatology_from_MPAS'

!format:
 201 format(3i9,10(1x,e15.8))

 end subroutine o3climatology_from_MPAS

!=================================================================================================================
 end module mpas_atmphys_nggps_o3climatology
!=================================================================================================================
