!Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_init_atm_static_nggps
 use mpas_kind_types
 use mpas_dmpar
 use mpas_pool_routines
 use mpas_init_atm_utilities
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_read_met

 implicit none
 private
 public:: init_atm_static_fields_forGFS, &
          init_atm_static_orogwd_forGFS


!>\brief MPAS initialization of the static fields using GFS binary data. This initialization is
!> needed for running the GFS physics suite.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-09-01.
!>


 contains

!=================================================================================================================
 subroutine init_atm_static_orogwd_forGFS(dminfo,mesh,dims,configs)
!=================================================================================================================

!input arguments:
 type(dm_info):: dminfo
 type(mpas_pool_type),intent(in):: dims
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: mesh

!local pointers:
 character(len=StrKIND),pointer:: config_geog_data_path

 integer,pointer:: nCells,nMtvr
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND),dimension(:,:),pointer:: orogwd

!local variables:
 character(len=StrKIND):: fname
 character(len=StrKIND+1):: geog_data_path !same as config_geog_data_path,but guaranteed to have a trailing slash

 integer,parameter:: nlon = 1760
 integer,parameter:: nlat =  880
 integer:: iCell
 integer:: i,ierr,ilat,ilon,nm

 real(kind=4),dimension(:,:,:),allocatable:: hprime

!-----------------------------------------------------------------------------------------------------------------

 write(0,*) '--- enter subroutine init_atm_static_orogwd_forGFS:'

 call mpas_pool_get_config(configs,'config_geog_data_path',config_geog_data_path)
 write(geog_data_path, '(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if (geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 end if

 call mpas_pool_get_dimension(dims,'nCells',nCells)
 call mpas_pool_get_dimension(dims,'nMtvr' ,nMtvr )

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'orogwd'  ,orogwd  )

 write(fname,'(a)') trim(geog_data_path)//'mtnvariance/global_mtnvar.t574.f77'
 write(0,*) trim(fname)
 open(unit=11,file=fname,form='unformatted',status='old',action='read',iostat=ierr)

!open(unit=11,file='/Users/LAURA/DATA/GFS_static/mntvariance/global_mtnvar.t574.f77', &
!     form='unformatted',status='old',action='read',iostat=ierr)
 if(.not.allocated(hprime)) allocate(hprime(nlon,nlat,nMtvr))
 read(11) hprime

 do iCell = 1,nCells
    ilon = min(1.0*nlon,nlon/2+0.5*nlon/nCells*(iCell-1))
    ilat = min(1.0*nlat,nlat/4+0.5*nlat/nCells*(iCell-1))
    do nm = 1,nMtvr
       orogwd(nm,iCell) = hprime(ilon,ilat,nm)
    enddo
 enddo 
 if(allocated(hprime)) deallocate(hprime)

 do nm = 1,nMtvr
    where(landmask(:) .eq. 0) orogwd(nm,:) = 0_RKIND 
    write(0,201) nm,maxval(orogwd(nm,:)),minval(orogwd(nm,:))
 enddo
 201 format(' ... nm = ',i3,'  max = ',e15.8,'  min = ',e15.8)

 write(0,*) '--- end subroutine init_atm_static_orogwd_forGFS'

 end subroutine init_atm_static_orogwd_forGFS

!=================================================================================================================
 subroutine init_atm_static_fields_forGFS(dminfo,mesh,dims,configs)
!=================================================================================================================

!input arguments:
 type(dm_info):: dminfo
 type (mpas_pool_type),intent(in):: dims
 type (mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type) ,intent(inout):: mesh

!local pointers:
 integer,pointer:: nCells
 integer,dimension(:),pointer:: landmask,lu_index,soilcat_bot,soilcat_top

 real(kind=RKIND),dimension(:),pointer  :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer  :: shdmax,shdmin,slope,snoalb,soiltemp,ter
 real(kind=RKIND),dimension(:),pointer  :: facsf,facwf
 real(kind=RKIND),dimension(:,:),pointer:: alnsf4s,alvsf4s,alnwf4s,alvwf4s,greenfrac,zorl12m

!local variables:
 type(proj_info):: proj

 character(len=StrKIND):: fname
 character(len=StrKIND):: dir_gwdo
 character(len=StrKIND+1):: geog_data_path !same as config_geog_data_path,but guaranteed to have a trailing slash
 character(len=StrKIND),pointer:: config_geog_data_path

 integer:: i,iiy,ix,iy,iz,j
 integer:: nx,ny,nz
 integer:: endian,isigned,istatus,wordsize
 integer:: iCell
 integer,dimension(5) :: interp_list

 real(kind=4):: scalefactor
 real(kind=4),dimension(:,:,:),allocatable:: rarray

 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND):: dx,dy,known_lat,known_lon,known_x,known_y
 real(kind=RKIND):: missing_value,rPoint
 real(kind=RKIND),dimension(:),allocatable:: temp
 real(kind=RKIND),dimension(:,:),allocatable:: xarray

 real(kind=RKIND):: lon2,x2,rPoint2

!additional local variables needed to interpolate the WPS_GEOG topo_30s data:
 integer:: iPoint,iTileStart,iTileEnd,jTileStart,jTileEnd
 integer,dimension(:),allocatable:: nhs

 integer,pointer:: maxEdges
 integer,dimension(:),pointer  :: nEdgesOnCell
 integer,dimension(:,:),pointer:: cellsOnCell

 real(kind=RKIND):: lat_pt,lon_pt


!-----------------------------------------------------------------------------------------------------------------

 write(0,*) '--- enter subroutine init_atm_static_fields_forGFS:'

 call mpas_pool_get_config(configs,'config_geog_data_path',config_geog_data_path)
 write(geog_data_path, '(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if (geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 end if

 call mpas_pool_get_dimension(dims,'nCells'  ,nCells  )
 call mpas_pool_get_dimension(dims,'maxEdges',maxEdges)
 
 call mpas_pool_get_array(mesh,'latCell'     ,latCell     )
 call mpas_pool_get_array(mesh,'lonCell'     ,lonCell     )
 call mpas_pool_get_array(mesh,'cellsOnCell' ,cellsOnCell )
 call mpas_pool_get_array(mesh,'nEdgesOnCell',nEdgesOnCell)

 call mpas_pool_get_array(mesh,'alnsf4s'    ,alnsf4s    )
 call mpas_pool_get_array(mesh,'alnwf4s'    ,alnwf4s    )
 call mpas_pool_get_array(mesh,'alvsf4s'    ,alvsf4s    )
 call mpas_pool_get_array(mesh,'alvwf4s'    ,alvwf4s    )
 call mpas_pool_get_array(mesh,'facsf'      ,facsf      )
 call mpas_pool_get_array(mesh,'facwf'      ,facwf      )
 call mpas_pool_get_array(mesh,'greenfrac'  ,greenfrac  )
 call mpas_pool_get_array(mesh,'landmask'   ,landmask   )
 call mpas_pool_get_array(mesh,'lu_index'   ,lu_index   )
 call mpas_pool_get_array(mesh,'shdmax'     ,shdmax     )
 call mpas_pool_get_array(mesh,'shdmin'     ,shdmin     )
 call mpas_pool_get_array(mesh,'slope'      ,slope      )
 call mpas_pool_get_array(mesh,'snoalb'     ,snoalb     )
 call mpas_pool_get_array(mesh,'soilcat_bot',soilcat_bot)
 call mpas_pool_get_array(mesh,'soilcat_top',soilcat_top)
 call mpas_pool_get_array(mesh,'soiltemp'   ,soiltemp   )
 call mpas_pool_get_array(mesh,'ter'        ,ter        )
 call mpas_pool_get_array(mesh,'zorl12m'    ,zorl12m    )

!
! Interpolate HGT
!
!nx = 3072
!ny = 1536
!nz = 1
!isigned     = 1
!endian      = 0
!wordsize    = 2
!scalefactor = 1._RKIND
!dx          = 0.117_RKIND
!dy          = 0.117_RKIND
!known_x     = 1._RKIND
!known_y     = 1._RKIND
!known_lat   = 89.91_RKIND
!known_lon   = 0.05859375
!ter(:)      = 0._RKIND

!write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
!            'orography/',1,'-',nx,'.',1,'-',ny
!write(0,*) trim(fname)

!allocate(xarray(nx,ny))
!allocate(rarray(nx,ny,nz))
!call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
!                  scalefactor,wordsize,istatus)
!call init_atm_check_read_error(istatus,fname)
!xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
!write(0,*) '... max orography = ',maxval(xarray(1:nx,1:ny))
!write(0,*) '... min orography = ',minval(xarray(1:nx,1:ny))

!call map_set(PROJ_GAUSS, proj,   &
!             nlat   = ny/2,      &
!             latinc = dy,        &
!             loninc = dx,        &
!             knowni = known_x,   &
!             knownj = known_y,   &
!             lat1   = known_lat, &
!             lon1   = known_lon)

!interp_list(1) = FOUR_POINT
!interp_list(2) = SEARCH
!interp_list(3) = 0

!do iCell = 1,nCells
!   lat = latCell(iCell)*DEG_PER_RAD
!   lon = lonCell(iCell)*DEG_PER_RAD
!   call latlon_to_ij(proj,lat,lon,x,y)
!   if(x .gt. nx) x = dble(nx)
!   ter(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,-1.e30_RKIND,interp_list,1)
!   if(ter(iCell) .eq. -1.e30_RKIND) write(0,101) iCell,lat,y,lon,x,ter(iCell) 
!enddo
!deallocate(rarray)
!deallocate(xarray)
!write(0,*) '--- end interpolate TER'
!write(0,*)


!
! Interpolate HGT using WPS_GEOG topo_30s
!
!nx = 126
!ny = 126
 nx = 1206
 ny = 1206
 nz = 1
 isigned  = 1
 endian   = 0
 wordsize = 2
 scalefactor = 1.0
 allocate(rarray(nx,ny,nz))
 allocate(nhs(nCells))
 nhs(:) = 0
 ter(:) = 0.0

 do jTileStart = 1,20401,ny-6
    jTileEnd = jTileStart + ny - 1 - 6

    do iTileStart=1,42001,nx-6
       iTileEnd = iTileStart + nx - 1 - 6
       write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'topo_30s/',iTileStart,'-',iTileEnd,'.',jTileStart,'-',jTileEnd
       write(0,*) trim(fname)

       call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                         scalefactor,wordsize,istatus)
       call init_atm_check_read_error(istatus, fname)

       iPoint = 1
       do j=4,ny-3
       do i=4,nx-3
          lat_pt = -89.99583  + (jTileStart + j - 5) * 0.0083333333
          lon_pt = -179.99583 + (iTileStart + i - 5) * 0.0083333333
          lat_pt = lat_pt * PI / 180.0
          lon_pt = lon_pt * PI / 180.0

          iPoint = nearest_cell(lat_pt,lon_pt,iPoint,nCells,maxEdges, &
                                nEdgesOnCell,cellsOnCell, &
                                latCell,lonCell)
          ter(iPoint) = ter(iPoint) + rarray(i,j,1)
          nhs(iPoint) = nhs(iPoint) + 1
       end do
       end do

    end do
 end do

 do iCell = 1,nCells
    ter(iCell) = ter(iCell) / real(nhs(iCell))
 end do
 deallocate(rarray)
 deallocate(nhs)
 write(0,*) '--- end interpolate TER'


!
! Interpolate IVGTYP (LU_INDEX):
!
 nx = 360
 ny = 180
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1._RKIND
 dx          = 1._RKIND
 dy          = 1._RKIND
 known_x     = 1._RKIND
 known_y     = 1._RKIND
 known_lat   =  -89.5_RKIND
 known_lon   = -179.5_RKIND
 lu_index(:) = 0

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'vegtype/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '--- max vegtype = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '--- min vegtype = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = SEARCH
 interp_list(2) = 0

 do iCell = 1,nCells
    lat = latCell(iCell)*DEG_PER_RAD
    lon = lonCell(iCell)*DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)
    if(lon >= 180._RKIND) then
       lon = lon - 360._RKIND
       call latlon_to_ij(proj,lat,lon,x,y)
    endif
    rPoint = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,-1.e30_RKIND,interp_list,1)
    lu_index(iCell) = nint(rPoint)
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate LU_INDEX'
 write(0,*)


!
! Interpolate ISLTYP (SOILCAT_TOP):
!
 nx = 360
 ny = 180
 nz = 1
 isigned        = 0
 endian         = 0
 wordsize       = 1
 scalefactor    = 1._RKIND
 dx             = 1._RKIND
 dy             = 1._RKIND
 known_x        = 1._RKIND
 known_y        = 1._RKIND
 known_lat      =  -89.5_RKIND
 known_lon      = -179.5_RKIND
 soilcat_top(:) = 0

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'soiltype/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '--- max soiltype = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '--- min soiltype = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = SEARCH
 interp_list(2) = 0

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = SEARCH
 interp_list(2) = 0

 do iCell = 1,nCells
    lat = latCell(iCell)*DEG_PER_RAD
    lon = lonCell(iCell)*DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)
    if(lon >= 180._RKIND) then
       lon = lon - 360._RKIND
       call latlon_to_ij(proj,lat,lon,x,y)
    endif
    rPoint = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,-1.e30_RKIND,interp_list,1)
    soilcat_top(iCell) = nint(rPoint)
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate SOILCAT_TOP'
 write(0,*)


! CORRECT INCONSISTENT SOIL AND LAND USE DATA AS WITH WPS DATA:
 do iCell = 1,nCells
    if(lu_index(iCell) == 0 .or. &
       soilcat_top(iCell) == 0) then
       if(lu_index(iCell) /= 0) then
          write(0,*) 'Turning lu_index into water at ', iCell
          lu_index(iCell) = 0
       endif
       if(soilcat_top(iCell) /= 0) then
          write(0,*) 'Turning soilcat_top into water at ', iCell
          soilcat_top(iCell) = 0
       endif
    endif
 enddo


!
! Derive LANDMASK
!
 landmask(:) = 0
 do iCell = 1, nCells
    if(lu_index(iCell) /= 0) landmask(iCell) = 1
 end do
 write(0,*) '--- end interpolate LANDMASK'
 write(0,*)

 write(0,*) '--- begin check landmask:'
 do iCell = 1, nCells
    if(landmask(iCell) == 0 .and. (lu_index(iCell)/=0 .or. soilcat_top(iCell)/=0)) &
       write(0,311) iCell,landmask(iCell),lu_index(iCell),soilcat_top(iCell)
 enddo
 311 format('1111 land 1111:',5i9)
 write(0,*) '--- end check landmask:'
 write(0,*)

!
! Interpolate SOILTEMP:
!
 nx = 138
 ny = 116
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 2
 scalefactor = 0.1_RKIND
 dx          = 2.60870_RKIND
 dy          = 1.55172_RKIND
 known_x     = 1.0_RKIND
 known_y     = 1.0_RKIND
 known_lat   = -89.22414_RKIND
 known_lon   =   1.30436_RKIND
 soiltemp(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'tg3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 do ix = 1,nx
 do iy = 1,ny
    iiy = ny-iy+1
    xarray(ix,iy) = rarray(ix,iiy,1)
 enddo
 enddo
 write(0,*) '--- max soiltemp  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '--- min soiltemp  = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells

    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj,lat,lon,x,y)
       if(y < dble(1) ) y = dble(1 )
       if(y > dble(ny)) y = dble(ny) 
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       endif
       soiltemp(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,-1.e30_RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate SOILTEMP'
 write(0,*)


!
! Interpolate SNOALB:
!
 nx = 360
 ny = 180
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1._RKIND
 dx          = 1._RKIND
 dy          = 1._RKIND
 known_x     = 1.0_RKIND
 known_y     = 1.0_RKIND
 known_lat   = -89.5_RKIND
 known_lon   = -179.5_RKIND
 snoalb(:)   = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'snoalb/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '--- max snoalb  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '--- min snoalb  = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj,lat,lon,x,y)
       if(x >= dble(nx)) then
          lon = lon - dble(nx)
          call latlon_to_ij(proj,lat,lon,x,y)
       end if
       if(y < dble(1) ) y = dble(1 )
       if(y > dble(ny)) y = dble(ny) 
       snoalb(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
       snoalb(iCell) = snoalb(iCell)*0.01_RKIND
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate SNOALB'
 write(0,*)


!
! Interpolate GREENFRAC:
!
 nx = 2500
 ny = 1250
 nz = 12
 isigned        = 0
 endian         = 0
 wordsize       = 1
 scalefactor    = 0.01_RKIND
 dx             = 0.144_RKIND
 dy             = 0.144_RKIND
 known_x        = 1._RKIND
 known_y        = 1._RKIND
 known_lat      =  -89.928_RKIND
 known_lon      = -179.928_RKIND
 greenfrac(:,:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'vegfrac12m/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(-2:nx+3,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iz = 1, nz
    xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
    xarray(0,1:ny) = rarray(nx,1:ny,iz)
    xarray(-1,1:ny) = rarray(nx-1,1:ny,iz)
    xarray(-2,1:ny) = rarray(nx-2,1:ny,iz)
    xarray(nx+1,1:ny) = rarray(1,1:ny,iz)
    xarray(nx+2,1:ny) = rarray(2,1:ny,iz)
    xarray(nx+3,1:ny) = rarray(3,1:ny,iz)
    write(0,*) iz,'--- max vegfra = ',maxval(xarray(1:nx,1:ny))
    write(0,*) iz,'--- min vegfra = ',minval(xarray(1:nx,1:ny))
    if(iz < nz) write(0,*)

    do iCell = 1,nCells
       if(landmask(iCell) == 1) then
          lat = latCell(iCell)*DEG_PER_RAD
          lon = lonCell(iCell)*DEG_PER_RAD
          call latlon_to_ij(proj, lat, lon, x, y)
          if(x < 0.5) then
             lon = lon + 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          else if(x >= 2500.5) then
             lon = lon - 360.0
             call latlon_to_ij(proj, lat, lon, x, y)
          end if
          if (y < 1.0) y = 1.0
          if (y > 1249.0) y = 1249.0
          greenfrac(iz,iCell) = interp_sequence(x,y,1,xarray,-2,nx+3,1,ny,1,1,-1.e30_RKIND,interp_list,1)
       endif
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate GREENFRAC'
 write(0,*)


!
! Interpolate SHDMAX:
!
 nx = 2500
 ny = 1250
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1._RKIND
 dx          = 0.144_RKIND
 dy          = 0.144_RKIND
 known_x     = 1._RKIND
 known_y     = 1._RKIND
 known_lat   =  -89.93_RKIND
 known_lon   = -180.00_RKIND
 shdmax(:)   = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'shdmax/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max shdmax = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min shdmax = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       else if(x >= 2500.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 1249.0) y = 1249.0
       shdmax(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,-1.e30_RKIND,interp_list,1)
       shdmax(iCell) = shdmax(iCell)*0.01_RKIND
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate SHDMAX'
 write(0,*)


!
! Interpolate SHDMIN:
!
 nx = 2500
 ny = 1250
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1._RKIND
 dx          = 0.144_RKIND
 dy          = 0.144_RKIND
 known_x     = 1._RKIND
 known_y     = 1._RKIND
 known_lat   =  -89.93_RKIND
 known_lon   = -180.00_RKIND
 shdmin(:)   = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'shdmin/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max shdmin = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min shdmin = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       else if(x >= 2500.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj, lat, lon, x, y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 1249.0) y = 1249.0
       shdmin(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,-1.e30_RKIND,interp_list,1)
       shdmin(iCell) = shdmin(iCell)*0.01_RKIND
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate SHDMIN'
 write(0,*)


!
! Interpolate FACSF:
!
 nx = 360
 ny = 180
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1._RKIND
 dx          = 1.0_RKIND
 dy          = 1.0_RKIND
 known_x     = 1.0_RKIND
 known_y     = 1.0_RKIND
 known_lat   = -89.5_RKIND
 known_lon   =   0.5_RKIND
 facsf(:)    = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'facsf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '--- max facsf  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '--- min facsf  = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = N_NEIGHBOR
 interp_list(2) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       y = 180._RKIND-y
       if(y < 1._RKIND)   y = 1._RKIND
       if(y > 180._RKIND) y = 180._RKIND
       facsf(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
       facsf(iCell) = facsf(iCell)*0.01_RKIND
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate FACSF'
 write(0,*)


!
! Interpolate FACWF:
!
 nx = 360
 ny = 180
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1._RKIND
 dx          = 1.0_RKIND
 dy          = 1.0_RKIND
 known_x     = 1.0_RKIND
 known_y     = 1.0_RKIND
 known_lat   = -89.5_RKIND
 known_lon   =   0.5_RKIND
 facwf(:)    = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'facwf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '--- max facwf  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '--- min facwf  = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(4) = N_NEIGHBOR
 interp_list(5) = 0 

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj, lat, lon, x, y)
       y = 180._RKIND-y
       if(y < 1._RKIND)   y = 1._RKIND
       if(y > 180._RKIND) y = 180._RKIND
       facwf(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
       facwf(iCell) = facwf(iCell)*0.01_RKIND
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate FACWF'
 write(0,*)


!
! Interpolate ALNSF4S:
!
 nx = 360
 ny = 180
 nz = 4
 isigned      = 0
 endian       = 0
 wordsize     = 1
 scalefactor  = 1._RKIND
 dx           = 1.0_RKIND
 dy           = 1.0_RKIND
 known_x      = 1.0_RKIND
 known_y      = 1.0_RKIND
 known_lat    = -89.5_RKIND
 known_lon    =   0.5_RKIND
 alnsf4s(:,:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'alnsf4s/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 
 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = N_NEIGHBOR
 interp_list(2) = 0

 do iz = 1,nz 
    xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
    write(0,*) '--- max alnsf4s  = ',maxval(xarray(1:nx,1:ny))
    write(0,*) '--- min alnsf4s  = ',minval(xarray(1:nx,1:ny))
    if(iz < nz) write(0,*)

    do iCell = 1,nCells
       if(landmask(iCell) == 1) then
          lat = latCell(iCell)*DEG_PER_RAD
          lon = lonCell(iCell)*DEG_PER_RAD
          call latlon_to_ij(proj,lat,lon,x,y)
          y = dble(ny)-y
          if(y < dble(1) ) y = dble(1)
          if(y > dble(ny)) y = dble(ny)
          alnsf4s(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
          alnsf4s(iz,iCell) = alnsf4s(iz,iCell)*0.01_RKIND
       endif
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALNSF4s'
 write(0,*)


!
! Interpolate ALNWF4S:
!
 nx = 360
 ny = 180
 nz = 4
 isigned      = 0
 endian       = 0
 wordsize     = 1
 scalefactor  = 1._RKIND
 dx           = 1.0_RKIND
 dy           = 1.0_RKIND
 known_x      = 1.0_RKIND
 known_y      = 1.0_RKIND
 known_lat    = -89.5_RKIND
 known_lon    =   0.5_RKIND
 alnwf4s(:,:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'alnwf4s/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 
 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = N_NEIGHBOR
 interp_list(2) = 0

 do iz = 1,nz 
    xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
    write(0,*) '--- max alnwf4s  = ',maxval(xarray(1:nx,1:ny))
    write(0,*) '--- min alnwf4s  = ',minval(xarray(1:nx,1:ny))
    if(iz < nz) write(0,*)

    do iCell = 1,nCells
       if(landmask(iCell) == 1) then
          lat = latCell(iCell)*DEG_PER_RAD
          lon = lonCell(iCell)*DEG_PER_RAD
          call latlon_to_ij(proj,lat,lon,x,y)
          y = dble(ny)-y
          if(y < dble(1) ) y = dble(1)
          if(y > dble(ny)) y = dble(ny)
          alnwf4s(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
          alnwf4s(iz,iCell) = alnwf4s(iz,iCell)*0.01_RKIND
       endif
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALNWF4S'
 write(0,*)


!
! Interpolate ALVSF4S:
!
 nx = 360
 ny = 180
 nz = 4
 isigned      = 0
 endian       = 0
 wordsize     = 1
 scalefactor  = 1._RKIND
 dx           = 1.0_RKIND
 dy           = 1.0_RKIND
 known_x      = 1.0_RKIND
 known_y      = 1.0_RKIND
 known_lat    = -89.5_RKIND
 known_lon    =   0.5_RKIND
 alvsf4s(:,:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'alvsf4s/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 
 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = N_NEIGHBOR
 interp_list(2) = 0

 do iz = 1,nz 
    xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
    write(0,*) '--- max alvsf4s  = ',maxval(xarray(1:nx,1:ny))
    write(0,*) '--- min alvsf4s  = ',minval(xarray(1:nx,1:ny))
    if(iz < nz) write(0,*)

    do iCell = 1,nCells
       if(landmask(iCell) == 1) then
          lat = latCell(iCell)*DEG_PER_RAD
          lon = lonCell(iCell)*DEG_PER_RAD
          call latlon_to_ij(proj,lat,lon,x,y)
          y = dble(ny)-y
          if(y < dble(1) ) y = dble(1)
          if(y > dble(ny)) y = dble(ny)
          alvsf4s(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
          alvsf4s(iz,iCell) = alvsf4s(iz,iCell)*0.01_RKIND
       endif
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALVSF4S'
 write(0,*)


!
! Interpolate ALVWF4S:
!
 nx = 360
 ny = 180
 nz = 4
 isigned      = 0
 endian       = 0
 wordsize     = 1
 scalefactor  = 1._RKIND
 dx           = 1.0_RKIND
 dy           = 1.0_RKIND
 known_x      = 1.0_RKIND
 known_y      = 1.0_RKIND
 known_lat    = -89.5_RKIND
 known_lon    =   0.5_RKIND
 alvwf4s(:,:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'alvwf4s/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 
 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = N_NEIGHBOR
 interp_list(2) = 0

 do iz = 1,nz 
    xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
    write(0,*) '--- max alvwf4s  = ',maxval(xarray(1:nx,1:ny))
    write(0,*) '--- min alvwf4s  = ',minval(xarray(1:nx,1:ny))
    if(iz < nz) write(0,*)

    do iCell = 1,nCells
       if(landmask(iCell) == 1) then
          lat = latCell(iCell)*DEG_PER_RAD
          lon = lonCell(iCell)*DEG_PER_RAD
          call latlon_to_ij(proj,lat,lon,x,y)
          y = dble(ny)-y
          if(y < dble(1) ) y = dble(1)
          if(y > dble(ny)) y = dble(ny)
          alvwf4s(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
          alvwf4s(iz,iCell) = alvwf4s(iz,iCell)*0.01_RKIND
       endif
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALVWF4S'
 write(0,*)


!
! Interpolate SLOPE:
!
 nx = 360
 ny = 180
 nz = 1
 isigned     = 0
 endian      = 0
 wordsize    = 1
 scalefactor = 1._RKIND
 dx          = 1.0_RKIND
 dy          = 1.0_RKIND
 known_x     = 1.0_RKIND
 known_y     = 1.0_RKIND
 known_lat   = -89.5_RKIND
 known_lon   = -179.5_RKIND
 slope(:)    = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'slope/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '--- max slope = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '--- min slope = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj,lat,lon,x,y)
       if(lon > 180._RKIND) then
          lon = lon - 360._RKIND
          call latlon_to_ij(proj,lat,lon,x,y)
       endif
       if(y < dble(1) ) y = dble(1)
       if(y > dble(ny)) y = dble(ny)
       slope(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,0._RKIND,interp_list,1)
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate SLOPE'
 write(0,*)

!
! Interpolate ZORL12M:
!
 nx = 360
 ny = 180
 nz = 12
 isigned      = 0
 endian       = 0
 wordsize     = 1
 scalefactor  = 0.1_RKIND
 dx           = 1.0_RKIND
 dy           = 1.0_RKIND
 known_x      = 1.0_RKIND
 known_y      = 1.0_RKIND
 known_lat    = -89.5_RKIND
 known_lon    = -179.5_RKIND
 zorl12m(:,:) = 10._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
             'zorl12m/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)

 do iz = 1,nz
    xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
    write(0,*) '--- max zorl12m = ',maxval(xarray(1:nx,1:ny))
    write(0,*) '--- min zorl12m = ',minval(xarray(1:nx,1:ny))
    if(iz < nz) write(0,*)

    do iCell = 1,nCells
       if(landmask(iCell) == 1) then
          lat = latCell(iCell)*DEG_PER_RAD
          lon = lonCell(iCell)*DEG_PER_RAD
          call latlon_to_ij(proj,lat,lon,x,y)
          if(lon > 180._RKIND) then
             lon = lon - 360._RKIND
             call latlon_to_ij(proj,lat,lon,x,y)
          endif
          if(y < dble(1) ) y = dble(1)
          if(y > dble(ny)) y = dble(ny)
          zorl12m(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1,10._RKIND,interp_list,1)
       endif
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ZORL12M'
 write(0,*)

 101 format(i9,8(1x,e15.8))
 102 format(i9,4(1x,e15.8),i9)
 103 format(i4,i9,8(1x,e15.8))
 104 format(2i9,8(1x,e15.8))

 write(0,*) '--- end subroutine init_atm_static_fields_forGFS'

 end subroutine init_atm_static_fields_forGFS

!=================================================================================================================
 end module mpas_init_atm_static_nggps
!=================================================================================================================
