! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_init_atm_static_gfs
 use mpas_configure
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_read_met

 implicit none
 private
 public:: init_atm_static_fields_forGFS

!>\brief MPAS initialization of the static fields using GFS binary data. This initialization is
!> needed for running the GFS physics suite.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2013-09-01.
!>


 contains

!==================================================================================================
 subroutine init_atm_static_fields_forGFS(mesh)
!==================================================================================================

!inout arguments:
 type(mesh_type),intent(inout):: mesh

!local variables:
 type(proj_info):: proj
 type(dm_info),pointer :: dminfo

 character(len=StrKIND):: fname
 character(len=StrKIND):: dir_gwdo

 integer:: i,iz,j
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

!--------------------------------------------------------------------------------------------------

 dminfo => mesh % block % domain % dminfo

 write(0,*)
 write(0,*) '--- enter subroutine init_atm_static_fields_forGFS:'

!
! Interpolate HGT
!
 nx = 3072
 ny = 1536
 nz = 1
 isigned     = 1
 endian      = 0
 wordsize    = 2
 scalefactor = 1._RKIND
 dx          = 0.117_RKIND
 dy          = 0.117_RKIND
 known_x     = 1._RKIND
 known_y     = 1._RKIND
 known_lat   = 89.91_RKIND
 known_lon   = 0.05859375
 mesh%ter%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'orography/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max orography = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min orography = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_GAUSS, proj,   &
              nlat   = ny/2,      &
              latinc = dy,        &
              loninc = dx,        &
              knowni = known_x,   &
              knownj = known_y,   &
              lat1   = known_lat, &
              lon1   = known_lon)

 interp_list(1) = FOUR_POINT
 interp_list(2) = SEARCH
 interp_list(3) = 0

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell) * DEG_PER_RAD
    lon = mesh%lonCell%array(iCell) * DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    mesh % ter % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                                0.0_RKIND,interp_list,1)
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate TER'
 write(0,*)
 go to 100


!
! Interpolate IVGTYP (LU_INDEX):
!
 nx = 360
 ny = 180
 nz = 1
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1._RKIND
 dy            = 1._RKIND
 known_x       = 1._RKIND
 known_y       = 1._RKIND
 known_lat     =  -89.5_RKIND
 known_lon     = -179.5_RKIND
 mesh%lu_index%array(:) = 0

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'vegtype/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max vegtype = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min vegtype = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(1) = SEARCH
 interp_list(2) = 0

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell) * DEG_PER_RAD
    lon = mesh%lonCell%array(iCell) * DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)
    if(lon > 180._RKIND) then
       lon = lon - 360._RKIND
       call latlon_to_ij(proj,lat,lon,x,y)
    endif
    rPoint = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                             -1.e30_RKIND,interp_list,1)
    mesh%lu_index%array(iCell) = nint(rPoint)
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate lu_index'
 write(0,*)


!
! Interpolate ISLTYP (SOILCAT_TOP):
!
 nx = 360
 ny = 180
 nz = 1
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1._RKIND
 dy            = 1._RKIND
 known_x       = 1._RKIND
 known_y       = 1._RKIND
 known_lat     =  -89.5_RKIND
 known_lon     = -179.5_RKIND
 mesh%soilcat_top%array(:) = 0

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'soiltype/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max soiltype = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min soiltype = ',minval(xarray(1:nx,1:ny))

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

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell) * DEG_PER_RAD
    lon = mesh%lonCell%array(iCell) * DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)
    if(lon > 180._RKIND) then
       lon = lon - 360._RKIND
       call latlon_to_ij(proj,lat,lon,x,y)
    endif
    rPoint = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                             -1.e30_RKIND,interp_list,1)
    mesh%soilcat_top%array(iCell) = nint(rPoint)
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate soilcat_top'
 write(0,*)


! CORRECT INCONSISTENT SOIL AND LAND USE DATA AS WITH WPS DATA:
 do iCell = 1,mesh%nCells
    if(mesh%lu_index%array(iCell) == 0 .or. &
       mesh%soilcat_top%array(iCell) == 0) then
       if(mesh%lu_index%array(iCell) /= 0) then
          write(0,*) 'Turning lu_index into water at ', iCell
          mesh%lu_index%array(iCell) = 0
       endif
       if(mesh%soilcat_top%array(iCell) /= 0) then
          write(0,*) 'Turning soilcat_top into water at ', iCell
          mesh%soilcat_top%array(iCell) = 0
       endif
    endif
 enddo


!
! Derive LANDMASK
!
 mesh%landmask%array(:) = 0
 do iCell = 1, mesh%nCells
    if (mesh%lu_index%array(iCell) /= 0) mesh%landmask%array(iCell) = 1
 end do
 write(0,*) '--- end interpolate LANDMASK'
 write(0,*)


!
! Interpolate SOILTEMP:
!
 nx = 138
 ny = 116
 nz = 1
 isigned       = 0
 endian        = 0
 wordsize      = 2
 scalefactor   = 0.1_RKIND
 dx            = 2.60870_RKIND
 dy            = 1.55172_RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.22414_RKIND
 known_lon     =   1.30436_RKIND
 mesh%soiltemp%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'tg3/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max soiltemp  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min soiltemp  = ',minval(xarray(1:nx,1:ny))

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

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell) * DEG_PER_RAD
    lon = mesh%lonCell%array(iCell) * DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    if(y < dble(1) ) y = dble(1 )
    if(y > dble(ny)) y = dble(ny) 
    mesh%soiltemp%array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                                 -1.e30_RKIND,interp_list,1)
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
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1._RKIND
 dy            = 1._RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.5_RKIND
 known_lon     = -179.5_RKIND
 mesh%snoalb%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'snoalb/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max snoalb  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min snoalb  = ',minval(xarray(1:nx,1:ny))

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

 do iCell = 1,mesh%nCells
    if(mesh%landmask%array(iCell) == 1) then
    lat = mesh%latCell%array(iCell)*DEG_PER_RAD
    lon = mesh%lonCell%array(iCell)*DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
!   if(x < 1._RKIND) then
!      write(0,101) iCell,lat,y,lon,x
!      lon = dble(nx) + lon
!      call latlon_to_ij(proj, lat, lon, x, y)
!   elseif(x >= dble(nx)) then
    if(x >= dble(nx)) then
       lon = lon - dble(nx)
       call latlon_to_ij(proj,lat,lon,x,y)
    end if
    if(y < dble(1) ) y = dble(1 )
    if(y > dble(ny)) y = dble(ny) 
    mesh%snoalb%array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                               -1.e30_RKIND,interp_list,1)
    if(mesh%snoalb%array(iCell) ==  -1.e30_RKIND) then                                        
!      write(0,101) iCell,lat,y,lon,x,mesh%snoalb%array(iCell)
       mesh%snoalb%array(iCell) = 0._RKIND
    endif
    endif
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate SNOALB'
 write(0,*)
!stop

!
! Interpolate GREENFRAC:
!
 nx = 2500
 ny = 1250
 nz = 12
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 0.01_RKIND
 dx            = 0.144_RKIND
 dy            = 0.144_RKIND
 known_x       = 1._RKIND
 known_y       = 1._RKIND
 known_lat     =  -89.928_RKIND
 known_lon     = -179.928_RKIND
 mesh%greenfrac%array(:,:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'vegfra/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

!allocate(xarray(nx,ny))
!allocate(rarray(nx,ny,nz))
!call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
!                  scalefactor,wordsize,istatus)
!call init_atm_check_read_error(istatus,fname,dminfo)
 
!do iz = 1, nz
!   xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
!   write(0,*) iz,'... max shdmax = ',maxval(xarray(1:nx,1:ny))
!   write(0,*) iz,'... min shdmax = ',minval(xarray(1:nx,1:ny))
!enddo

!call map_set(PROJ_LATLON, proj,    &
!             latinc = dy,          &
!             loninc = dx,          &
!             knowni = known_x,     &
!             knownj = known_y,     &
!             lat1   = known_lat,   &
!             lon1   = known_lon)

!interp_list(1) = FOUR_POINT
!interp_list(2) = W_AVERAGE4
!interp_list(3) = W_AVERAGE16
!interp_list(4) = SEARCH
!interp_list(5) = 0

!do iCell = 1,mesh%nCells
!   lat = mesh % latCell % array(iCell) * DEG_PER_RAD
!   lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
!   call latlon_to_ij(proj, lat, lon, x, y)
!   if(x < 0.5) then
!      lon = lon + 360.0
!      call latlon_to_ij(proj, lat, lon, x, y)
!   else if(x >= 2500.5) then
!      lon = lon - 360.0
!      call latlon_to_ij(proj, lat, lon, x, y)
!   end if
!   if (y < 1.0) y = 1.0
!   if (y > 1249.0) y = 1249.0
!   mesh % shdmax % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
!                                         -1.e30_RKIND,interp_list,1)
!enddo
!deallocate(rarray)
!deallocate(xarray)
!write(0,*) '--- end interpolate SHDMAX'
!write(0,*)


!
! Interpolate SHDMAX:
!
 nx = 2500
 ny = 1250
 nz = 1
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 0.144_RKIND
 dy            = 0.144_RKIND
 known_x       = 1._RKIND
 known_y       = 1._RKIND
 known_lat     =  -89.93_RKIND
 known_lon     = -180.00_RKIND
 mesh%shdmax%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'shdmax/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max shdmax = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min shdmax = ',minval(xarray(1:nx,1:ny))

!call map_set(PROJ_LATLON, proj,    &
!             latinc = dy,          &
!             loninc = dx,          &
!             knowni = known_x,     &
!             knownj = known_y,     &
!             lat1   = known_lat,   &
!             lon1   = known_lon)

!interp_list(1) = FOUR_POINT
!interp_list(2) = W_AVERAGE4
!interp_list(3) = W_AVERAGE16
!interp_list(4) = SEARCH
!interp_list(5) = 0

!do iCell = 1,mesh%nCells
!   lat = mesh % latCell % array(iCell) * DEG_PER_RAD
!   lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
!   call latlon_to_ij(proj, lat, lon, x, y)
!   if(x < 0.5) then
!      lon = lon + 360.0
!      call latlon_to_ij(proj, lat, lon, x, y)
!   else if(x >= 2500.5) then
!      lon = lon - 360.0
!      call latlon_to_ij(proj, lat, lon, x, y)
!   end if
!   if (y < 1.0) y = 1.0
!   if (y > 1249.0) y = 1249.0
!   mesh % shdmax % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
!                                         -1.e30_RKIND,interp_list,1)
!enddo
!deallocate(rarray)
!deallocate(xarray)
!write(0,*) '--- end interpolate SHDMAX'
!write(0,*)

!
! Interpolate SHDMIN:
!
!nx = 2500
!ny = 1250
!nz = 1
!isigned       = 0
!endian        = 0
!wordsize      = 1
!scalefactor   = 1._RKIND
!dx            = 0.144_RKIND
!dy            = 0.144_RKIND
!known_x       = 1._RKIND
!known_y       = 1._RKIND
!known_lat     =  -89.93_RKIND
!known_lon     = -180.00_RKIND
!mesh%shdmin%array(:) = 0._RKIND

!write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
!            'shdmin/',1,'-',nx,'.',1,'-',ny
!write(0,*) trim(fname)

!allocate(xarray(nx,ny))
!allocate(rarray(nx,ny,nz))
!call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
!                  scalefactor,wordsize,istatus)
!call init_atm_check_read_error(istatus,fname,dminfo)
!xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
!write(0,*) '... max shdmin = ',maxval(xarray(1:nx,1:ny))
!write(0,*) '... min shdmin = ',minval(xarray(1:nx,1:ny))

!call map_set(PROJ_LATLON, proj,    &
!             latinc = dy,          &
!             loninc = dx,          &
!             knowni = known_x,     &
!             knownj = known_y,     &
!             lat1   = known_lat,   &
!             lon1   = known_lon)

!interp_list(1) = FOUR_POINT
!interp_list(2) = W_AVERAGE4
!interp_list(3) = W_AVERAGE16
!interp_list(4) = SEARCH
!interp_list(5) = 0

!do iCell = 1,mesh%nCells
!   lat = mesh % latCell % array(iCell) * DEG_PER_RAD
!   lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
!   call latlon_to_ij(proj, lat, lon, x, y)
!   if(x < 0.5) then
!      lon = lon + 360.0
!      call latlon_to_ij(proj, lat, lon, x, y)
!   else if(x >= 2500.5) then
!      lon = lon - 360.0
!      call latlon_to_ij(proj, lat, lon, x, y)
!   end if
!   if (y < 1.0) y = 1.0
!   if (y > 1249.0) y = 1249.0
!   mesh % shdmin % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
!                                         -1.e30_RKIND,interp_list,1)
!enddo
!deallocate(rarray)
!deallocate(xarray)
!write(0,*) '--- end interpolate SHDMIN'
!write(0,*)

!
! Interpolate FACSF:
!
 nx = 360
 ny = 180
 nz = 1
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1.0_RKIND
 dy            = 1.0_RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.5_RKIND
 known_lon     =   0.5_RKIND
 mesh%facsf%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'facsf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max facsf  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min facsf  = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(2) = FOUR_POINT
 interp_list(3) = AVERAGE4
 interp_list(1) = SEARCH
 interp_list(4) = 0
 interp_list(5) = 0 

 do iCell = 1,mesh%nCells
    lat = mesh % latCell % array(iCell) * DEG_PER_RAD
    lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
    call latlon_to_ij(proj, lat, lon, x, y)
    y = 180._RKIND-y
    if(y < 1._RKIND)   y = 1._RKIND
    if(y > 180._RKIND) y = 180._RKIND
    mesh % facsf % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          -1.e30_RKIND,interp_list,1)
!write(0,101) iCell,lat,lon,y,x,mesh%facsf%array(iCell)
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
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1.0_RKIND
 dy            = 1.0_RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.5_RKIND
 known_lon     =   0.5_RKIND
 mesh%facwf%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'facwf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 xarray(1:nx,1:ny) = rarray(1:nx,1:ny,1)
 write(0,*) '... max facwf  = ',maxval(xarray(1:nx,1:ny))
 write(0,*) '... min facwf  = ',minval(xarray(1:nx,1:ny))

 call map_set(PROJ_LATLON, proj,    &
              latinc = dy,          &
              loninc = dx,          &
              knowni = known_x,     &
              knownj = known_y,     &
              lat1   = known_lat,   &
              lon1   = known_lon)

 interp_list(2) = FOUR_POINT
 interp_list(3) = AVERAGE4
 interp_list(1) = SEARCH
 interp_list(4) = 0
 interp_list(5) = 0 

 do iCell = 1,mesh%nCells
    lat = mesh % latCell % array(iCell) * DEG_PER_RAD
    lon = mesh % lonCell % array(iCell) * DEG_PER_RAD
    call latlon_to_ij(proj, lat, lon, x, y)
    y = 180._RKIND-y
    if(y < 1._RKIND)   y = 1._RKIND
    if(y > 180._RKIND) y = 180._RKIND
    mesh % facwf % array(iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                          -1.e30_RKIND,interp_list,1)
!write(0,101) iCell,lat,lon,y,x,mesh%facwf%array(iCell)
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate FACWF'


 100 continue
!
! Interpolate ALNSF:
!
 nx = 360
 ny = 180
 nz = 4
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1.0_RKIND
 dy            = 1.0_RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.5_RKIND
 known_lon     =   0.5_RKIND
 mesh%alnsf4%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'alnsf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 
 do i = 1,nz 
    write(0,*) '... max alnsf4  = ',maxval(rarray(1:nx,1:ny,i))
    write(0,*) '... min alnsf4  = ',minval(rarray(1:nx,1:ny,i))
 enddo

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

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell)*DEG_PER_RAD
    lon = mesh%lonCell%array(iCell)*DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    y = dble(ny)-y
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)

    do iz = 1,nz
       xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
       mesh%alnsf4%array(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                                    -1.e30_RKIND,interp_list,1)
!      write(0,101) iCell,lat,lon,y,x,mesh%alnsf4%array(iz,iCell)
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALNSF4'
 write(0,*)


!
! Interpolate ALNWF:
!
 nx = 360
 ny = 180
 nz = 4
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1.0_RKIND
 dy            = 1.0_RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.5_RKIND
 known_lon     =   0.5_RKIND
 mesh%alnwf4%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'alnwf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 
 do i = 1,nz 
    write(0,*) '... max alnwf4  = ',maxval(rarray(1:nx,1:ny,i))
    write(0,*) '... min alnwf4  = ',minval(rarray(1:nx,1:ny,i))
 enddo

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

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell)*DEG_PER_RAD
    lon = mesh%lonCell%array(iCell)*DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    y = dble(ny)-y
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)

    do iz = 1,nz
       xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
       mesh%alnwf4%array(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                                    -1.e30_RKIND,interp_list,1)
!      write(0,101) iCell,lat,lon,y,x,mesh%alnwf4%array(iz,iCell)
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALNWF4'
 write(0,*)

 
 !
! Interpolate ALVSF:
!
 nx = 360
 ny = 180
 nz = 4
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1.0_RKIND
 dy            = 1.0_RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.5_RKIND
 known_lon     =   0.5_RKIND
 mesh%alvsf4%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'alvsf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 
 do i = 1,nz 
    write(0,*) '... max alvsf4  = ',maxval(rarray(1:nx,1:ny,i))
    write(0,*) '... min alvsf4  = ',minval(rarray(1:nx,1:ny,i))
 enddo

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

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell)*DEG_PER_RAD
    lon = mesh%lonCell%array(iCell)*DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    y = dble(ny)-y
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)

    do iz = 1,nz
       xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
       mesh%alvsf4%array(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                                    -1.e30_RKIND,interp_list,1)
!      write(0,101) iCell,lat,lon,y,x,mesh%alvsf4%array(iz,iCell)
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALVSF4'
 write(0,*)


!
! Interpolate ALVWF:
!
 nx = 360
 ny = 180
 nz = 4
 isigned       = 0
 endian        = 0
 wordsize      = 1
 scalefactor   = 1._RKIND
 dx            = 1.0_RKIND
 dy            = 1.0_RKIND
 known_x       = 1.0_RKIND
 known_y       = 1.0_RKIND
 known_lat     = -89.5_RKIND
 known_lon     =   0.5_RKIND
 mesh%alvwf4%array(:) = 0._RKIND

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(config_geog_data_path)// &
             'alvwf/',1,'-',nx,'.',1,'-',ny
 write(0,*) trim(fname)

 allocate(xarray(nx,ny))
 allocate(rarray(nx,ny,nz))
 call read_geogrid(fname,len_trim(fname),rarray,nx,ny,nz,isigned,endian, &
                   scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname,dminfo)
 
 do i = 1,nz 
    write(0,*) '... max alvwf4  = ',maxval(rarray(1:nx,1:ny,i))
    write(0,*) '... min alvwf4  = ',minval(rarray(1:nx,1:ny,i))
 enddo

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

 do iCell = 1,mesh%nCells
    lat = mesh%latCell%array(iCell)*DEG_PER_RAD
    lon = mesh%lonCell%array(iCell)*DEG_PER_RAD
    call latlon_to_ij(proj,lat,lon,x,y)
    y = dble(ny)-y
    if(y < dble(1) ) y = dble(1)
    if(y > dble(ny)) y = dble(ny)

    do iz = 1,nz
       xarray(1:nx,1:ny) = rarray(1:nx,1:ny,iz)
       mesh%alvwf4%array(iz,iCell) = interp_sequence(x,y,1,xarray,1,nx,1,ny,1,1, &
                                                    -1.e30_RKIND,interp_list,1)
!      write(0,101) iCell,lat,lon,y,x,mesh%alvwf4%array(iz,iCell)
    enddo
 enddo
 deallocate(rarray)
 deallocate(xarray)
 write(0,*) '--- end interpolate ALNWF4'
 write(0,*)

 101 format(i9,8(1x,e15.8))
 102 format(i9,4(1x,e15.8),i9)

 write(0,*) '--- end subroutine init_atm_static_fields_forGFS'

 end subroutine init_atm_static_fields_forGFS

!==================================================================================================
 end module mpas_init_atm_static_gfs
!==================================================================================================
