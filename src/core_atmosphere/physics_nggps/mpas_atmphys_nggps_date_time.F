! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_date_time
 use mpas_kind_types

 use mpas_atmphys_date_time, only: get_julgmt

 implicit none
 private
 public:: season_interp_to_date

!>\brief MPAS utility module for time management.
!>\author Laura D. Fowler (send comments to laura@ucar.edu)
!>\date 2013-01-27.
!>
!>\details
!> subroutines called in mpas_atmphys_gfs_date_time:
!> -------------------------------------------------
!> season_interp_to_date : interpolates seasonal-mean data to current julian day.


 contains


!=================================================================================================================
 subroutine season_interp_to_date(npoints,date_str,field_in,field_out)
!=================================================================================================================

!input arguments:
 character(len=StrKIND),intent(in):: date_str
 integer,intent(in):: npoints
 real(kind=R8KIND),intent(in) ,dimension(4,npoints):: field_in

!output arguments:
 real(kind=R8KIND),intent(out),dimension(npoints):: field_out

!local variables:
 character(len=2):: day15,mon

 integer:: l,n
 integer:: julyr,julday,int_season,season1,season2
 integer:: target_julyr,target_julday,target_date
 integer,dimension(0:5) :: season
 integer,dimension(0:13):: middle

 real(kind=R8KIND):: gmt

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine season_interp_to_date:'
!write(0,*) '--- current_date  = ',trim(date_str)

 write(day15,fmt='(I2.2)') 15
 do l = 1,12
    write(mon,fmt='(I2.2)') l
    call get_julgmt(date_str(1:4)//'-'//mon//'-'//day15//'_'//'00:00:00.0000',julyr,julday,gmt)
    middle(l) = julyr*1000 + julday
 enddo
 middle(0)  = middle(1)  - 31
 middle(13) = middle(12) + 31

 season(1) = middle(1)  !December-January-February.
 season(2) = middle(4)  !March-April-May.
 season(3) = middle(7)  !June-July-August.
 season(4) = middle(10) !September-October-November.
 season(0) = season(1) - 3*31
 season(5) = season(4) + 3*31

 call get_julgmt(date_str,target_julyr,target_julday,gmt)
 target_date = target_julyr * 1000 + target_julday

 find_season: do l = 0,4
    if(season(l).lt.target_date .and. season(l+1).ge.target_date) then
       int_season = l
       if(int_season.eq.0 .or. int_season.eq.4 ) then
          season1 = 4
          season2 = 1
       else
          season1 = int_season
          season2 = season1+1 
       endif
!      write(0,*)
!      write(0,*) '--- season1 = ',season1
!      write(0,*) '--- season2 = ',season2
       do n = 1,npoints
          field_out(n) = ( field_in(season2,n) * (target_date - season(l))    &
                       +   field_in(season1,n) * (season(l+1) - target_date)) &
                       / (season(l+1) - season(l))
       enddo
    endif
 enddo find_season

!write(0,*) '--- end season_interp_to_date'

 end subroutine season_interp_to_date

!=================================================================================================================
 end module mpas_atmphys_nggps_date_time
!=================================================================================================================


