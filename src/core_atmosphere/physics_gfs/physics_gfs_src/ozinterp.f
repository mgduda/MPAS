      SUBROUTINE setindxoz(latd,nlats,global_lats_r,
     &                     jindx1,jindx2,ddy)
!
      USE MACHINE , ONLY : kind_phys
      use resol_def
      use layout1
      use gg_def
      use ozne_def , only : jo3 => latsozp, pl_lat
!
      implicit none
!
      integer              global_lats_r(latr)
      integer latd,nlats,j,lat,i
      real(kind=kind_phys) pi
      real(kind=kind_phys) GAUL(latd),DDY(latd)
      integer JINDX1(latd),JINDX2(latd)
 
cyt      if(me.eq.0) print*,'begin setindxoz nlats=latd=',nlats,latd
 
         PI=ACOS(-1.)
         DO J=1,nlats
          lat = global_lats_r(ipt_lats_node_r-1+J)
          if (lat.le.latr2) then
            GAUL(J) = 90.0 - colrad_r(lat)*180.0/PI
          else
            GAUL(J) = -(90.0 - colrad_r(lat)*180.0/PI)
          endif
!cselaif(me.eq.0) print*,'gau(j,1) gau(j,2)',gaul(j,1),gaul(j,2)
         ENDDO
!
         DO J=1,nlats
           lat = global_lats_r(ipt_lats_node_r-1+J)
           jindx2(j) = jo3 + 1
           do i=1,jo3
             if (gaul(j) .lt. pl_lat(i)) then
               jindx2(j) = i
               exit
             endif
           enddo
           jindx1(j) = max(jindx2(j)-1,1)
           jindx2(j) = min(jindx2(j),jo3)
           if (jindx2(j) .ne. jindx1(j)) then
             DDY(j) = (gaul(j)           - pl_lat(jindx1(j)))
     &              / (pl_lat(jindx2(j)) - pl_lat(jindx1(j)))
           else
             ddy(j) = 1.0
           endif
!     print *,' j=',j,' gaul=',gaul(j),' jindx12=',jindx1(j),
!    &jindx2(j),' pl_lat=',pl_lat(jindx1(j)),pl_lat(jindx2(j))
!    &,' ddy=',ddy(j)
!csela if(me.eq.0) print*,'1st ddy(j,1) ddy(j,2),j=',ddy(j,1),ddy(j,2),j
 
         ENDDO
 
csela do j=1,nlats
csela if(me.eq.0) print*,'x1(j,1) jindx1(j,2)',jindx1(j,1),jindx1(j,2),j
csela if(me.eq.0) print*,'x2(j,1) jindx2(j,2)',jindx2(j,1),jindx2(j,2),j
csela enddo
csela do j=1,nlats
csela  if(me.eq.0) print*,'ddy(j,1) ddy(j,2)',ddy(j,1),ddy(j,2)
csela enddo
cyt   if(me.eq.0) print*,'completed setindxoz for nasa prod. and diss'
 
      RETURN
      END
!
!**********************************************************************
!
      SUBROUTINE ozinterpol(me,latd,nlats,IDATE,FHOUR,
     & jindx1,jindx2,ozplin,ozplout,ddy)
!
      USE MACHINE , ONLY : kind_phys
      use ozne_def
      implicit none
      integer             iday,j,j1,j2,l,latd,nc,n1,n2
      real(kind=kind_phys) fhour,tem, tx1, tx2
!
 
      integer  JINDX1(LATD), JINDX2(LATD)
      integer  me,idate(4),nlats
      integer  IDAT(8),JDAT(8)
!
      real(kind=kind_phys) ozplin(latsozp,levozp,pl_coeff,timeoz)
      real(kind=kind_phys) DDY(LATD)
      real(kind=kind_phys) ozplout(levozp,LATD,pl_coeff)
      real(kind=kind_phys) RINC(5), rjday
      integer jdow, jdoy, jday
!
      IDAT=0
      IDAT(1)=IDATE(4)
      IDAT(2)=IDATE(2)
      IDAT(3)=IDATE(3)
      IDAT(5)=IDATE(1)
      RINC=0.
      RINC(2)=FHOUR
      CALL W3MOVDAT(RINC,IDAT,JDAT)
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY .LT. PL_time(1)) RJDAY = RJDAY+365.
!
      n2 = timeoz + 1
      do j=1,timeoz
        if (rjday .lt. pl_time(j)) then
          n2 = j
          exit
        endif
      enddo
      n1 = n2 - 1
      if (n1 <= 0)     n1 = n1 + timeoz
      if (n2 > timeoz) n2 = n2 - timeoz

!
!     if (me .eq. 0) print *,' n1=',n1,' n2=',n2,' rjday=',rjday
!    &,'pl_time=',pl_time(n1),pl_time(n2)
!

      tx1 = (pl_time(n2) - rjday) / (pl_time(n2) - pl_time(n1))
      tx2 = 1.0 - tx1
!
      do nc=1,pl_coeff
        DO L=1,levozp
          DO J=1,nlats
            J1  = JINDX1(J)
            J2  = JINDX2(J)
            TEM = 1.0 - DDY(J)
            ozplout(L,j,nc) =
     &      tx1*(TEM*ozplin(J1,L,nc,n1)+DDY(J)*ozplin(J2,L,nc,n1))
     &    + tx2*(TEM*ozplin(J1,L,nc,n2)+DDY(J)*ozplin(J2,L,nc,n2))
          ENDDO
        ENDDO
      enddo
!
      RETURN
      END

!------------------------------------------------------
!---interpoltion on unstructured grid array 
!--Fanglin Yang, September 2012

      SUBROUTINE ozintp_pnt(me,lonr,lats_node_r,xlon,xlat,
     &                 idate,fhour,ozplin,ozplout)
!
      USE MACHINE , ONLY : kind_phys
      use ozne_def
      implicit none

!--inout variables
      integer             :: me,lonr,lats_node_r,idate(4)
      real(kind=kind_phys):: fhour
      real(kind=kind_phys):: ozplin(latsozp,levozp,pl_coeff,timeoz)
      real(kind=kind_phys):: xlon(lonr,lats_node_r)
      real(kind=kind_phys):: xlat(lonr,lats_node_r)
      real(kind=kind_phys):: ozplout(lonr,levozp,pl_coeff,lats_node_r)
!
!--local variables
      integer :: IDAT(8),JDAT(8)
      integer :: iday,i,j,jo,j1,j2,L,latd,nc,n1,n2
      integer :: jdow, jdoy, jday
      real(kind=kind_phys) RINC(5), rjday
      real(kind=kind_phys) :: ylat,dy2y1,ry2y,ryy1, tx1, tx2
      real(kind=kind_phys) :: pi
!
      pi=4.0*atan(1.0)
!
      IDAT=0
      IDAT(1)=IDATE(4)
      IDAT(2)=IDATE(2)
      IDAT(3)=IDATE(3)
      IDAT(5)=IDATE(1)
      RINC=0.
      RINC(2)=FHOUR
      CALL W3MOVDAT(RINC,IDAT,JDAT)
!
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jdat,jdow,jdoy,jday)
      rjday = jdoy + jdat(5) / 24.
      IF (RJDAY .LT. PL_time(1)) RJDAY = RJDAY+365.
!
      n2 = timeoz + 1
      do j=1,timeoz
        if (rjday .lt. pl_time(j)) then
          n2 = j
          exit
        endif
      enddo
      n1 = n2 - 1
      if (n1 <= 0)     n1 = n1 + timeoz
      if (n2 > timeoz) n2 = n2 - timeoz

!
!     if (me .eq. 0) print *,' n1=',n1,' n2=',n2,' rjday=',rjday
!    &,'pl_time=',pl_time(n1),pl_time(n2)
!

      tx1 = (pl_time(n2) - rjday) / (pl_time(n2) - pl_time(n1))
      tx2 = 1.0 - tx1
!
      do j=1,lats_node_r
      do i=1,lonr              !unstructured array
         ylat=xlat(i,j)*180/pi
       do jo=1,latsozp-1      !source data lat points
         if( ylat.ge.pl_lat(jo) .and. ylat.lt.pl_lat(jo+1)) then
           j1=jo
           j2=jo+1
           dy2y1=pl_lat(j2)-pl_lat(j1)
           ry2y=(pl_lat(j2)-ylat)/dy2y1       
           ryy1=(ylat-pl_lat(j1))/dy2y1       
         else if ( ylat.lt.pl_lat(1) ) then 
           j1=1
           j2=1
           ry2y=1.0
           ryy1=0.0
         else if ( ylat.ge.pl_lat(latsozp) ) then 
           j1=latsozp
           j2=latsozp
           ry2y=0.0
           ryy1=1.0
         endif
       enddo
       do nc=1,pl_coeff
       do L=1,levozp
        ozplout(i,L,nc,j) =
     &   tx1*(ry2y*ozplin(j1,L,nc,n1)+ryy1*ozplin(j2,L,nc,n1))
     & + tx2*(ry2y*ozplin(j1,L,nc,n2)+ryy1*ozplin(j2,L,nc,n2))
       enddo
       enddo
      enddo
      enddo
!
      RETURN
      END
