! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================
 module mpas_atmphys_interp_diagnostics
 use mpas_kind_types
 use mpas_grid_types
 use mpas_atmphys_constants

 implicit none
 private
 public:: init_pFixed,    &
          interp_pFixed,  &
          interp_pFixed_mp


!>\brief Interpolation of physics tendencies, potential temperature, and water species to fixed
!>\pressure levels.
!>\author Laura D. Fowler (send comments to laura@ucar.edu).
!>\date 2014-09-29.
!>
!>\details
!>
!> subroutines in mpas_atmphys_interp_diagnostics:
!> -----------------------------------------------
!> interp_pFixed          : interpolate physics tendencies, except cloud microphysics tendencies.
!> interp_pFixed_mp       : interpolate cloud microphysics tendencies, and vertical soundings.
!> init_pFixed            : read fixed pressure levels. initialize diagnostics.
!> interp_tofixed_pressure: interpolation subroutine from model to fixed pressure levels (in hPa).

!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * in subroutine init_pFixed, changed the pressure level 550. hPa listed twice to 500. hPa.
!>      Laura D. Fowler (laura@ucar.edu) / 2014-10-03.


 contains


!==================================================================================================
 subroutine interp_pFixed(dt_dyn,mesh,diag,diag_physics,tend_physics)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: tend_physics

 real(kind=RKIND),intent(in):: dt_dyn

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell,k,kk

!local pointers:
 integer,pointer:: nCells,nVertLevels,nIntP

 real(kind=RKIND),dimension(:),pointer  :: pFixed
 real(kind=RKIND),dimension(:),pointer  :: raincv
 real(kind=RKIND),dimension(:,:),pointer:: pressure_p,pressure_b,exner
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac
 real(kind=RKIND),dimension(:,:),pointer:: rthcuten,rqvcuten,rqccuten,rqrcuten,rqicuten,rqscuten
 real(kind=RKIND),dimension(:,:),pointer:: rthblten,rqvblten,rqcblten,rqiblten
 real(kind=RKIND),dimension(:,:),pointer:: rthratenlw,rthratensw

 real(kind=RKIND),dimension(:,:),pointer:: cldf
 real(kind=RKIND),dimension(:,:),pointer:: dtcup,dqvcup,dqccup,dqrcup,dqicup,dqscup 
 real(kind=RKIND),dimension(:,:),pointer:: dtpbl,dqvpbl,dqcpbl,dqipbl
 real(kind=RKIND),dimension(:,:),pointer:: dtradlw,dtradsw

 real(kind=RKIND),dimension(:,:),pointer:: accldf
 real(kind=RKIND),dimension(:,:),pointer:: acdtcup,acdqvcup,acdqccup,acdqrcup,acdqicup,acdqscup 
 real(kind=RKIND),dimension(:,:),pointer:: acdtpbl,acdqvpbl,acdqcpbl,acdqipbl
 real(kind=RKIND),dimension(:,:),pointer:: acdtradlw,acdtradsw

!local interpolated fields:
 real(kind=RKIND),dimension(:,:),allocatable:: field_in,press_in
 real(kind=RKIND),dimension(:,:),allocatable:: field_interp,press_interp

!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine interp_pFixed:'

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nIntP'      ,nIntP      )

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)

 call mpas_pool_get_array(diag_physics,'pFixed'    ,pFixed    )
 call mpas_pool_get_array(diag_physics,'raincv'    ,raincv    )
 call mpas_pool_get_array(diag_physics,'cldfrac'   ,cldfrac   )
 call mpas_pool_get_array(tend_physics,'rthcuten'  ,rthcuten  )
 call mpas_pool_get_array(tend_physics,'rqvcuten'  ,rqvcuten  )
 call mpas_pool_get_array(tend_physics,'rqccuten'  ,rqccuten  )
 call mpas_pool_get_array(tend_physics,'rqrcuten'  ,rqrcuten  )
 call mpas_pool_get_array(tend_physics,'rqicuten'  ,rqicuten  )
 call mpas_pool_get_array(tend_physics,'rqscuten'  ,rqscuten  )
 call mpas_pool_get_array(tend_physics,'rthblten'  ,rthblten  )
 call mpas_pool_get_array(tend_physics,'rqvblten'  ,rqvblten  )
 call mpas_pool_get_array(tend_physics,'rqcblten'  ,rqcblten  )
 call mpas_pool_get_array(tend_physics,'rqiblten'  ,rqiblten  )
 call mpas_pool_get_array(tend_physics,'rthratenlw',rthratenlw)
 call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)

 call mpas_pool_get_array(diag_physics,'cld_interp'    ,cldf   )
 call mpas_pool_get_array(diag_physics,'dtcup_interp'  ,dtcup  )
 call mpas_pool_get_array(diag_physics,'dqvcup_interp' ,dqvcup )
 call mpas_pool_get_array(diag_physics,'dqccup_interp' ,dqccup )
 call mpas_pool_get_array(diag_physics,'dqrcup_interp' ,dqrcup )
 call mpas_pool_get_array(diag_physics,'dqicup_interp' ,dqicup )
 call mpas_pool_get_array(diag_physics,'dqscup_interp' ,dqscup )
 call mpas_pool_get_array(diag_physics,'dtpbl_interp'  ,dtpbl  )
 call mpas_pool_get_array(diag_physics,'dqvpbl_interp' ,dqvpbl )
 call mpas_pool_get_array(diag_physics,'dqcpbl_interp' ,dqcpbl )
 call mpas_pool_get_array(diag_physics,'dqipbl_interp' ,dqipbl )
 call mpas_pool_get_array(diag_physics,'dtradlw_interp',dtradlw)
 call mpas_pool_get_array(diag_physics,'dtradsw_interp',dtradsw)

 call mpas_pool_get_array(diag_physics,'accld_interp'    ,accldf   )
 call mpas_pool_get_array(diag_physics,'acdtcup_interp'  ,acdtcup  )
 call mpas_pool_get_array(diag_physics,'acdqvcup_interp' ,acdqvcup )
 call mpas_pool_get_array(diag_physics,'acdqccup_interp' ,acdqccup )
 call mpas_pool_get_array(diag_physics,'acdqrcup_interp' ,acdqrcup )
 call mpas_pool_get_array(diag_physics,'acdqicup_interp' ,acdqicup )
 call mpas_pool_get_array(diag_physics,'acdqscup_interp' ,acdqscup )
 call mpas_pool_get_array(diag_physics,'acdtpbl_interp'  ,acdtpbl  )
 call mpas_pool_get_array(diag_physics,'acdqvpbl_interp' ,acdqvpbl )
 call mpas_pool_get_array(diag_physics,'acdqcpbl_interp' ,acdqcpbl )
 call mpas_pool_get_array(diag_physics,'acdqipbl_interp' ,acdqipbl )
 call mpas_pool_get_array(diag_physics,'acdtradlw_interp',acdtradlw)
 call mpas_pool_get_array(diag_physics,'acdtradsw_interp',acdtradsw)

 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       press_in(iCell,kk) = (pressure_p(k,iCell)+pressure_b(k,iCell)) / 100._RKIND 
    enddo
 enddo

!interpolation to fixed pressure levels for fields located at cells centers and at mass points:
 if(.not.allocated(field_interp)) allocate(field_interp(nCells,nIntP))
 if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIntP))
 do k = 1, nIntP
    kk = nIntP+1-k
    do iCell = 1, nCells
       press_interp(iCell,kk) = pFixed(k)
    enddo
 enddo

!write(0,*) '--- before call to interp_to_fixed_pressure:'
!do iCell = 1, nCells
!   if(raincv(iCell) .gt. 0._RKIND) then
!      do k = 1, nVertLevels
!         write(0,201) iCell,k,(pressure_p(k,iCell)+pressure_b(k,iCell)),cldfrac(k,iCell), &
!            rthcuten(k,iCell)*exner(k,iCell),rqvcuten(k,iCell),rqccuten(k,iCell), &
!            rqrcuten(k,iCell),rqicuten(k,iCell),rqscuten(k,iCell)
!      enddo
!      write(0,*)
!   endif
!enddo

!... dtcup:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rthcuten(k,iCell)*exner(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dtcup(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqvcup:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqvcuten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqvcup(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqccup:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqccuten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqccup(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqrcup:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqrcuten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqrcup(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqicup:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqicuten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqicup(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqscup:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqscuten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqscup(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo

!... dtpbl:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rthblten(k,iCell)*exner(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dtpbl(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqvpbl:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqvblten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqvpbl(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqcpbl:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqcblten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqcpbl(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqipbl:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqiblten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqipbl(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo

!... dtradlw:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rthratenlw(k,iCell)*exner(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dtradlw(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dtradsw:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rthratensw(k,iCell)*exner(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dtradsw(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo

!... cloud fraction
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = cldfrac(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       cldf(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo

!... accumulated diagnostics:
 do k = 1, nIntP
    do iCell = 1, nCells
       accldf(k,iCell) = accldf(k,iCell) + cldf(k,iCell)
       acdtcup(k,iCell)   = acdtcup(k,iCell)   + dtcup(k,iCell)  *dt_dyn
       acdqvcup(k,iCell)  = acdqvcup(k,iCell)  + dqvcup(k,iCell) *dt_dyn
       acdqccup(k,iCell)  = acdqccup(k,iCell)  + dqccup(k,iCell) *dt_dyn
       acdqrcup(k,iCell)  = acdqrcup(k,iCell)  + dqrcup(k,iCell) *dt_dyn
       acdqicup(k,iCell)  = acdqicup(k,iCell)  + dqicup(k,iCell) *dt_dyn
       acdqscup(k,iCell)  = acdqscup(k,iCell)  + dqscup(k,iCell) *dt_dyn
       acdtpbl(k,iCell)   = acdtpbl(k,iCell)   + dtpbl(k,iCell)  *dt_dyn
       acdqvpbl(k,iCell)  = acdqvpbl(k,iCell)  + dqvpbl(k,iCell) *dt_dyn
       acdqcpbl(k,iCell)  = acdqcpbl(k,iCell)  + dqcpbl(k,iCell) *dt_dyn
       acdqipbl(k,iCell)  = acdqipbl(k,iCell)  + dqipbl(k,iCell) *dt_dyn
       acdtradlw(k,iCell) = acdtradlw(k,iCell) + dtradlw(k,iCell)*dt_dyn
       acdtradsw(k,iCell) = acdtradsw(k,iCell) + dtradsw(k,iCell)*dt_dyn
    enddo
 enddo

!write(0,*)
!write(0,*) '--- after call to interp_to_fixed_pressure:' 
!do iCell = 1, nCells
!   if(raincv(iCell) .gt. 0._RKIND) then
!      do k = 1, nIntP
!         write(0,201) iCell,k,pFixed(k),cldf(k,iCell),dtcup(k,iCell),dqvcup(k,iCell),dqccup(k,iCell), &
!                      dqrcup(k,iCell),dqicup(k,iCell),dqscup(k,iCell)
!      enddo
!      write(0,*)
!   endif
!enddo

!deallocate local arrays:
 if(allocated(press_in)) deallocate(press_in)
 if(allocated(field_in)) deallocate(field_in)

 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)

 write(0,*) '--- end subroutine interp_pFixed:'
 201 format(2i7,10(1x,e15.8))

 end subroutine interp_pFixed

!==================================================================================================
 subroutine interp_pFixed_mp(dt_dyn,configs,mesh,state,time_lev,diag,diag_physics)
!==================================================================================================

!input arguments:
 integer,intent(in):: time_lev
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: state

 real(kind=RKIND),intent(in):: dt_dyn

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables:
 integer:: iCell,k,kk

!local pointers:
 integer,pointer:: nCells,nVertLevels,nIntP
 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,dimension(:,:),pointer:: i_acthi

 real(kind=RKIND),pointer:: bucket_radt
 real(kind=RKIND),dimension(:),pointer    :: pFixed
 real(kind=RKIND),dimension(:),pointer    :: rainncv
 real(kind=RKIND),dimension(:,:),pointer  :: theta_m,pressure_p,pressure_b,exner
 real(kind=RKIND),dimension(:,:),pointer  :: rthmpten,rqvmpten,rqcmpten,rqrmpten, &
                                             rqimpten,rqsmpten,rqgmpten
 real(kind=RKIND),dimension(:,:),pointer  :: dtmp,dqvmp,dqcmp,dqrmp,dqimp,dqsmp,dqgmp
 real(kind=RKIND),dimension(:,:),pointer  :: thm,qv,qc,qr,qi,qs,qg
 real(kind=RKIND),dimension(:,:),pointer  :: thi,qvi,qci,qri,qii,qsi,qgi
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND),dimension(:,:),pointer  :: acdtmp,acdqvmp,acdqcmp,acdqrmp,acdqimp,acdqsmp,acdqgmp
 real(kind=RKIND),dimension(:,:),pointer  :: acthi,acqvi,acqci,acqri,acqii,acqsi,acqgi

!local interpolated fields:
 real(kind=RKIND),dimension(:,:),allocatable:: field_in,press_in
 real(kind=RKIND),dimension(:,:),allocatable:: field_interp,press_interp

!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine interp_pFixed_mp:'

 call mpas_pool_get_config(configs,'config_bucket_radt',bucket_radt)

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nIntP'      ,nIntP      )

 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_dimension(state,'index_qc',index_qc)
 call mpas_pool_get_dimension(state,'index_qr',index_qr)
 call mpas_pool_get_dimension(state,'index_qi',index_qi)
 call mpas_pool_get_dimension(state,'index_qs',index_qs)
 call mpas_pool_get_dimension(state,'index_qg',index_qg)

 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)
 qc => scalars(index_qc,:,:)
 qr => scalars(index_qr,:,:)
 qi => scalars(index_qi,:,:)
 qs => scalars(index_qs,:,:)
 qg => scalars(index_qg,:,:)

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)

 call mpas_pool_get_array(diag_physics,'pFixed'  ,pFixed  )
 call mpas_pool_get_array(diag_physics,'rainncv' ,rainncv )
 call mpas_pool_get_array(diag_physics,'rthmpten',rthmpten)
 call mpas_pool_get_array(diag_physics,'rqvmpten',rqvmpten)
 call mpas_pool_get_array(diag_physics,'rqcmpten',rqcmpten)
 call mpas_pool_get_array(diag_physics,'rqrmpten',rqrmpten)
 call mpas_pool_get_array(diag_physics,'rqimpten',rqimpten)
 call mpas_pool_get_array(diag_physics,'rqsmpten',rqsmpten)
 call mpas_pool_get_array(diag_physics,'rqsmpten',rqsmpten)
 call mpas_pool_get_array(diag_physics,'rqgmpten',rqgmpten)

 call mpas_pool_get_array(diag_physics,'dtmp_interp' ,dtmp )
 call mpas_pool_get_array(diag_physics,'dqvmp_interp',dqvmp)
 call mpas_pool_get_array(diag_physics,'dqcmp_interp',dqcmp)
 call mpas_pool_get_array(diag_physics,'dqrmp_interp',dqrmp)
 call mpas_pool_get_array(diag_physics,'dqimp_interp',dqimp)
 call mpas_pool_get_array(diag_physics,'dqsmp_interp',dqsmp)
 call mpas_pool_get_array(diag_physics,'dqgmp_interp',dqgmp)
 call mpas_pool_get_array(diag_physics,'th_interp'   ,thi  )
 call mpas_pool_get_array(diag_physics,'qv_interp'   ,qvi  )
 call mpas_pool_get_array(diag_physics,'qc_interp'   ,qci  )
 call mpas_pool_get_array(diag_physics,'qr_interp'   ,qri  )
 call mpas_pool_get_array(diag_physics,'qi_interp'   ,qii  )
 call mpas_pool_get_array(diag_physics,'qs_interp'   ,qsi  )
 call mpas_pool_get_array(diag_physics,'qg_interp'   ,qgi  )

 call mpas_pool_get_array(diag_physics,'acdtmp_interp' ,acdtmp )
 call mpas_pool_get_array(diag_physics,'acdqvmp_interp',acdqvmp)
 call mpas_pool_get_array(diag_physics,'acdqcmp_interp',acdqcmp)
 call mpas_pool_get_array(diag_physics,'acdqrmp_interp',acdqrmp)
 call mpas_pool_get_array(diag_physics,'acdqimp_interp',acdqimp)
 call mpas_pool_get_array(diag_physics,'acdqsmp_interp',acdqsmp)
 call mpas_pool_get_array(diag_physics,'acdqgmp_interp',acdqgmp)
 call mpas_pool_get_array(diag_physics,'acth_interp'   ,acthi  )
 call mpas_pool_get_array(diag_physics,'acqv_interp'   ,acqvi  )
 call mpas_pool_get_array(diag_physics,'acqc_interp'   ,acqci  )
 call mpas_pool_get_array(diag_physics,'acqr_interp'   ,acqri  )
 call mpas_pool_get_array(diag_physics,'acqi_interp'   ,acqii  )
 call mpas_pool_get_array(diag_physics,'acqs_interp'   ,acqsi  )
 call mpas_pool_get_array(diag_physics,'acqg_interp'   ,acqgi  )
 call mpas_pool_get_array(diag_physics,'i_acth_interp' ,i_acthi)

 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       press_in(iCell,kk) = (pressure_p(k,iCell)+pressure_b(k,iCell)) / 100._RKIND
    enddo
 enddo

!interpolation to fixed pressure levels for fields located at cells centers and
!at mass points:
 if(.not.allocated(field_interp)) allocate(field_interp(nCells,nIntP))
 if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIntP))
 do k = 1, nIntP
    kk = nIntP+1-k
    do iCell = 1, nCells
       press_interp(iCell,kk) = pFixed(k)
    enddo
 enddo

!write(0,*) '--- before call to interp_to_fixed_pressure:'
!do iCell = 1, nCells
!   if(rainncv(iCell) .gt. 0._RKIND) then
!      do k = 1, nVertLevels
!         write(0,201) iCell,k,(pressure_p(k,iCell)+pressure_b(k,iCell))/100., &
!            rthmpten(k,iCell)*exner(k,iCell),rqvmpten(k,iCell),rqcmpten(k,iCell), &
!            rqrmpten(k,iCell),rqimpten(k,iCell),rqsmpten(k,iCell)
!      enddo
!      write(0,*)
!   endif
!enddo

!... dtmp:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rthmpten(k,iCell)*exner(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dtmp(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqvmp:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqvmpten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqvmp(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqcmp:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqcmpten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqcmp(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqrmp:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqrmpten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqrmp(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqimp:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqimpten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqimp(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqsmp:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqsmpten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqsmp(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... dqgmp:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = rqgmpten(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       dqgmp(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo


!... temperature:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = theta_m(k,iCell)/(1._RKIND+R_v/R_d* max(0._RKIND,qv(k,iCell)))
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       thi(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... qv:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = qv(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       qvi(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... qc:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = qc(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       qci(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... qr:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = qr(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       qri(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... qi:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = qi(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       qii(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... qs:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = qs(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       qsi(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... qg:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = qg(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       qgi(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo

!... accumulated diagnostics:
 do k = 1, nIntP
    do iCell = 1, nCells
       acdtmp(k,iCell)  = acdtmp(k,iCell)  + dtmp(k,iCell) *dt_dyn
       acdqvmp(k,iCell) = acdqvmp(k,iCell) + dqvmp(k,iCell)*dt_dyn
       acdqcmp(k,iCell) = acdqcmp(k,iCell) + dqcmp(k,iCell)*dt_dyn
       acdqrmp(k,iCell) = acdqrmp(k,iCell) + dqrmp(k,iCell)*dt_dyn
       acdqimp(k,iCell) = acdqimp(k,iCell) + dqimp(k,iCell)*dt_dyn
       acdqsmp(k,iCell) = acdqsmp(k,iCell) + dqsmp(k,iCell)*dt_dyn

       acthi(k,iCell) = acthi(k,iCell) + thi(k,iCell)
       acqvi(k,iCell) = acqvi(k,iCell) + qvi(k,iCell)
       acqci(k,iCell) = acqci(k,iCell) + qci(k,iCell)
       acqri(k,iCell) = acqri(k,iCell) + qri(k,iCell)
       acqii(k,iCell) = acqii(k,iCell) + qii(k,iCell)
       acqsi(k,iCell) = acqsi(k,iCell) + qsi(k,iCell)
       acqgi(k,iCell) = acqgi(k,iCell) + qgi(k,iCell)

       if(acthi(k,iCell) .gt. bucket_radt) then
          i_acthi(k,iCell) = i_acthi(k,iCell) + 1
          acthi(k,iCell) = acthi(k,iCell) - bucket_radt
       endif
    enddo
 enddo

!write(0,*)
!write(0,*) '--- after call to interp_to_fixed_pressure:'
!do iCell = 1, nCells
!   if(rainncv(iCell) .gt. 0._RKIND) then
!      do k = 1, nIntP
!         write(0,201) iCell,k,pFixed(k),dtmp(k,iCell),dqvmp(k,iCell),dqcmp(k,iCell), &
!                      dqrmp(k,iCell),dqimp(k,iCell),dqsmp(k,iCell),dqgmp(k,iCell)
!      enddo
!      write(0,*)
!   endif
!enddo

!deallocate local arrays:
 if(allocated(press_in)) deallocate(press_in)
 if(allocated(field_in)) deallocate(field_in)

 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)

 write(0,*) '--- end subroutine interp_pFixed_mp:'
 201 format(2i7,10(1x,e15.8))

 end subroutine interp_pFixed_mp

!==================================================================================================
 subroutine init_pFixed(mesh,state,time_lev,diag,diag_physics)
!==================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and pointers:
 integer:: iCell,k,kk
 integer,pointer:: nCells,nIntP,nVertLevels
 integer,pointer:: index_qv
 integer,dimension(:,:),pointer:: i_acthi

 real(kind=RKIND),dimension(:),pointer  :: pFixed
 real(kind=RKIND),dimension(:,:),pointer:: pressure_p,pressure_b
 real(kind=RKIND),dimension(:,:),pointer:: accldf,acthi,acqvi,acqci,acqri,acqii,acqsi,acqgi
 real(kind=RKIND),dimension(:,:),pointer:: acdtcup,acdqvcup,acdqccup,acdqrcup,acdqicup,acdqscup
 real(kind=RKIND),dimension(:,:),pointer:: acdtpbl,acdqvpbl,acdqcpbl,acdqipbl
 real(kind=RKIND),dimension(:,:),pointer:: acdtradlw,acdtradsw
 real(kind=RKIND),dimension(:,:),pointer:: acdtmp,acdqvmp,acdqcmp,acdqrmp,acdqimp,acdqsmp,acdqgmp
 real(kind=RKIND),dimension(:,:),pointer:: theta_m,qv
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local interpolated fields:
 real(kind=RKIND),dimension(:,:),allocatable:: field_in,press_in
 real(kind=RKIND),dimension(:,:),allocatable:: field_interp,press_interp

!--------------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) '--- enter subroutine init_pFixed:'

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nIntP'      ,nIntp      )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(diag_physics,'pFixed',pFixed)
!pFixed(:) = (/1000._RKIND,950._RKIND,900._RKIND,850._RKIND,800._RKIND, &
!               750._RKIND,700._RKIND,650._RKIND,600._RKIND,550._RKIND, &
!               500._RKIND,450._RKIND,400._RKIND,350._RKIND,300._RKIND, &
!               250._RKIND,200._RKIND,150._RKIND,100._RKIND, 70._RKIND, &
!                50._RKIND, 30._RKIND, 20._RKIND, 10._RKIND/)
!... GFS levels:
 pFixed(:) = (/1000._RKIND,975._RKIND,950._RKIND,925._RKIND,900._RKIND, &
                875._RKIND,850._RKIND,825._RKIND,800._RKIND,775._RKIND, &
                750._RKIND,700._RKIND,650._RKIND,600._RKIND,550._RKIND, &
                500._RKIND,450._RKIND,400._RKIND,350._RKIND,300._RKIND, &
                250._RKIND,225._RKIND,200._RKIND,175._RKIND,150._RKIND, &
                125._RKIND,100._RKIND, 70._RKIND, 50._RKIND, 30._RKIND, &
                 20._RKIND, 10._RKIND,  7._RKIND,  5._RKIND/)

 do k = 1,nIntP
    write(0,201) k,pFixed(k)
 enddo
 201 format(i4,f6.0)

!... initialize accumulated diagnostics interpolate to the fixed pressure levels:
 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)

 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)

 call mpas_pool_get_array(diag_physics,'acdtcup_interp'  ,acdtcup  )
 call mpas_pool_get_array(diag_physics,'acdqvcup_interp' ,acdqvcup )
 call mpas_pool_get_array(diag_physics,'acdqccup_interp' ,acdqccup )
 call mpas_pool_get_array(diag_physics,'acdqrcup_interp' ,acdqrcup )
 call mpas_pool_get_array(diag_physics,'acdqicup_interp' ,acdqicup )
 call mpas_pool_get_array(diag_physics,'acdqscup_interp' ,acdqscup )
 call mpas_pool_get_array(diag_physics,'acdtmp_interp'   ,acdtmp   )
 call mpas_pool_get_array(diag_physics,'acdqvmp_interp'  ,acdqvmp  )
 call mpas_pool_get_array(diag_physics,'acdqcmp_interp'  ,acdqcmp  )
 call mpas_pool_get_array(diag_physics,'acdqrmp_interp'  ,acdqrmp  )
 call mpas_pool_get_array(diag_physics,'acdqimp_interp'  ,acdqimp  )
 call mpas_pool_get_array(diag_physics,'acdqsmp_interp'  ,acdqsmp  )
 call mpas_pool_get_array(diag_physics,'acdqgmp_interp'  ,acdqgmp  )
 call mpas_pool_get_array(diag_physics,'acdtpbl_interp'  ,acdtpbl  )
 call mpas_pool_get_array(diag_physics,'acdqvpbl_interp' ,acdqvpbl )
 call mpas_pool_get_array(diag_physics,'acdqcpbl_interp' ,acdqcpbl )
 call mpas_pool_get_array(diag_physics,'acdqipbl_interp' ,acdqipbl )
 call mpas_pool_get_array(diag_physics,'acdtradlw_interp',acdtradlw)
 call mpas_pool_get_array(diag_physics,'acdtradsw_interp',acdtradsw)

 call mpas_pool_get_array(diag_physics,'accld_interp' ,accldf)
 call mpas_pool_get_array(diag_physics,'acth_interp'  ,acthi )
 call mpas_pool_get_array(diag_physics,'acqv_interp'  ,acqvi )
 call mpas_pool_get_array(diag_physics,'acqc_interp'  ,acqci )
 call mpas_pool_get_array(diag_physics,'acqr_interp'  ,acqri )
 call mpas_pool_get_array(diag_physics,'acqi_interp'  ,acqii )
 call mpas_pool_get_array(diag_physics,'acqs_interp'  ,acqsi )
 call mpas_pool_get_array(diag_physics,'acqg_interp'  ,acqgi )
 call mpas_pool_get_array(diag_physics,'i_acth_interp',i_acthi)

 do k = 1, nIntP
    do iCell = 1, nCells
       accldf(k,iCell)    = 0._RKIND
       acdtcup(k,iCell)   = 0._RKIND
       acdqvcup(k,iCell)  = 0._RKIND
       acdqccup(k,iCell)  = 0._RKIND
       acdqrcup(k,iCell)  = 0._RKIND
       acdqicup(k,iCell)  = 0._RKIND
       acdqscup(k,iCell)  = 0._RKIND
       acdtmp(k,iCell)    = 0._RKIND
       acdqvmp(k,iCell)   = 0._RKIND
       acdqcmp(k,iCell)   = 0._RKIND
       acdqrmp(k,iCell)   = 0._RKIND
       acdqimp(k,iCell)   = 0._RKIND
       acdqsmp(k,iCell)   = 0._RKIND
       acdqgmp(k,iCell)   = 0._RKIND
       acdtpbl(k,iCell)   = 0._RKIND
       acdqvpbl(k,iCell)  = 0._RKIND
       acdqcpbl(k,iCell)  = 0._RKIND
       acdqipbl(k,iCell)  = 0._RKIND
       acdtradlw(k,iCell) = 0._RKIND
       acdtradsw(k,iCell) = 0._RKIND

       acqci(k,iCell) = 0._RKIND
       acqri(k,iCell) = 0._RKIND
       acqii(k,iCell) = 0._RKIND
       acqsi(k,iCell) = 0._RKIND
       acqgi(k,iCell) = 0._RKIND

       i_acthi(k,iCell) = 0
    enddo
 enddo

 if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
 if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       press_in(iCell,kk) = (pressure_p(k,iCell)+pressure_b(k,iCell)) / 100._RKIND
    enddo
 enddo

!interpolation to fixed pressure levels for fields located at cells centers and
!at mass points:
 if(.not.allocated(field_interp)) allocate(field_interp(nCells,nIntP))
 if(.not.allocated(press_interp)) allocate(press_interp(nCells,nIntP))
 do k = 1, nIntP
    kk = nIntP+1-k
    do iCell = 1, nCells
       press_interp(iCell,kk) = pFixed(k)
    enddo
 enddo

!... interpolation of the initial potential temperature:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = theta_m(k,iCell)/(1._RKIND+R_v/R_d*max(0._RKIND,qv(k,iCell)))
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       acthi(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo
!... interpolation of the initial water vapor mixing ratio:
 do k = 1, nVertLevels
    kk = nVertLevels+1-k
    do iCell = 1, nCells
       field_in(iCell,kk) = qv(k,iCell)
    enddo
 enddo
 call interp_tofixed_pressure(nCells,nVertLevels,nIntP,press_in,field_in,press_interp,field_interp)
 do k = 1, nIntP
    kk = nIntp+1-k
    do iCell = 1, nCells
       acqvi(k,iCell) = field_interp(iCell,kk)
    enddo
 enddo

!write(0,*)
!do iCell = 1, 100
!   do k = 1, nIntP
!      write(0,202) iCell,k,pFixed(k),acthi(k,iCell),acqvi(k,iCell)
!   enddo
!   write(0,*)
!enddo
!202 format(2i7,10(1x,e15.8))

!deallocate local arrays:
 if(allocated(press_in)) deallocate(press_in)
 if(allocated(field_in)) deallocate(field_in)

 if(allocated(field_interp)) deallocate(field_interp)
 if(allocated(press_interp)) deallocate(press_interp)

 write(0,*) '--- end subroutine init_pFixed'

 end subroutine init_pFixed

!==================================================================================================
 subroutine interp_tofixed_pressure(ncol,nlev_in,nlev_out,pres_in,field_in,pres_out,field_out)
!==================================================================================================

!input arguments:
 integer,intent(in):: ncol,nlev_in,nlev_out

 real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: pres_in,field_in
 real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: pres_out

!output arguments:
 real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out

!local variables:
! integer:: i1,i2,icol,k,kk
 integer:: icol,k,kk
 integer:: kkstart,kount
 integer,dimension(ncol):: kupper

 real(kind=RKIND):: dpl,dpu

!--------------------------------------------------------------------------------------------------

!formats:
! 201 format(i5,8(1x,e15.8))

!write(0,*)
!write(0,*) '--- enter subroutine interp_tofixed_pressure:'
!write(0,*) '... ncol     = ',ncol
!write(0,*) '... nlev_in  = ',nlev_in
!write(0,*) '... nlev_out = ',nlev_out
!i1=1 ; i2=ncol
!do k = 1, nlev_in
!   write(0,201) k,pres_in(i1,k),field_in(i1,k),pres_in(i2,k),field_in(i2,k)
!enddo
!write(0,*)

 do icol = 1, ncol
    kupper(icol) = 1
 enddo

 do k = 1, nlev_out 

    kkstart = nlev_in
    do icol = 1, ncol
       kkstart = min0(kkstart,kupper(icol))
    enddo
    kount = 0

    do kk = kkstart, nlev_in-1
       do icol = 1, ncol
          if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
             kupper(icol) = kk
             kount = kount + 1
!            write(0,201) kupper(icol),pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)
          endif
       enddo

       if(kount.eq.ncol) then
          do icol = 1, ncol
             dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
             dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
             field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                               + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
          end do
          goto 35
        end if
    enddo

    do icol = 1, ncol
       if(pres_out(icol,k) .lt. pres_in(icol,1)) then
          field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)
       elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in)) then
          field_out(icol,k) = field_in(icol,nlev_in)
       else
          dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
          dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
          field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                            + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
       endif
    enddo

 35 continue
!   write(0,201) kupper(i1),pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1), &
!                field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)
!   write(0,201) kupper(i2),pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1), &
!                field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)

 enddo

 end subroutine interp_tofixed_pressure

!==================================================================================================
 end module mpas_atmphys_interp_diagnostics
!==================================================================================================


