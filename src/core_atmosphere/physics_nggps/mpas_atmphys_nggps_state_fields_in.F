! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_nggps_state_fields_in
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use mpas_atmphys_nggps_constants
 use mpas_atmphys_nggps_vars
 use machine
 use nuopc_physics

 implicit none
 private
 public:: mpas_state_fields_in

 type,extends(state_fields_in):: mpas_state_fields_in

    contains
       procedure:: mpas_set       => mpas_state_fld_in_set
       procedure:: allocate_mem   => mpas_state_fld_in_allocate
       procedure:: deallocate_mem => mpas_state_fld_in_deallocate

 end type


!modifications to sourcecode:
!----------------------------
! * in subroutine mpas_state_fld_in_set, corrected the calculation of vvel. vvel is defined at the interface
!   between layers and is in Pa per second.
!   Laura D. Fowler (laura@ucar.edu) / 2015-12-01. 
! * commented out some print statements.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-14.
! * in the GFS physics, total cloud condensate and ozone are "specific" mixing ratios and not "strict" mixing
!   ratios. the arrays qgrs and tracers are corrected accordingly.
!   Laura D. Fowler (laura@ucar.edu) / 2016-01-24.

 contains


!=================================================================================================================
 subroutine mpas_state_fld_in_allocate(this,mesh)
!=================================================================================================================

 type(mpas_pool_type),intent(in):: mesh
 class(mpas_state_fields_in),intent(inout):: this

!local pointers:
 integer,pointer:: nCellsSolve
 integer,pointer:: nVertLevels

 integer:: ntrac

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_state_fld_allocate:'

 call mpas_pool_get_dimension(mesh,"nCellsSolve",nCellsSolve)
 call mpas_pool_get_dimension(mesh,"nVertLevels",nVertLevels)

 ntrac = mdl_parameters%ntrac

!... FOR PHYSICS AND RADIATION:
 if(.not.allocated(prsi)     ) allocate(prsi(1:nCellsSolve,1:nVertLevels+1)      )
 if(.not.allocated(prsl)     ) allocate(prsl(1:nCellsSolve,1:nVertLevels)        )
 if(.not.allocated(prslk)    ) allocate(prslk(1:nCellsSolve,1:nVertLevels)       )
 if(.not.allocated(tgrs)     ) allocate(tgrs(1:nCellsSolve,1:nVertLevels)        )

!... FOR RADIATION ONLY:
 if(.not.allocated(qgrs_rad) ) allocate(qgrs_rad(1:nCellsSolve,1:nVertLevels)    ) 
 if(.not.allocated(vvl)      ) allocate(vvl(1:nCellsSolve,1:nVertLevels)         )
 if(.not.allocated(tracer)   ) allocate(tracer(1:nCellsSolve,1:nVertLevels,ntrac))

!... FOR PHYSICS ONLY:
 if(.not.allocated(adjtrc)   ) allocate(adjtrc(ntrac)                            ) 
 if(.not.allocated(pgr)      ) allocate(pgr(1:nCellsSolve)                       )
 if(.not.allocated(ugrs)     ) allocate(ugrs(1:nCellsSolve,1:nVertLevels)        )
 if(.not.allocated(vgrs)     ) allocate(vgrs(1:nCellsSolve,1:nVertLevels)        )
 if(.not.allocated(prsik)    ) allocate(prsik(1:nCellsSolve,1:nVertLevels+1)     )
 if(.not.allocated(phii)     ) allocate(phii(1:nCellsSolve,1:nVertLevels+1)      )
 if(.not.allocated(phil)     ) allocate(phil(1:nCellsSolve,1:nVertLevels)        )
 if(.not.allocated(qgrs)     ) allocate(qgrs(1:nCellsSolve,1:nVertLevels,ntrac)  )

!... additional arrays needed for the initialization of the state_fields_in arrays:
 if(.not.allocated(rhos)     )  allocate(rhos(1:nCellsSolve,1:nVertLevels)       )
 if(.not.allocated(dzgrid)   ) allocate(dzgrid(1:nCellsSolve,1:nVertLevels)      )
 if(.not.allocated(fzm)      ) allocate(fzm(1:nCellsSolve,1:nVertLevels)         )
 if(.not.allocated(fzp)      ) allocate(fzp(1:nCellsSolve,1:nVertLevels)         )
 if(.not.allocated(tgrsi)    ) allocate(tgrsi(1:nCellsSolve,1:nVertLevels+1)     )
 if(.not.allocated(rhosi)    ) allocate(rhosi(1:nCellsSolve,1:nVertLevels+1)     )

!... hydrostatic pressure sounding:
 if(.not.allocated(pgr_hyd)  ) allocate(pgr_hyd(1:nCellsSolve)                   )
 if(.not.allocated(prsl_hyd) ) allocate(prsl_hyd(1:nCellsSolve,1:nVertLevels)    )
 if(.not.allocated(prsi_hyd) ) allocate(prsi_hyd(1:nCellsSolve,1:nVertLevels+1)  )

!write(0,*) '--- end subroutine mpas_state_fld_allocate'

 end subroutine mpas_state_fld_in_allocate

!=================================================================================================================
 subroutine mpas_state_fld_in_deallocate(this)
!=================================================================================================================

 class(mpas_state_fields_in),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
!write(0,*)
!write(0,*) '--- enter subroutine mpas_state_fld_in_deallocate:'

!... FOR PHYSICS AND RADIATION:
 if(allocated(prsi)    ) deallocate(prsi    )
 if(allocated(prsl)    ) deallocate(prsl    )
 if(allocated(prslk)   ) deallocate(prslk   )
 if(allocated(tgrs)    ) deallocate(tgrs    )

!... FOR RADIATION ONLY:
 if(allocated(qgrs_rad)) deallocate(qgrs_rad) 
 if(allocated(vvl)     ) deallocate(vvl     )
 if(allocated(tracer)  ) deallocate(tracer  )

!... FOR PHYSICS ONLY:
 if(allocated(pgr)     ) deallocate(pgr     )
 if(allocated(adjtrc)  ) deallocate(adjtrc  ) 

 if(allocated(ugrs)    ) deallocate(ugrs    )
 if(allocated(vgrs)    ) deallocate(vgrs    )
 if(allocated(prsik)   ) deallocate(prsik   )
 if(allocated(phii)    ) deallocate(phii    )
 if(allocated(phil)    ) deallocate(phil    )
 if(allocated(qgrs)    ) deallocate(qgrs    )

!... additional arrays needed for the initialization of the state_fields_in arrays:
 if(allocated(rhos)    ) deallocate(rhos    )
 if(allocated(dzgrid)  ) deallocate(dzgrid  )
 if(allocated(fzm)     ) deallocate(fzm     )
 if(allocated(fzp)     ) deallocate(fzp     )
 if(allocated(tgrsi)   ) deallocate(tgrsi   )
 if(allocated(rhosi)   ) deallocate(rhosi   )

!... hydrostatic pressure sounding:
 if(allocated(pgr_hyd) ) deallocate(pgr_hyd )
 if(allocated(prsl_hyd)) deallocate(prsl_hyd)
 if(allocated(prsi_hyd)) deallocate(prsi_hyd)

!... pressure for radiation parameterizations:
!if(.not.allocated(pgr_rad)  ) allocate(pgr_rad(1:nCellsSolve)                   )
!if(.not.allocated(prsl_rad) ) allocate(prsl_rad(1:nCellsSolve,1:nVertLevels)    )
!if(.not.allocated(prsi_rad) ) allocate(prsi_rad(1:nCellsSolve,1:nVertLevels+1)  )

!if(.not.allocated(qtrcs)    ) allocate(qtrcs(1:nCellsSolve,1:nVertLevels,ntrac-1))
!if(.not.allocated(qg0)      ) allocate(qg0(1:nCellsSolve,1:nVertLevels,ntrac)    )

!write(0,*) '--- end subroutine mpas_state_fld_in_deallocate'

 end subroutine mpas_state_fld_in_deallocate

!=================================================================================================================
 subroutine mpas_state_fld_in_set(this,mesh,state,time_lev,diag)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag

 integer,intent(in):: time_lev

!inout variables:
 class(mpas_state_fields_in),intent(inout):: this

!local variables:
 integer:: k,kk,iCell,i1,i2,nt,it,izgrid,n,ntcw,ntoz,ntrac
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: index_qv,index_qc,index_qo3

 real(kind=RKIND),dimension(:,:),pointer  :: uReconstructZonal,uReconstructMeridional
 real(kind=RKIND),dimension(:,:),pointer  :: zz,zgrid,pressure_p,pressure_b,exner
 real(kind=RKIND),dimension(:,:),pointer  :: rho,theta_m,w
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qo3
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 
 real(kind=RKIND):: tem1,tem2,rho1,rho2,z0,z1,z2,w1,w2
 real(kind=RKIND):: sph

!-----------------------------------------------------------------------------------------------------------------
!write(0,*) '--- enter subroutine mpas_state_fld_setphys_in:'

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'zz'   ,zz   )
 call mpas_pool_get_array(mesh,'zgrid',zgrid)

 call mpas_pool_get_array(diag,'exner'                 ,exner                 )
 call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p            )
 call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b            )
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,uReconstructZonal     )
 call mpas_pool_get_array(diag,'uReconstructMeridional',uReconstructMeridional)

 call mpas_pool_get_array(state,'rho_zz' ,rho    ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w'      ,w      ,time_lev)
 
 call mpas_pool_get_dimension(state,'index_qv' ,index_qv)
 call mpas_pool_get_dimension(state,'index_qc' ,index_qc)
 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv  => scalars(index_qv,:,:)
 qc  => scalars(index_qc,:,:)
 qo3 => scalars(index_qo3,:,:)
 
 ntcw  = mdl_parameters%ntcw
 ntoz  = mdl_parameters%ntoz
 ntrac = mdl_parameters%ntrac

!... geopotential height and layer thickness:
 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    phil(iCell,k)   = 0._kind_phys
    phii(iCell,k)   = 0._kind_phys
    dzgrid(iCell,k) = 0._kind_phys
 enddo
 enddo
 k = nVertLevels+1
 do iCell = 1, nCellsSolve
    phii(iCell,k) = 0._kind_phys
 enddo

 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    dzgrid(iCell,k) = zgrid(k+1,iCell) - zgrid(k,iCell)
    phii(iCell,k+1) = phii(iCell,k) + gravity*dzgrid(iCell,k)
 enddo
 enddo

 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    phil(iCell,k) = 0.5_kind_phys*(phii(iCell,k) + phii(iCell,k+1))
 enddo
 enddo

!... layer mean temperature, layer mean pressure, layer mean exner function, and layer mean
!... geopotential height:
 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    rhos(iCell,k)  = rho(k,iCell) * zz(k,iCell) * (1._kind_phys + qv(k,iCell)) 
    tgrs(iCell,k)  = theta_m(k,iCell) * exner(k,iCell) / (1. + R_v/R_d * qv(k,iCell))
    prsl(iCell,k)  = pressure_b(k,iCell) + pressure_p(k,iCell)
    prslk(iCell,k) = exner(k,iCell)
 enddo
 enddo
!write(0,*) '--- end rhos,tgrs,prsl,prslk,phil,dzgrid:'
!do iCell = 1, nCellsSolve
!   do k = 1, nVertLevels
!      write(0,201) iCell,k,prsl(iCell,k)/100.,prslk(iCell,k),phil(iCell,k),tgrs(iCell,k)
!   enddo
!   write(0,*)
!enddo
!stop

!... temperature,pressure,exner function,and geopotential height at layer interfaces:
 do k = 2, nVertLevels
 do iCell = 1, nCellsSolve
    tem1 = 1./(zgrid(k+1,iCell)-zgrid(k-1,iCell))
    fzm(iCell,k) = (zgrid(k,iCell)-zgrid(k-1,iCell)) * tem1
    fzp(iCell,k) = (zgrid(k+1,iCell)-zgrid(k,iCell)) * tem1

    tgrsi(iCell,k) = fzm(iCell,k)*tgrs(iCell,k) + fzp(iCell,k)*tgrs(iCell,k-1)
    rhosi(iCell,k) = fzm(iCell,k)*rhos(iCell,k) + fzp(iCell,k)*rhos(iCell,k-1)
    prsi(iCell,k)  = fzm(iCell,k)*prsl(iCell,k) + fzp(iCell,k)*prsl(iCell,k-1)
 enddo
 enddo
!... extrapolation of temperature and pressure to the top-of-the-model:
 k = nVertLevels+1
 do iCell = 1, nCellsSolve
    z0 = zgrid(k,iCell)
    z1 = 0.5_kind_phys*(zgrid(k,iCell)+zgrid(k-1,iCell)) 
    z2 = 0.5_kind_phys*(zgrid(k-1,iCell)+zgrid(k-2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    tgrsi(iCell,k) = w1*tgrs(iCell,k-1) + w2*tgrs(iCell,k-2)
    rhosi(iCell,k) = w1*rhos(iCell,k-1) + w2*rhos(iCell,k-2)
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    prsi(iCell,k) = exp(w1*log(prsl(iCell,k-1)) + w2*log(prsl(iCell,k-2)))
 enddo
!... extrapolation of temperature and pressure to the surface; calculation of the surface pressure:
 k = 1
 do iCell = 1, nCellsSolve
    z0 = zgrid(k,iCell)
    z1 = 0.5*(zgrid(k,iCell)+zgrid(k+1,iCell)) 
    z2 = 0.5*(zgrid(k+1,iCell)+zgrid(k+2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    tgrsi(iCell,k) = w1*tgrs(iCell,k) + w2*tgrs(iCell,k+1)
    rhosi(iCell,k) = w1*rhos(iCell,k) + w2*rhos(iCell,k+1)
!   prsi(iCell,k)  = w1*prsl(iCell,k) + w2*prsl(iCell,k+1)

    tem1 = zgrid(2,iCell)-zgrid(1,iCell)
    tem2 = zgrid(3,iCell)-zgrid(2,iCell)
    rho1 = rho(1,iCell)*zz(1,iCell)*(1._RKIND + qv(1,iCell))
    rho2 = rho(2,iCell)*zz(2,iCell)*(1._RKIND + qv(2,iCell))
    prsi(iCell,k) = 0.5*gravity*(zgrid(2,iCell)-zgrid(1,iCell)) &
                  * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
    prsi(iCell,k) = prsi(iCell,k) + prsl(iCell,k)
 enddo
!write(0,*) '--- end prsi,phii,tgrsi:'
!do iCell = 1, nCellsSolve
!   do k = 1, nVertLevels
!      write(0,201) iCell,k,prsl(iCell,k)/100.,prsi(iCell,k)/100.,phil(iCell,k),phii(iCell,k),&
!                   tgrs(iCell,k),tgrsi(iCell,k)
!   enddo
!   write(0,*)
!enddo

!... surface pressure and exner function at interface between layers:
 do iCell = 1,nCellsSolve
    pgr(iCell)    = prsi(iCell,1)
    do k = 1,nVertLevels+1
       prsik(iCell,k) = (prsi(iCell,k)/pgr(iCell))**rcp
    enddo
 enddo

!... check that the pressure at the surface and layer adjacent next to the surface are not negative:
 do k = 1,nVertLevels
    do iCell = 1,nCellsSolve
       if(prsl(iCell,k) .le. 0._RKIND .or. prsi(iCell,k) .le. 0._RKIND) then
          do kk = 1, nVertLevels
             write(0,202) iCell,k,kk,prsl(iCell,kk),prsi(iCell,kk)
          enddo
       endif
    enddo
 enddo

!... zonal and meridional winds, and vertical velocities:
 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    ugrs(iCell,k) = uReconstructZonal(k,iCell)
    vgrs(iCell,k) = uReconstructMeridional(k,iCell)
    vvl(iCell,k) = - gravity * rhosi(iCell,k) * w(k,iCell)
 enddo
 enddo
!write(0,*) '--- end ugrs, vgrs, vvl:'
!do iCell = 1, nCellsSolve
!   do k = 1, nVertLevels
!      write(0,201) iCell,k,w(k,iCell),rhosi(iCell,k),vvl(iCell,k),rhos(iCell,k), &
!                   ugrs(iCell,k),vgrs(iCell,k)
!   enddo
!   k = nVertLevels+1
!   write(0,201) iCell,k,w(k,iCell)
!   write(0,*)
!enddo
!stop

!... layer specific humidity and "tracer" mixing ratios used in radiation, and mixing ratios used in other
!... physics parameterizations:
 do k = 1, ntrac
    adjtrc(k) = 1._RKIND
 enddo

 do k = 1, nVertLevels
 do iCell = 1,nCellsSolve
    tracer(iCell,k,:) = 0._RKIND
    qgrs(iCell,k,:)   = 0._RKIND
 enddo
 enddo

 do k = 1, nVertLevels
    do iCell = 1, nCellsSolve
       sph = 1._RKIND/(1._RKIND + qv(k,iCell))

       !... for radiation:
       qgrs_rad(iCell,k) = qv(k,iCell)  * sph
       tracer(iCell,k,1) = qo3(k,iCell) * sph
       tracer(iCell,k,2) = qc(k,iCell)  * sph
       qgrs_rad(iCell,k) = max(qmin,qgrs_rad(iCell,k))

       !... for other physics:
       qgrs(iCell,k,1)    = qv(k,iCell)  * sph
       qgrs(iCell,k,ntoz) = qo3(k,iCell) * sph
       qgrs(iCell,k,ntcw) = qc(k,iCell)  * sph
    enddo
 enddo
!write(0,*) '--- end qgrs_rad, qgrs, and tracers:'
!do iCell = 1, 50
!   do k = 1, nVertLevels
!      write(0,201) iCell,k,qgrs(iCell,k,1),qgrs(iCell,k,ntoz),qgrs(iCell,k,ntcw),&
!         qgrs_rad(iCell,k),tracer(iCell,k,1),tracer(iCell,k,2)
!   enddo
!   write(0,*)
!enddo

!... fills type "mpas_state_fields_in":
 call this%setphys(prsi , prsl , prslk , tgrs , qgrs , vvl , pgr, &
                   ugrs , vgrs , prsik , phii , phil , adjtrc)

 call this%setrad(prsi,prsl,prslk,tgrs,qgrs_rad,tracer,vvl)

!write(0,*) '--- end subroutine mpas_state_fld_setphys_in'

!formats:
 201 format(2i5,10(1x,e15.8))
 202 format(3i5,10(1x,e15.8))

 end subroutine mpas_state_fld_in_set

!=================================================================================================================
 end module mpas_atmphys_nggps_state_fields_in
!=================================================================================================================
