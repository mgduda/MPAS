! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_gfs_interface
 use mpas_grid_types

 use mpas_atmphys_gfs_constants, grav => g
 use mpas_atmphys_gfs_vars

 implicit none
 private
 public:: allocate_forall_physics,   &
          deallocate_forall_physics, &
          MPAS_to_physics


!>\brief interface for conversion between variables used in the MPAS dynamical core and variables
!> needed in the GFS physics parameterizations.
!>\author Laura D. Fowler (send comments to laura@mmm.ucar.edu).
!>\date 2013-11-21.
!>
!>\details
!>
!> subroutines in mpas_atmphys_interface_nhyd:
!> -------------------------------------------
!> allocate_forall_physics  : allocation of all meteorological variables needed in the physics.
!> deallocate_forall_physics: deallocation of all meteorological variables needed in the physics.
!> MPAS_to_physics          : conversion of input "MPAS" variables to "physics" variables.
!> add-ons and modifications to sourcecode:
!> ----------------------------------------
!>    * Added diagnostic for the surface pressure.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-03-20.
!>    * Modified sourcecode to use pools.
!>      Laura D. Fowler (birch.mmm.ucar.edu) / 2014-06-06.


 contains

!=================================================================================================================
 subroutine allocate_forall_physics(mesh)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!local pointers:
 integer,pointer:: nCellsSolve,nVertLevels

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 if(.not.allocated(dumphys)  ) allocate(dumphys(1:nCellsSolve) )
 
 if(.not.allocated(pgr)      ) allocate(pgr(1:nCellsSolve) )
 if(.not.allocated(fzm)      ) allocate(fzm(1:nCellsSolve,1:nVertLevels)       )
 if(.not.allocated(fzp)      ) allocate(fzp(1:nCellsSolve,1:nVertLevels)       )
 if(.not.allocated(ugrs)     ) allocate(ugrs(1:nCellsSolve,1:nVertLevels)      )
 if(.not.allocated(vgrs)     ) allocate(vgrs(1:nCellsSolve,1:nVertLevels)      )
 if(.not.allocated(vvel)     ) allocate(vvel(1:nCellsSolve,1:nVertLevels)      )

 if(.not.allocated(dzgrid)   ) allocate(dzgrid(1:nCellsSolve,1:nVertLevels)    )
 if(.not.allocated(qvrs)     ) allocate(qvrs(1:nCellsSolve,1:nVertLevels)      )
 if(.not.allocated(rhos)     ) allocate(rhos(1:nCellsSolve,1:nVertLevels)      )
 if(.not.allocated(tgrs)     ) allocate(tgrs(1:nCellsSolve,1:nVertLevels)      )
 if(.not.allocated(prsl)     ) allocate(prsl(1:nCellsSolve,1:nVertLevels)      )
 if(.not.allocated(prslk)    ) allocate(prslk(1:nCellsSolve,1:nVertLevels)     )
 if(.not.allocated(phil)     ) allocate(phil(1:nCellsSolve,1:nVertLevels)      )

 if(.not.allocated(tgrsi)    ) allocate(tgrsi(1:nCellsSolve,1:nVertLevels+1)   )
 if(.not.allocated(prsi)     ) allocate(prsi(1:nCellsSolve,1:nVertLevels+1)    )
 if(.not.allocated(prsik)    ) allocate(prsik(1:nCellsSolve,1:nVertLevels+1)   )
 if(.not.allocated(phii)     ) allocate(phii(1:nCellsSolve,1:nVertLevels+1)    )

!... hydrostatic pressure sounding:
 if(.not.allocated(pgr_hyd)  ) allocate(pgr_hyd(1:nCellsSolve) )
 if(.not.allocated(prsl_hyd) ) allocate(prsl_hyd(1:nCellsSolve,1:nVertLevels)  )
 if(.not.allocated(prsi_hyd) ) allocate(prsi_hyd(1:nCellsSolve,1:nVertLevels+1))

!... pressure for radiation parameterizations:
 if(.not.allocated(pgr_rad)  ) allocate(pgr_rad(1:nCellsSolve) )
 if(.not.allocated(prsl_rad) ) allocate(prsl_rad(1:nCellsSolve,1:nVertLevels)  )
 if(.not.allocated(prsi_rad) ) allocate(prsi_rad(1:nCellsSolve,1:nVertLevels+1))
 
!... updated arrays after call to gbphys:
 if(.not.allocated(gt0)      ) allocate(gt0(1:nCellsSolve,1:nVertLevels)       )
 if(.not.allocated(gu0)      ) allocate(gu0(1:nCellsSolve,1:nVertLevels)       )
 if(.not.allocated(gv0)      ) allocate(gv0(1:nCellsSolve,1:nVertLevels)       )

 end subroutine allocate_forall_physics

!=================================================================================================================
 subroutine deallocate_forall_physics
!=================================================================================================================

 if(allocated(dumphys) ) deallocate(dumphys   )

 if(allocated(pgr)      ) deallocate(pgr      )
 if(allocated(fzm)      ) deallocate(fzm      )
 if(allocated(fzp)      ) deallocate(fzp      )
 if(allocated(ugrs)     ) deallocate(ugrs     )
 if(allocated(vgrs)     ) deallocate(vgrs     )
 if(allocated(vvel)     ) deallocate(vvel     )

 if(allocated(dzgrid)   ) deallocate(dzgrid   )
 if(allocated(qvrs)     ) deallocate(qvrs     )
 if(allocated(rhos)     ) deallocate(rhos     )
 if(allocated(tgrs)     ) deallocate(tgrs     )
 if(allocated(prsl)     ) deallocate(prsl     )
 if(allocated(prslk)    ) deallocate(prslk    )
 if(allocated(phil)     ) deallocate(phil     )

 if(allocated(tgrsi)    ) deallocate(tgrsi    )
 if(allocated(prsi)     ) deallocate(prsi     )
 if(allocated(prsik)    ) deallocate(prsik    )
 if(allocated(phii)     ) deallocate(phii     )

 if(allocated(qgrs)     ) deallocate(qgrs     )
 if(allocated(qtrcs)    ) deallocate(qtrcs    )

!... hydrostatic pressure sounding:
 if(allocated(pgr_hyd)  ) deallocate(pgr_hyd  )
 if(allocated(prsl_hyd) ) deallocate(prsl_hyd )
 if(allocated(prsi_hyd) ) deallocate(prsi_hyd )

!... pressure for radiation parameterizations:
 if(allocated(pgr_rad)  ) deallocate(pgr_rad  )
 if(allocated(prsl_rad) ) deallocate(prsl_rad )
 if(allocated(prsi_rad) ) deallocate(prsi_rad )
 
!... updated arrays after call to gbphys:
 if(allocated(gt0)      ) deallocate(gt0      )
 if(allocated(gu0)      ) deallocate(gu0      )
 if(allocated(gv0)      ) deallocate(gv0      )

 if(allocated(gq0)      ) deallocate(gq0      )

 end subroutine deallocate_forall_physics

!=================================================================================================================
 subroutine MPAS_to_physics(mesh,state,time_lev,diag)
!=================================================================================================================

!input variables:
 integer,intent(in):: time_lev
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state

!inout variables:
 type(mpas_pool_type),intent(inout):: diag

!local pointers:
 integer,pointer:: nCells,nCellsSolve,nVertLevels
 integer,pointer:: index_qv,index_qc,index_qo3

 real(kind=RKIND),dimension(:),pointer    :: surface_pressure
 real(kind=RKIND),dimension(:,:),pointer  :: zz,zgrid,pressure_p,pressure_b,exner
 real(kind=RKIND),dimension(:,:),pointer  :: rho,theta_m,w
 real(kind=RKIND),dimension(:,:),pointer  :: qv,qc,qo3
 real(kind=RKIND),dimension(:,:),pointer  :: uReconstructZonal,uReconstructMeridional
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!local variables:
 integer:: k,kk,iCell,i1,i2,nt,it,izgrid

 real(kind=RKIND):: tem1,z0,z1,z2,w1,w2
 real(kind=RKIND),dimension(:,:),allocatable:: height

!-----------------------------------------------------------------------------------------------------------------
 write(0,*)
 write(0,*) '--- enter subroutine MPAS_to_physics:'

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(mesh,'zz'   ,zz   )
 call mpas_pool_get_array(mesh,'zgrid',zgrid)

 call mpas_pool_get_array(diag,'exner'                 ,exner                 )
 call mpas_pool_get_array(diag,'pressure_p'            ,pressure_p            )
 call mpas_pool_get_array(diag,'pressure_base'         ,pressure_b            )
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,uReconstructZonal     )
 call mpas_pool_get_array(diag,'uReconstructMeridional',uReconstructMeridional)
 call mpas_pool_get_array(diag,'surface_pressure'      ,surface_pressure      )

 call mpas_pool_get_array(state,'rho_zz' ,rho    ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)
 call mpas_pool_get_array(state,'w'      ,w      ,time_lev)
 
 call mpas_pool_get_dimension(state,'index_qv' ,index_qv )
 call mpas_pool_get_dimension(state,'index_qc' ,index_qc )
 call mpas_pool_get_dimension(state,'index_qo3',index_qo3)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv  => scalars(index_qv,:,:)
 qc  => scalars(index_qc,:,:)
 qo3 => scalars(index_qo3,:,:)

 if(.not.allocated(height)) allocate(height(nVertLevels+1,nCellsSolve))
 do k = 1, nVertLevels + 1
 do iCell = 1, nCellsSolve
!    height(k,iCell) = max(0._RKIND,zgrid(k,iCell))
    height(k,iCell) = zgrid(k,iCell)
 enddo
 enddo

!... zonal and meridional winds, and vertical velocities:
 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    ugrs(iCell,k) = uReconstructZonal(k,iCell)
    vgrs(iCell,k) = uReconstructMeridional(k,iCell)
    vvel(iCell,k) = 0.5_RKIND * (w(k,iCell) + w(k+1,iCell))
    vvel(iCell,k) = - grav * rho(k,iCell) * zz(k,iCell) * vvel(iCell,k)
 enddo
 enddo

!... layer mean temperature, layer mean pressure, layer mean exner function, and layer mean
!... geopotential height:
 do k = 1, nVertLevels
 do iCell = 1, nCellsSolve
    qvrs(iCell,k)   = qv(k,iCell)
    rhos(iCell,k)   = rho(k,iCell) * zz(k,iCell) * (1._RKIND + qv(k,iCell)) 
    tgrs(iCell,k)   = theta_m(k,iCell) * exner(k,iCell) / (1. + R_v/R_d * qv(k,iCell))
    prsl(iCell,k)   = pressure_b(k,iCell) + pressure_p(k,iCell)
    prslk(iCell,k)  = exner(k,iCell)
    phil(iCell,k)   = 0.5_RKIND * grav * (height(k,iCell)+height(k+1,iCell))
    dzgrid(iCell,k) = height(k+1,iCell) - height(k,iCell) 
 enddo
 enddo
 k = 1
 do iCell = 1, nCellsSolve
    surface_pressure(iCell) = prsl(iCell,k)
 enddo

!... temperature, pressure, exner function, and geopotential height at layer interfaces:
 do k = 2, nVertLevels
 do iCell = 1, nCellsSolve
    tem1 = 1./(height(k+1,iCell)-height(k-1,iCell))
    fzm(iCell,k) = (height(k,iCell)-height(k-1,iCell)) * tem1
    fzp(iCell,k) = (height(k+1,iCell)-height(k,iCell)) * tem1

    tgrsi(iCell,k) = fzm(iCell,k)*tgrs(iCell,k) + fzp(iCell,k)*tgrs(iCell,k-1)
    prsi(iCell,k)  = fzm(iCell,k)*prsl(iCell,k) + fzp(iCell,k)*prsl(iCell,k-1)
    phii(iCell,k)  = grav * height(k,iCell)
 enddo
 enddo
!... extrapolation of temperature and pressure to the surface:
 k = 1
 do iCell = 1, nCellsSolve
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k+1,iCell)) 
    z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    tgrsi(iCell,k) = w1*tgrs(iCell,k) + w2*tgrs(iCell,k+1)
    prsi(iCell,k)  = w1*prsl(iCell,k) + w2*prsl(iCell,k+1)
    phii(iCell,k)  = grav * height(k,iCell)
 enddo
!... extrapolation of temperature and pressure to the top-of-the-model:
 k = nVertLevels+1
 do iCell = 1, nCellsSolve
    z0 = height(k,iCell)
    z1 = 0.5*(height(k,iCell)+height(k-1,iCell)) 
    z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
    w1 = (z0-z2)/(z1-z2)
    w2 = 1.-w1
    tgrsi(iCell,k) = w1*tgrs(iCell,k-1) + w2*tgrs(iCell,k-2)
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    prsi(iCell,k) = exp(w1*log(prsl(iCell,k-1)) + w2*log(prsl(iCell,k-2)))
    phii(iCell,k) = grav * height(k,iCell)
 enddo
!... surface pressure and exner function:
 do iCell = 1,nCellsSolve
    do k = 1,nVertLevels
       if(prsl(iCell,k) .le. 0._RKIND .or. prsi(iCell,k) .le. 0._RKIND) then
          do kk = 1, nVertLevels
             write(0,201) iCell,kk,fzm(iCell,kk),fzp(iCell,kk), &
                prsl(iCell,kk),prsi(iCell,kk)
          enddo
       endif
    enddo
 enddo

 do iCell = 1,nCellsSolve
    pgr(iCell) = prsi(iCell,1)
    do k = 1,nVertLevels
       prslk(iCell,k) = (prsl(iCell,k)/pgr(iCell))**rcp
       prsik(iCell,k) = (prsi(iCell,k)/pgr(iCell))**rcp
    enddo
    k = nVertLevels+1
    prsik(iCell,k) = (prsi(iCell,k)/pgr(iCell))**rcp
 enddo

!... calculation of the hydrostatic pressure:
 do iCell = 1,nCellsSolve
    k = nVertLevels+1
    prsi_hyd(iCell,k) = prsi(iCell,k)
    do k = nVertLevels,1,-1
       prsi_hyd(iCell,k) = prsi_hyd(iCell,k+1) + grav * rhos(iCell,k) * (height(k+1,iCell)-height(k,iCell))
    enddo
    !surface pressure:
    pgr_hyd(iCell) = prsi_hyd(iCell,1)

    !hydrostatic pressure at theta points:
    do k = 1, nVertLevels
       prsl_hyd(iCell,k) = 0.5_RKIND * (prsi_hyd(iCell,k) + prsi_hyd(iCell,k+1))
    enddo 
 enddo

!... conversion of all pressures from Pa to centibars for radiation parameterizations:
#if defined(do_hydrostatic_pressure)
 do iCell = 1, nCellsSolve
    pgr_rad(iCell) = pgr_hyd(iCell) * 0.001_RKIND

    do k = 1, nVertLevels
       prsl_rad(iCell,k) = prsl_hyd(iCell,k) * 0.001_RKIND
       prsi_rad(iCell,k) = prsi_hyd(iCell,k) * 0.001_RKIND
    enddo

    k = nVertLevels+1
    prsi_rad(iCell,k) = prsi_hyd(iCell,k) * 0.001_RKIND
 enddo
#else
 do iCell = 1, nCellsSolve
    pgr_rad(iCell) = pgr(iCell) * 0.001_RKIND

    do k = 1, nVertLevels
       prsl_rad(iCell,k) = prsl(iCell,k) * 0.001_RKIND
       prsi_rad(iCell,k) = prsi(iCell,k) * 0.001_RKIND
    enddo

    k = nVertLevels+1
    prsi_rad(iCell,k) = prsi(iCell,k) * 0.001_RKIND
 enddo
#endif

 if(allocated(height)) deallocate(height)

!... fill the tracer array:
 if(.not. allocated(qgrs) ) allocate(qgrs(1:nCellsSolve,1:nVertLevels,ntrac)   )
 if(.not. allocated(qtrcs)) allocate(qtrcs(1:nCellsSolve,1:nVertLevels,ntrac-1))
 if(.not. allocated(gq0)  ) allocate(gq0(1:nCellsSolve,1:nVertLevels,ntrac)    )

 do iCell = 1,nCellsSolve
    do k = 1,nVertLevels
       qgrs(iCell,k,1)    = qv(k,iCell)
       qgrs(iCell,k,2)    = qo3(k,iCell)
       qgrs(iCell,k,ntcw) = qc(k,iCell) 

       qtrcs(iCell,k,1)   = qo3(k,iCell)
       qtrcs(iCell,k,2)   = qc(k,iCell)
    enddo
 enddo

 write(0,*) '--- end subroutine MPAS_to_physics'

!formats:
 201 format(2i5,10(1x,e15.8))

 end subroutine MPAS_to_physics

!=================================================================================================================
 end module mpas_atmphys_gfs_interface
!=================================================================================================================
