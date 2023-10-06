!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> \brief Contains the definitions of the physical fluxes of the equation system.
!>
!> The routine EvalFlux3D will compute the advection (Euler) part only, and can be called either for a single point or for
!> a volume cell. The fluxes are computed in three spatial dimension - for 2D computations, the fluxes in the third dimension
!> will always be set to 0.
!> EvalDiffFlux3D will do the same thing, but compute only the diffusive part of the fluxes. Additionally, a routine to compute
!> the fluxes on a single side is provided (used in the riemann routines).
!> The EvalFlux1D routines are used in the Riemann solver, where only a flux in one spatial dimension is needed.
!>
!> The flux definitions are only done once in the single point routines, all other (side, volume) routines will simply wrap
!> to this definition.
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE EvalFlux3D
  MODULE PROCEDURE EvalFlux3D
  MODULE PROCEDURE EvalFlux3D_Elem
  MODULE PROCEDURE EvalFlux3D_Elem_CUDA
  MODULE PROCEDURE EvalFlux3D_Elems_CUDA
END INTERFACE

INTERFACE EvalEulerFlux1D
  MODULE PROCEDURE EvalEulerFlux1D
END INTERFACE

INTERFACE EvalEulerFlux1D_fast
  MODULE PROCEDURE EvalEulerFlux1D_fast
END INTERFACE

#if PARABOLIC
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D
  !MODULE PROCEDURE EvalDiffFlux3D_Point
  MODULE PROCEDURE EvalDiffFlux3D_Surface
  MODULE PROCEDURE EvalDiffFlux3D_Volume
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalFlux3D, EvalEulerFlux1D, EvalEulerFlux1D_fast,EvalTransformedFlux3D_Kernel
#if PARABOLIC
PUBLIC::EvalDiffFlux3D
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute advection part of the Navier-Stokes fluxes in all space dimensions using the conservative and primitive variables
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalFlux3D(U,UPrim,f,g,h)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U        !< Conservative solution
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: f,g,h    !> Physical fluxes in x/y/z direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Ep
!==================================================================================================================================
! auxiliary variables
Ep   = U(ENER) + UPrim(PRES)
#if PP_dim==3
! Euler part
! Euler fluxes x-direction
f(DENS) = U(MOM1)                             ! rho*u
f(MOM1) = U(MOM1) * UPrim(VEL1) + UPrim(PRES) ! rho*u²+p
f(MOM2) = U(MOM1) * UPrim(VEL2)               ! rho*u*v
f(MOM3) = U(MOM1) * UPrim(VEL3)               ! rho*u*w
f(ENER) = Ep * UPrim(VEL1)                    ! (rho*e+p)*u
! Euler fluxes y-direction
g(DENS) = U(MOM2)                             ! rho*v
g(MOM1) = f(MOM2)                             ! rho*u*v
g(MOM2) = U(MOM2) * UPrim(VEL2) + UPrim(PRES) ! rho*v²+p
g(MOM3) = U(MOM2) * UPrim(VEL3)               ! rho*v*w
g(ENER) = Ep * UPrim(VEL2)                    ! (rho*e+p)*v
! Euler fluxes z-direction
h(DENS) = U(MOM3)                             ! rho*v
h(MOM1) = f(MOM3)                             ! rho*u*w
h(MOM2) = g(MOM3)                             ! rho*v*w
h(MOM3) = U(MOM3) * UPrim(VEL3) + UPrim(PRES) ! rho*v²+p
h(ENER) = Ep * UPrim(VEL3)                    ! (rho*e+p)*w
#else

! Euler part
! Euler fluxes x-direction
f(DENS) = U(MOM1)                             ! rho*u
f(MOM1) = U(MOM1)*UPrim(VEL1)+UPrim(PRES)     ! rho*u²+p
f(MOM2) = U(MOM1)*UPrim(VEL2)                 ! rho*u*v
f(MOM3) = 0.
f(ENER) = Ep*UPrim(VEL1)                      ! (rho*e+p)*u
! Euler fluxes y-direction
g(DENS)= U(MOM2)                              ! rho*v
g(MOM1)= f(MOM2)                              ! rho*u*v
g(MOM2)= U(MOM2)*UPrim(VEL2)+UPrim(PRES)      ! rho*v²+p
g(MOM3)= 0.
g(ENER)= Ep*UPrim(VEL2)                       ! (rho*e+p)*v
! Euler fluxes z-direction
h   = 0.
#endif
END SUBROUTINE EvalFlux3D

!==================================================================================================================================
!> Compute advection part of the Navier-Stokes fluxes in all space dimensions using the conservative and primitive variables
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalTransformedFlux3D(U,UPrim,f,g,h,Mf,Mg,Mh)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U        !< Conservative solution
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: f,g,h    !> Physical fluxes in x/y/z direction
REAL,DIMENSION(3          ),INTENT(IN)  :: Mf,Mg,Mh !> Metrics in x/y/z direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Ep
REAL                :: Mmom
!==================================================================================================================================
!! auxiliary variables
!Ep = (U(ENER) + UPrim(PRES))/U(DENS)
!
!f(DENS) = Mf(1)* U(MOM1) &
!        + Mf(2)* U(MOM2) &
!        + Mf(3)* U(MOM3)
!f(MOM1) = Mf(1)*(U(MOM1) * UPrim(VEL1) + UPrim(PRES)) &
!        + Mf(2)*(U(MOM1) * UPrim(VEL2)              ) &
!        + Mf(3)*(U(MOM1) * UPrim(VEL3)              )
!f(MOM2) = Mf(1)*(U(MOM2) * UPrim(VEL1)              ) &
!        + Mf(2)*(U(MOM2) * UPrim(VEL2) + UPrim(PRES)) &
!        + Mf(3)*(U(MOM2) * UPrim(VEL3)              )
!f(MOM3) = Mf(1)*(U(MOM3) * UPrim(VEL1)              ) &
!        + Mf(2)*(U(MOM3) * UPrim(VEL2)              ) &
!        + Mf(3)*(U(MOM3) * UPrim(VEL3) + UPrim(PRES))
!f(ENER) = Mf(1)*(Ep*UPrim(VEL1)) &
!        + Mf(2)*(Ep*UPrim(VEL2)) &
!        + Mf(3)*(Ep*UPrim(VEL3))
!
!g(DENS) = Mg(1)* U(MOM1) &
!        + Mg(2)* U(MOM2) &
!        + Mg(3)* U(MOM3)
!g(MOM1) = Mg(1)*(U(MOM1) * UPrim(VEL1) + UPrim(PRES)) &
!        + Mg(2)*(U(MOM1) * UPrim(VEL2)              ) &
!        + Mg(3)*(U(MOM1) * UPrim(VEL3)              )
!g(MOM2) = Mg(1)*(U(MOM2) * UPrim(VEL1)              ) &
!        + Mg(2)*(U(MOM2) * UPrim(VEL2) + UPrim(PRES)) &
!        + Mg(3)*(U(MOM2) * UPrim(VEL3)              )
!g(MOM3) = Mg(1)*(U(MOM3) * UPrim(VEL1)              ) &
!        + Mg(2)*(U(MOM3) * UPrim(VEL2)              ) &
!        + Mg(3)*(U(MOM3) * UPrim(VEL3) + UPrim(PRES))
!g(ENER) = Mg(1)*(Ep*UPrim(VEL1)) &
!        + Mg(2)*(Ep*UPrim(VEL2)) &
!        + Mg(3)*(Ep*UPrim(VEL3))
!
!h(DENS) = Mh(1)* U(MOM1) &
!        + Mh(2)* U(MOM2) &
!        + Mh(3)* U(MOM3)
!h(MOM1) = Mh(1)*(U(MOM1) * UPrim(VEL1) + UPrim(PRES)) &
!        + Mh(2)*(U(MOM1) * UPrim(VEL2)              ) &
!        + Mh(3)*(U(MOM1) * UPrim(VEL3)              )
!h(MOM2) = Mh(1)*(U(MOM2) * UPrim(VEL1)              ) &
!        + Mh(2)*(U(MOM2) * UPrim(VEL2) + UPrim(PRES)) &
!        + Mh(3)*(U(MOM2) * UPrim(VEL3)              )
!h(MOM3) = Mh(1)*(U(MOM3) * UPrim(VEL1)              ) &
!        + Mh(2)*(U(MOM3) * UPrim(VEL2)              ) &
!        + Mh(3)*(U(MOM3) * UPrim(VEL3) + UPrim(PRES))
!h(ENER) = Mh(1)*(Ep*UPrim(VEL1)) &
!        + Mh(2)*(Ep*UPrim(VEL2)) &
!        + Mh(3)*(Ep*UPrim(VEL3))

! Heavily optimized variant. But seems to be not much faster...
! auxiliary variables
Ep = (U(ENER) + UPrim(PRES))/U(DENS)

Mmom    = DOT_PRODUCT(Mf(:),U(MOMV)) ! metrics times momentum vector
f(DENS) =  Mmom
f(MOM1) =  Mmom*UPrim(VEL1) &
        + Mf(1)*UPrim(PRES)
f(MOM2) =  Mmom*UPrim(VEL2) &
        + Mf(2)*UPrim(PRES)
f(MOM3) =  Mmom*UPrim(VEL3) &
        + Mf(3)*UPrim(PRES)
f(ENER) = Mmom*Ep

Mmom    = DOT_PRODUCT(Mg(:),U(MOMV)) ! metrics times momentum vector
g(DENS) =  Mmom
g(MOM1) =  Mmom*UPrim(VEL1) &
        + Mg(1)*UPrim(PRES)
g(MOM2) =  Mmom*UPrim(VEL2) &
        + Mg(2)*UPrim(PRES)
g(MOM3) =  Mmom*UPrim(VEL3) &
        + Mg(3)*UPrim(PRES)
g(ENER) =  Mmom*Ep

Mmom    = DOT_PRODUCT(Mh(:),U(MOMV)) ! metrics times momentum vector
h(DENS) =  Mmom
h(MOM1) =  Mmom*UPrim(VEL1) &
        + Mh(1)*UPrim(PRES)
h(MOM2) =  Mmom*UPrim(VEL2) &
        + Mh(2)*UPrim(PRES)
h(MOM3) =  Mmom*UPrim(VEL3) &
        + Mh(3)*UPrim(PRES)
h(ENER) =  Mmom*Ep

END SUBROUTINE EvalTransformedFlux3D

!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
PPURE SUBROUTINE EvalFlux3D_Elem(Nloc,U,UPrim,f,g,h)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER                                               ,INTENT(IN)  :: Nloc     !< Polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U        !< Conservative solution
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h    !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,ZDIM(Nloc);  DO j=0,Nloc; DO i=0,Nloc
  CALL EvalFlux3D(U(:,i,j,k),UPrim(:,i,j,k),f(:,i,j,k),g(:,i,j,k),h(:,i,j,k))
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalFlux3D_Elem

!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE EvalFlux3D_CUDA_Kernel(nDOF,U,UPrim,f,g,h)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE                            ,INTENT(IN)  :: nDOF   !< number of degrees of freedom in arrays
REAL,DEVICE,DIMENSION(PP_nVar    ,1:nDOF),INTENT(IN)  :: U      !< Conservative solution
REAL,DEVICE,DIMENSION(PP_nVarPrim,1:nDOF),INTENT(IN)  :: UPrim  !< Primitive solution
REAL,DEVICE,DIMENSION(PP_nVar    ,1:nDOF),INTENT(OUT) :: f,g,h  !> Physical fluxes in x,y,z
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
!==================================================================================================================================
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF (i.LE.nDOF) CALL EvalFlux3D(U(:,i),UPrim(:,i),f(:,i),g(:,i),h(:,i))
END SUBROUTINE EvalFlux3D_CUDA_Kernel

!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE EvalTransformedFlux3D_Kernel(nDOF,U,UPrim,f,g,h,Mf,Mg,Mh)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE                            ,INTENT(IN)  :: nDOF   !< number of degrees of freedom in arrays
REAL,DEVICE,DIMENSION(PP_nVar    ,1:nDOF),INTENT(IN)  :: U      !< Conservative solution
REAL,DEVICE,DIMENSION(PP_nVarPrim,1:nDOF),INTENT(IN)  :: UPrim  !< Primitive solution
REAL,DEVICE,DIMENSION(PP_nVar    ,1:nDOF),INTENT(OUT) :: f,g,h  !> Physical fluxes in x,y,z
REAL,DEVICE,DIMENSION(3          ,1:nDOF),INTENT(IN)  :: Mf,Mg,Mh  !< Primitive solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
!==================================================================================================================================
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF (i.LE.nDOF) CALL EvalTransformedFlux3D(U(:,i),UPrim(:,i),f(:,i),g(:,i),h(:,i),Mf(:,i),Mg(:,i),Mh(:,i))
END SUBROUTINE EvalTransformedFlux3D_Kernel

!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
PPURE SUBROUTINE EvalFlux3D_Elem_CUDA(Nloc,U,UPrim,f,g,h)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE                                                ,INTENT(IN)  :: Nloc     !< Polynomial degree
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U        !< Conservative solution
REAL,DEVICE,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h    !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: nDOF
INTEGER,PARAMETER   :: nThreads=256
!==================================================================================================================================
nDOF=(Nloc+1)*(Nloc+1)*(ZDIM(Nloc)+1)
CALL EvalFlux3D_CUDA_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,U,UPrim,f,g,h)
END SUBROUTINE EvalFlux3D_Elem_CUDA

!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
PPURE SUBROUTINE EvalFlux3D_Elems_CUDA(Nloc,nElems,U,UPrim,f,g,h)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE                                                       ,INTENT(IN)  :: Nloc     !< Polynomial degree
INTEGER,VALUE                                                       ,INTENT(IN)  :: nElems   !< Polynomial degree
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems),INTENT(IN)  :: U        !< Conservative solution
REAL,DEVICE,DIMENSION(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems),INTENT(OUT) :: f,g,h    !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: nDOF
INTEGER,PARAMETER   :: nThreads=256
!==================================================================================================================================
nDOF=(Nloc+1)*(Nloc+1)*(ZDIM(Nloc)+1)*nElems
CALL EvalFlux3D_CUDA_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,U,UPrim,f,g,h)
END SUBROUTINE EvalFlux3D_Elems_CUDA

#if PARABOLIC
!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux using the primitive variables and derivatives.
!==================================================================================================================================
PPURE SUBROUTINE ATTRIBUTES(DEVICE,HOST) EvalDiffFlux3D(UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim                 !< Solution vector
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx,gradUy,gradUz  !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: f,g,h                 !> Physical fluxes in x,y,z directions
REAL,INTENT(IN)  :: mu  !< viscosity of fluid
REAL,INTENT(IN)  :: mu  !< thermal conductivity of fluid
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: tau_zz,tau_xz,tau_yz
#endif
REAL,PARAMETER      :: s23=2./3.
REAL,PARAMETER      :: s43=4./3.
!==================================================================================================================================
!! ideal gas law
!muS    = VISCOSITY_PRIM(UPrim)
!lambda = THERMAL_CONDUCTIVITY_H(muS)
!#if EDDYVISCOSITY
!!Add turbulent sub grid scale viscosity to mu
!muS    = muS    + muSGS
!lambda = lambda + muSGS*cp/PrSGS
!#endif

#if PP_dim==3
! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

! viscous fluxes in x-direction
tau_xx = muS * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = muS * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_zz = muS * (-s23 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) + s43 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
tau_xy = muS * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))  ! mu*(u_y+v_x)
tau_xz = muS * (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))  ! mu*(u_z+w_x)
tau_yz = muS * (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))  ! mu*(y_z+w_y)

f(DENS) = 0.
f(MOM1) = -tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(MOM2) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
f(MOM3) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
f(ENER) = -tau_xx*UPrim(VEL1)-tau_xy*UPrim(VEL2)-tau_xz*UPrim(VEL3) &
          -lambda*gradUx(LIFT_TEMP)                     ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
g(MOM2) = -tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(MOM3) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
g(ENER) = -tau_xy*UPrim(VEL1)-tau_yy*UPrim(VEL2)-tau_yz*UPrim(VEL3) &
          -lambda*gradUy(LIFT_TEMP)                     ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
! viscous fluxes in z-direction
h(DENS) = 0.
h(MOM1) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
h(MOM2) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
h(MOM3) = -tau_zz                                       ! F_euler-4/3*mu*w_z+2/3*mu*(u_x+v_y)
h(ENER) = -tau_xz*UPrim(VEL1)-tau_yz*UPrim(VEL2)-tau_zz*UPrim(VEL3) &
          -lambda*gradT3                                ! F_euler-(tau_zx*u+tau_zy*v+tau_zz*w-q_z) q_z=-lambda*T_z
#else
! gradients of primitive variables are directly available gradU = (/ drho, dv1, dv2, dv3, dT /)

! viscous fluxes in x-direction
tau_xx = muS * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2))  ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = muS * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2))  !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_xy = muS * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))               ! mu*(u_y+v_x)

f(DENS) = 0.
f(MOM1) = -tau_xx                                                         ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(MOM2) = -tau_xy                                                         ! F_euler-mu*(u_y+v_x)
f(MOM3) = 0.
f(ENER) = -tau_xx*UPrim(VEL1)-tau_xy*UPrim(VEL2)-lambda*gradUx(LIFT_TEMP) ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) q_x=-lambda*T_x
! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -tau_xy                                                         ! F_euler-mu*(u_y+v_x)
g(MOM2) = -tau_yy                                                         ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(MOM3) = 0.
g(ENER) = -tau_xy*UPrim(VEL1)-tau_yy*UPrim(VEL2)-lambda*gradUy(LIFT_TEMP) ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) q_y=-lambda*T_y
! viscous fluxes in z-direction
h    = 0.
#endif
END SUBROUTINE EvalDiffFlux3D

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Volume(UPrim,gradUx,gradUy,gradUz,f,g,h,iElem)
! MODULES
USE MOD_PreProc
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)   :: UPrim                !< Solution vector
REAL,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)   :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT)  :: f,g,h                !> Physical fluxes in x,y,z directions
INTEGER,INTENT(IN)      :: iElem                !< element index in global array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,PP_NZ;  DO j=0,PP_N; DO i=0,PP_N
  CALL EvalDiffFlux3D_Point(Uprim(:,i,j,k),gradUx(:,i,j,k),gradUy(:,i,j,k),gradUz(:,i,j,k), &
                                                f(:,i,j,k),     g(:,i,j,k),     h(:,i,j,k)  &
                           )
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalDiffFlux3D_Volume

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single side
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Surface(Nloc,UPrim,gradUx,gradUy,gradUz,f,g,h &
                                 )
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: Nloc  !< Polynomial degree of input solution
REAL,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h                !> Physical fluxes in x,y,z directions
#if EDDYVISCOSITY
REAL,DIMENSION(1             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: muSGS                !< SGS viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j
!==================================================================================================================================
DO j=0,ZDIM(Nloc); DO i=0,PP_N
  CALL EvalDiffFlux3D_Point(Uprim(:,i,j),gradUx(:,i,j),gradUy(:,i,j),gradUz(:,i,j), &
                                              f(:,i,j),     g(:,i,j),     h(:,i,j)  &
                           )
END DO; END DO ! i,j
END SUBROUTINE EvalDiffFlux3D_Surface
#endif /*PARABOLIC*/

!==================================================================================================================================
!> Computes 1D Euler flux using the conservative variables.
!==================================================================================================================================
PPURE SUBROUTINE EvalEulerFlux1D(U,F)
! MODULES
USE MOD_EOS_Vars ,ONLY:KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_nVar)  !< vector of conservative variables
REAL,INTENT(OUT)    :: F(PP_nVar)  !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: UE(PP_2Var)  ! auxiliary variables
!==================================================================================================================================
! auxiliary variables
! TODO: ATTENTION: Temperature of UE not filled!!!
UE(EXT_CONS)=U
UE(EXT_SRHO)=1./UE(EXT_DENS)
UE(EXT_VELV)=VELOCITY_HE(UE)
UE(EXT_PRES)=PRESSURE_HE(UE)
! Euler fluxes x-direction
F(DENS)= U(MOM1)                             ! rho*u
F(MOM1)= U(MOM1)*UE(EXT_VEL1)+UE(EXT_PRES)   ! rho*u²+p
F(MOM2)= U(MOM1)*UE(EXT_VEL2)                ! rho*u*v
#if PP_dim==3
F(MOM3)=U(MOM1)*UE(EXT_VEL3)                 ! rho*u*w
#else
F(MOM3)=0.
#endif
F(ENER)=(U(ENER)+UE(EXT_PRES))*UE(EXT_VEL1)  ! (rho*e+p)*u
END SUBROUTINE EvalEulerFlux1D

!==================================================================================================================================
!> Computes 1D Euler flux using the conservative and primitive variables (for better performance)
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalEulerFlux1D_fast(U,F)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_2Var)  !< vector of conservative and primitive variables
REAL,INTENT(OUT)    :: F(PP_nVar)  !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Euler fluxes x-direction
F(DENS)= U(EXT_MOM1)                         ! rho*u
F(MOM1)= U(EXT_MOM1)*U(EXT_VEL1)+U(EXT_PRES) ! rho*u²+p
F(MOM2)= U(EXT_MOM1)*U(EXT_VEL2)             ! rho*u*v
#if PP_dim==3
F(MOM3)= U(EXT_MOM1)*U(EXT_VEL3)             ! rho*u*w
#else
F(MOM3)= 0.
#endif
F(ENER)=(U(EXT_ENER)+U(EXT_PRES))*U(EXT_VEL1)! (rho*e+p)*u
END SUBROUTINE EvalEulerFlux1D_fast

END MODULE MOD_Flux
