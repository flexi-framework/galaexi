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
  MODULE PROCEDURE EvalNormalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D_Point
  MODULE PROCEDURE EvalDiffFlux3D_Surface
  MODULE PROCEDURE EvalDiffFlux3D_Elem
  MODULE PROCEDURE EvalDiffFlux3D_Surface_CUDA
  MODULE PROCEDURE EvalDiffFlux3D_Sides_CUDA
  MODULE PROCEDURE EvalDiffFlux3D_Elem_CUDA
  MODULE PROCEDURE EvalDiffFlux3D_Elems_CUDA
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalFlux3D, EvalEulerFlux1D, EvalEulerFlux1D_fast,EvalTransformedFlux3D
#if PARABOLIC
PUBLIC::EvalDiffFlux3D
PUBLIC::EvalTransformedDiffFlux3D_CUDA
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
!> Heavily optimized variant.
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalTransformedEulerFlux3D_fast(U,UPrim,f,g,h,Mf,Mg,Mh)
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
! auxiliary variables
Ep = (U(ENER) + UPrim(PRES))/U(DENS)

! fluxes in x-direction
Mmom    = DOT_PRODUCT(Mf(:),U(MOMV)) ! metrics times momentum vector
f(DENS) =  Mmom
f(MOM1) =  Mmom*UPrim(VEL1) &
        + Mf(1)*UPrim(PRES)
f(MOM2) =  Mmom*UPrim(VEL2) &
        + Mf(2)*UPrim(PRES)
f(MOM3) =  Mmom*UPrim(VEL3) &
        + Mf(3)*UPrim(PRES)
f(ENER) = Mmom*Ep

! fluxes in y-direction
Mmom    = DOT_PRODUCT(Mg(:),U(MOMV)) ! metrics times momentum vector
g(DENS) =  Mmom
g(MOM1) =  Mmom*UPrim(VEL1) &
        + Mg(1)*UPrim(PRES)
g(MOM2) =  Mmom*UPrim(VEL2) &
        + Mg(2)*UPrim(PRES)
g(MOM3) =  Mmom*UPrim(VEL3) &
        + Mg(3)*UPrim(PRES)
g(ENER) =  Mmom*Ep

! fluxes in z-direction
Mmom    = DOT_PRODUCT(Mh(:),U(MOMV)) ! metrics times momentum vector
h(DENS) =  Mmom
h(MOM1) =  Mmom*UPrim(VEL1) &
        + Mh(1)*UPrim(PRES)
h(MOM2) =  Mmom*UPrim(VEL2) &
        + Mh(2)*UPrim(PRES)
h(MOM3) =  Mmom*UPrim(VEL3) &
        + Mh(3)*UPrim(PRES)
h(ENER) =  Mmom*Ep
END SUBROUTINE EvalTransformedEulerFlux3D_fast

!==================================================================================================================================
!> Compute advection part of the Navier-Stokes fluxes in all space dimensions using the conservative and primitive variables
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalTransformedEulerFlux3D(U,UPrim,f,g,h,Mf,Mg,Mh)
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
! auxiliary variables
Ep = (U(ENER) + UPrim(PRES))

! fluxes in x-direction
f(DENS) = Mf(1)* U(MOM1) &
        + Mf(2)* U(MOM2) &
        + Mf(3)* U(MOM3)
f(MOM1) = Mf(1)*(U(MOM1) * UPrim(VEL1) + UPrim(PRES)) &
        + Mf(2)*(U(MOM1) * UPrim(VEL2)              ) &
        + Mf(3)*(U(MOM1) * UPrim(VEL3)              )
f(MOM2) = Mf(1)*(U(MOM2) * UPrim(VEL1)              ) &
        + Mf(2)*(U(MOM2) * UPrim(VEL2) + UPrim(PRES)) &
        + Mf(3)*(U(MOM2) * UPrim(VEL3)              )
f(MOM3) = Mf(1)*(U(MOM3) * UPrim(VEL1)              ) &
        + Mf(2)*(U(MOM3) * UPrim(VEL2)              ) &
        + Mf(3)*(U(MOM3) * UPrim(VEL3) + UPrim(PRES))
f(ENER) = Mf(1)*(Ep*UPrim(VEL1)) &
        + Mf(2)*(Ep*UPrim(VEL2)) &
        + Mf(3)*(Ep*UPrim(VEL3))

! fluxes in y-direction
g(DENS) = Mg(1)* U(MOM1) &
        + Mg(2)* U(MOM2) &
        + Mg(3)* U(MOM3)
g(MOM1) = Mg(1)*(U(MOM1) * UPrim(VEL1) + UPrim(PRES)) &
        + Mg(2)*(U(MOM1) * UPrim(VEL2)              ) &
        + Mg(3)*(U(MOM1) * UPrim(VEL3)              )
g(MOM2) = Mg(1)*(U(MOM2) * UPrim(VEL1)              ) &
        + Mg(2)*(U(MOM2) * UPrim(VEL2) + UPrim(PRES)) &
        + Mg(3)*(U(MOM2) * UPrim(VEL3)              )
g(MOM3) = Mg(1)*(U(MOM3) * UPrim(VEL1)              ) &
        + Mg(2)*(U(MOM3) * UPrim(VEL2)              ) &
        + Mg(3)*(U(MOM3) * UPrim(VEL3) + UPrim(PRES))
g(ENER) = Mg(1)*(Ep*UPrim(VEL1)) &
        + Mg(2)*(Ep*UPrim(VEL2)) &
        + Mg(3)*(Ep*UPrim(VEL3))

! fluxes in z-direction
h(DENS) = Mh(1)* U(MOM1) &
        + Mh(2)* U(MOM2) &
        + Mh(3)* U(MOM3)
h(MOM1) = Mh(1)*(U(MOM1) * UPrim(VEL1) + UPrim(PRES)) &
        + Mh(2)*(U(MOM1) * UPrim(VEL2)              ) &
        + Mh(3)*(U(MOM1) * UPrim(VEL3)              )
h(MOM2) = Mh(1)*(U(MOM2) * UPrim(VEL1)              ) &
        + Mh(2)*(U(MOM2) * UPrim(VEL2) + UPrim(PRES)) &
        + Mh(3)*(U(MOM2) * UPrim(VEL3)              )
h(MOM3) = Mh(1)*(U(MOM3) * UPrim(VEL1)              ) &
        + Mh(2)*(U(MOM3) * UPrim(VEL2)              ) &
        + Mh(3)*(U(MOM3) * UPrim(VEL3) + UPrim(PRES))
h(ENER) = Mh(1)*(Ep*UPrim(VEL1)) &
        + Mh(2)*(Ep*UPrim(VEL2)) &
        + Mh(3)*(Ep*UPrim(VEL3))
END SUBROUTINE EvalTransformedEulerFlux3D

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
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE EvalTransformedFlux3D_Kernel(nDOF,U,UPrim &
#if PARABOLIC
                                                                ,gradUx,gradUy,gradUz,mu,lambda &
#endif
                                                                ,f,g,h,Mf,Mg,Mh)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)  :: nDOF     !< number of degrees of freedom in arrays
REAL,DEVICE,DIMENSION(PP_nVar       ,1:nDOF),INTENT(IN)  :: U        !< Conservative solution
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,1:nDOF),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DEVICE,DIMENSION(3             ,1:nDOF),INTENT(IN)  :: Mf,Mg,Mh !< Metrics in x,y,z
REAL,DEVICE,DIMENSION(PP_nVar       ,1:nDOF),INTENT(OUT) :: f,g,h    !> Physical fluxes in x,y,z
#if PARABOLIC
REAL,DEVICE,DIMENSION(PP_nVarLifting,1:nDOF),INTENT(IN)  :: gradUx,gradUy,gradUz !> gradients in x,y,z
REAL,VALUE,INTENT(IN) :: mu,lambda
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
!==================================================================================================================================
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF (i.LE.nDOF) THEN
#if PARABOLIC
  CALL EvalTransformedEulerDiffFlux3D(U(:,i),UPrim(:,i),gradUx(:,i),gradUy(:,i),gradUz(:,i) &
                                     ,f(:,i),g(:,i),h(:,i),Mf(:,i),Mg(:,i),Mh(:,i),mu,lambda)
#else
  CALL EvalTransformedEulerFlux3D_fast(U(:,i),UPrim(:,i),f(:,i),g(:,i),h(:,i),Mf(:,i),Mg(:,i),Mh(:,i))
#endif
END IF
END SUBROUTINE EvalTransformedFlux3D_Kernel

#if PARABOLIC
!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE EvalTransformedDiffFlux3D_Kernel(nDOF,UPrim &
                                                           ,gradUx,gradUy,gradUz,mu,lambda &
                                                           ,f,g,h,Mf,Mg,Mh)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)  :: nDOF     !< number of degrees of freedom in arrays
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,1:nDOF),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DEVICE,DIMENSION(3             ,1:nDOF),INTENT(IN)  :: Mf,Mg,Mh !< Metrics in x,y,z
REAL,DEVICE,DIMENSION(PP_nVar       ,1:nDOF),INTENT(OUT) :: f,g,h    !> Physical fluxes in x,y,z
REAL,DEVICE,DIMENSION(PP_nVarLifting,1:nDOF),INTENT(IN)  :: gradUx,gradUy,gradUz !> gradients in x,y,z
REAL,VALUE,INTENT(IN) :: mu,lambda
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
!==================================================================================================================================
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF (i.LE.nDOF) THEN
  CALL EvalTransformedDiffFlux3D( UPrim(:,i) &
                                ,gradUx(:,i),gradUy(:,i),gradUz(:,i) &
                                ,     f(:,i),     g(:,i),     h(:,i) &
                                ,    Mf(:,i),    Mg(:,i),    Mh(:,i) &
                                ,mu,lambda)
END IF
END SUBROUTINE EvalTransformedDiffFlux3D_Kernel
#endif

!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes
!==================================================================================================================================
PPURE SUBROUTINE EvalTransformedFlux3D(nDOF,U,UPrim &
#if PARABOLIC
                                      ,gradUx,gradUy,gradUz &
#endif
                                      ,f,g,h,Mf,Mg,Mh)
! MODULES
#if PARABOLIC
USE MOD_EOS_Vars,ONLY: mu0,cp,Pr
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)  :: nDOF     !< number of degrees of freedom in arrays
REAL,DEVICE,DIMENSION(PP_nVar       ,1:nDOF),INTENT(IN)  :: U        !< Conservative solution
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,1:nDOF),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DEVICE,DIMENSION(3             ,1:nDOF),INTENT(IN)  :: Mf,Mg,Mh !< Metrics in x,y,z
REAL,DEVICE,DIMENSION(PP_nVar       ,1:nDOF),INTENT(OUT) :: f,g,h    !> Physical fluxes in x,y,z
#if PARABOLIC
REAL,DEVICE,DIMENSION(PP_nVarLifting,1:nDOF),INTENT(IN)  :: gradUx,gradUy,gradUz !> gradients in x,y,z
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: nThreads=256
#if PARABOLIC
REAL                :: mu,lambda
#endif
!==================================================================================================================================
! TODO: Impelment cleanly!
#if PARABOLIC
mu     = VISCOSITY_PRIM(UPrim(:,1))
lambda = THERMAL_CONDUCTIVITY_H(mu)
#endif
CALL EvalTransformedFlux3D_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,U,UPrim &
#if PARABOLIC
                                                               ,gradUx,gradUy,gradUz,mu,lambda &
#endif
                                                               ,f,g,h,Mf,Mg,Mh)
END SUBROUTINE EvalTransformedFlux3D

#if PARABOLIC
!==================================================================================================================================
!> Wrapper routine to compute the advection part of the Navier-Stokes fluxes
!==================================================================================================================================
PPURE SUBROUTINE EvalTransformedDiffFlux3D_CUDA(nDOF,UPrim &
                                          ,gradUx,gradUy,gradUz &
                                          ,f,g,h,Mf,Mg,Mh)
! MODULES
USE MOD_EOS_Vars,ONLY: mu0,cp,Pr
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)  :: nDOF     !< number of degrees of freedom in arrays
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,1:nDOF),INTENT(IN)  :: UPrim    !< Primitive solution
REAL,DEVICE,DIMENSION(3             ,1:nDOF),INTENT(IN)  :: Mf,Mg,Mh !< Metrics in x,y,z
REAL,DEVICE,DIMENSION(PP_nVar       ,1:nDOF),INTENT(OUT) :: f,g,h    !> Physical fluxes in x,y,z
REAL,DEVICE,DIMENSION(PP_nVarLifting,1:nDOF),INTENT(IN)  :: gradUx,gradUy,gradUz !> gradients in x,y,z
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER   :: nThreads=256
#if PARABOLIC
REAL                :: mu,lambda
#endif
!==================================================================================================================================
! TODO: Impelment cleanly!
mu     = VISCOSITY_PRIM(UPrim(:,1))
lambda = THERMAL_CONDUCTIVITY_H(mu)
CALL EvalTransformedDiffFlux3D_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,UPrim &
                                                                   ,gradUx,gradUy,gradUz,mu,lambda &
                                                                   ,f,g,h,Mf,Mg,Mh)
END SUBROUTINE EvalTransformedDiffFlux3D_CUDA
#endif

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
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalTransformedEulerDiffFlux3D(U,UPrim,gradUx,gradUy,gradUz,f,g,h,Mf,Mg,Mh,mu,lambda)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar       ),INTENT(IN)  :: U                     !< Solution vector in conservative variables
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim                 !< Solution vector in primitive variables
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx,gradUy,gradUz  !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: f,g,h                 !> Physical fluxes in x,y,z directions
REAL,DIMENSION(3             ),INTENT(IN)  :: Mf,Mg,Mh              !> Metrics in x/y/z direction
REAL,INTENT(IN)  :: mu      !< viscosity of fluid
REAL,INTENT(IN)  :: lambda  !< thermal conductivity of fluid
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Ep,Mmom
REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: tau_zz,tau_xz,tau_yz
#endif
REAL                :: tauXVelPlusQ,tauYVelPlusQ,tauZVelPlusQ
REAL,PARAMETER      :: s23=2./3.
REAL,PARAMETER      :: s43=4./3.
!==================================================================================================================================
! auxiliary variables
Ep = (U(ENER) + UPrim(PRES))/U(DENS)
! Precompute entries of shear-stress tensor
tau_xx = mu * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = mu * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_zz = mu * (-s23 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) + s43 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
tau_xy = mu * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))  ! mu*(u_y+v_x)
tau_xz = mu * (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))  ! mu*(u_z+w_x)
tau_yz = mu * (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))  ! mu*(y_z+w_y)
! Precompute tau*Vel+q (shear stress tensor times velocity vector plus heat flux vector)
tauXVelPlusQ = tau_xx*UPrim(VEL1)+tau_xy*UPrim(VEL2)+tau_xz*UPrim(VEL3)+lambda*gradUx(LIFT_TEMP)
tauYVelPlusQ = tau_xy*UPrim(VEL1)+tau_yy*UPrim(VEL2)+tau_yz*UPrim(VEL3)+lambda*gradUy(LIFT_TEMP)
tauZVelPlusQ = tau_xz*UPrim(VEL1)+tau_yz*UPrim(VEL2)+tau_zz*UPrim(VEL3)+lambda*gradUz(LIFT_TEMP)

! viscous fluxes in x-direction
Mmom    = DOT_PRODUCT(Mf(:),U(MOMV)) ! metrics times momentum vector
f(DENS) =  Mmom
f(MOM1) =  Mmom*UPrim(VEL1)  &
        + Mf(1)*UPrim(PRES)  &
        - Mf(1)*tau_xx       &
        - Mf(2)*tau_xy       &
        - Mf(3)*tau_xz
f(MOM2) =  Mmom*UPrim(VEL2)  &
        + Mf(2)*UPrim(PRES)  &
        - Mf(1)*tau_xy       &
        - Mf(2)*tau_yy       &
        - Mf(3)*tau_yz
f(MOM3) =  Mmom*UPrim(VEL3)  &
        + Mf(3)*UPrim(PRES)  &
        - Mf(1)*tau_xz       &
        - Mf(2)*tau_yz       &
        - Mf(3)*tau_zz
f(ENER) =  Mmom*Ep           &
        - Mf(1)*tauXVelPlusQ &
        - Mf(2)*tauYVelPlusQ &
        - Mf(3)*tauZVelPlusQ

! viscous fluxes in y-direction
Mmom    = DOT_PRODUCT(Mg(:),U(MOMV)) ! metrics times momentum vector
g(DENS) =  Mmom
g(MOM1) =  Mmom*UPrim(VEL1)  &
        + Mg(1)*UPrim(PRES)  &
        - Mg(1)*tau_xx       &
        - Mg(2)*tau_xy       &
        - Mg(3)*tau_xz
g(MOM2) =  Mmom*UPrim(VEL2)  &
        + Mg(2)*UPrim(PRES)  &
        - Mg(1)*tau_xy       &
        - Mg(2)*tau_yy       &
        - Mg(3)*tau_yz
g(MOM3) =  Mmom*UPrim(VEL3)  &
        + Mg(3)*UPrim(PRES)  &
        - Mg(1)*tau_xz       &
        - Mg(2)*tau_yz       &
        - Mg(3)*tau_zz
g(ENER) =  Mmom*Ep           &
        - Mg(1)*tauXVelPlusQ &
        - Mg(2)*tauYVelPlusQ &
        - Mg(3)*tauZVelPlusQ

! viscous fluxes in z-direction
Mmom    = DOT_PRODUCT(Mh(:),U(MOMV)) ! metrics times momentum vector
h(DENS) =  Mmom
h(MOM1) =  Mmom*UPrim(VEL1)  &
        + Mh(1)*UPrim(PRES)  &
        - Mh(1)*tau_xx       &
        - Mh(2)*tau_xy       &
        - Mh(3)*tau_xz
h(MOM2) =  Mmom*UPrim(VEL2)  &
        + Mh(2)*UPrim(PRES)  &
        - Mh(1)*tau_xy       &
        - Mh(2)*tau_yy       &
        - Mh(3)*tau_yz
h(MOM3) =  Mmom*UPrim(VEL3)  &
        + Mh(3)*UPrim(PRES)  &
        - Mh(1)*tau_xz       &
        - Mh(2)*tau_yz       &
        - Mh(3)*tau_zz
h(ENER) =  Mmom*Ep           &
        - Mh(1)*tauXVelPlusQ &
        - Mh(2)*tauYVelPlusQ &
        - Mh(3)*tauZVelPlusQ
END SUBROUTINE EvalTransformedEulerDiffFlux3D

!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux using the primitive variables and derivatives.
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalTransformedDiffFlux3D(UPrim,gradUx,gradUy,gradUz,f,g,h,Mf,Mg,Mh,mu,lambda)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim                 !< Solution vector
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx,gradUy,gradUz  !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: f,g,h                 !> Physical fluxes in x,y,z directions
REAL,DIMENSION(3             ),INTENT(IN)  :: Mf,Mg,Mh              !> Metrics in x/y/z direction
REAL,INTENT(IN)  :: mu      !< viscosity of fluid
REAL,INTENT(IN)  :: lambda  !< thermal conductivity of fluid
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: tau_zz,tau_xz,tau_yz
#endif
REAL                :: tauXVelPlusQ,tauYVelPlusQ,tauZVelPlusQ
REAL,PARAMETER      :: s23=2./3.
REAL,PARAMETER      :: s43=4./3.
!==================================================================================================================================
! Precompute entries of shear-stress tensor
tau_xx = mu * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = mu * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_zz = mu * (-s23 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) + s43 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
tau_xy = mu * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))  ! mu*(u_y+v_x)
tau_xz = mu * (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))  ! mu*(u_z+w_x)
tau_yz = mu * (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))  ! mu*(y_z+w_y)

! Precompute tau*Vel+q (shear stress tensor times velocity vector plus heat flux vector)
tauXVelPlusQ = tau_xx*UPrim(VEL1)+tau_xy*UPrim(VEL2)+tau_xz*UPrim(VEL3)+lambda*gradUx(LIFT_TEMP)
tauYVelPlusQ = tau_xy*UPrim(VEL1)+tau_yy*UPrim(VEL2)+tau_yz*UPrim(VEL3)+lambda*gradUy(LIFT_TEMP)
tauZVelPlusQ = tau_xz*UPrim(VEL1)+tau_yz*UPrim(VEL2)+tau_zz*UPrim(VEL3)+lambda*gradUz(LIFT_TEMP)

! viscous fluxes in x-direction
f(DENS) = 0.
f(MOM1) = -Mf(1)*tau_xx &
          -Mf(2)*tau_xy &
          -Mf(3)*tau_xz
f(MOM2) = -Mf(1)*tau_xy &
          -Mf(2)*tau_yy &
          -Mf(3)*tau_yz
f(MOM3) = -Mf(1)*tau_xz &
          -Mf(2)*tau_yz &
          -Mf(3)*tau_zz
f(ENER) = -Mf(1)*tauXVelPlusQ &
          -Mf(2)*tauYVelPlusQ &
          -Mf(3)*tauZVelPlusQ

! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -Mg(1)*tau_xx &
          -Mg(2)*tau_xy &
          -Mg(3)*tau_xz
g(MOM2) = -Mg(1)*tau_xy &
          -Mg(2)*tau_yy &
          -Mg(3)*tau_yz
g(MOM3) = -Mg(1)*tau_xz &
          -Mg(2)*tau_yz &
          -Mg(3)*tau_zz
g(ENER) = -Mg(1)*tauXVelPlusQ &
          -Mg(2)*tauYVelPlusQ &
          -Mg(3)*tauZVelPlusQ

! viscous fluxes in z-direction
h(DENS) = 0.
h(MOM1) = -Mh(1)*tau_xx &
          -Mh(2)*tau_xy &
          -Mh(3)*tau_xz
h(MOM2) = -Mh(1)*tau_xy &
          -Mh(2)*tau_yy &
          -Mh(3)*tau_yz
h(MOM3) = -Mh(1)*tau_xz &
          -Mh(2)*tau_yz &
          -Mh(3)*tau_zz
h(ENER) = -Mh(1)*tauXVelPlusQ &
          -Mh(2)*tauYVelPlusQ &
          -Mh(3)*tauZVelPlusQ
END SUBROUTINE EvalTransformedDiffFlux3D

!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux using the primitive variables and derivatives.
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalDiffFlux3D(UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim                 !< Solution vector
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx,gradUy,gradUz  !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: f,g,h                 !> Physical fluxes in x,y,z directions
REAL,INTENT(IN)  :: mu      !< viscosity of fluid
REAL,INTENT(IN)  :: lambda  !< thermal conductivity of fluid
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: tau_zz,tau_xz,tau_yz
#endif
REAL,PARAMETER      :: s23=2./3.
REAL,PARAMETER      :: s43=4./3.
!==================================================================================================================================
#if PP_dim==3
! Precompute entries of shear-stress tensor
tau_xx = mu * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = mu * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_zz = mu * (-s23 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) + s43 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
tau_xy = mu * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))  ! mu*(u_y+v_x)
tau_xz = mu * (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))  ! mu*(u_z+w_x)
tau_yz = mu * (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))  ! mu*(y_z+w_y)

! viscous fluxes in x-direction
f(DENS) = 0.
f(MOM1) = -tau_xx                                       ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(MOM2) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
f(MOM3) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
f(ENER) = -tau_xx*UPrim(VEL1)-tau_xy*UPrim(VEL2)-tau_xz*UPrim(VEL3) &
          -lambda*gradUx(LIFT_TEMP)                     ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) with q_x=-lambda*T_x
! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -tau_xy                                       ! F_euler-mu*(u_y+v_x)
g(MOM2) = -tau_yy                                       ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(MOM3) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
g(ENER) = -tau_xy*UPrim(VEL1)-tau_yy*UPrim(VEL2)-tau_yz*UPrim(VEL3) &
          -lambda*gradUy(LIFT_TEMP)                     ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) with q_y=-lambda*T_y
! viscous fluxes in z-direction
h(DENS) = 0.
h(MOM1) = -tau_xz                                       ! F_euler-mu*(u_z+w_x)
h(MOM2) = -tau_yz                                       ! F_euler-mu*(y_z+w_y)
h(MOM3) = -tau_zz                                       ! F_euler-4/3*mu*w_z+2/3*mu*(u_x+v_y)
h(ENER) = -tau_xz*UPrim(VEL1)-tau_yz*UPrim(VEL2)-tau_zz*UPrim(VEL3) &
          -lambda*gradUz(LIFT_TEMP)                     ! F_euler-(tau_zx*u+tau_zy*v+tau_zz*w-q_z) with q_z=-lambda*T_z
#else
! Precompute entries of shear-stress tensor
tau_xx = mu * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2))  ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = mu * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2))  !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_xy = mu * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))               ! mu*(u_y+v_x)

! viscous fluxes in x-direction
f(DENS) = 0.
f(MOM1) = -tau_xx                                    ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(MOM2) = -tau_xy                                    ! F_euler-mu*(u_y+v_x)
f(MOM3) = 0.
f(ENER) = -tau_xx*UPrim(VEL1)-tau_xy*UPrim(VEL2) &   ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) with q_x=-lambda*T_x
          -lambda*gradUx(LIFT_TEMP)
! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -tau_xy                                    ! F_euler-mu*(u_y+v_x)
g(MOM2) = -tau_yy                                    ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(MOM3) = 0.
g(ENER) = -tau_xy*UPrim(VEL1)-tau_yy*UPrim(VEL2) &   ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) with q_y=-lambda*T_y
          -lambda*gradUy(LIFT_TEMP)
! viscous fluxes in z-direction
h    = 0.
#endif
END SUBROUTINE EvalDiffFlux3D

!==================================================================================================================================
!> Compute Navier-Stokes diffusive flux using the primitive variables and derivatives.
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE EvalNormalDiffFlux3D(UPrim,gradUx,gradUy,gradUz,f,nv,mu,lambda)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim                 !< Solution vector
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx,gradUy,gradUz  !> Gradients in x,y,z directions
REAL,DIMENSION(3             ),INTENT(IN)  :: nv                    !> Normal vector
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: f                     !> Physical fluxes in normal directions
REAL,INTENT(IN)  :: mu      !< viscosity of fluid
REAL,INTENT(IN)  :: lambda  !< thermal conductivity of fluid
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: tau_xx,tau_yy,tau_xy
#if PP_dim==3
REAL                :: tau_zz,tau_xz,tau_yz
#endif
REAL,PARAMETER      :: s23=2./3.
REAL,PARAMETER      :: s43=4./3.
!==================================================================================================================================
#if PP_dim==3
! Precompute entries of shear-stress tensor
tau_xx = mu * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = mu * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2) - s23 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_zz = mu * (-s23 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2) + s43 * gradUz(LIFT_VEL3)) !-2/3*mu*u_x-2/3*mu*v_y +4/3*mu*w*z
tau_xy = mu * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))  ! mu*(u_y+v_x)
tau_xz = mu * (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))  ! mu*(u_z+w_x)
tau_yz = mu * (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))  ! mu*(y_z+w_y)

! viscous fluxes in x-direction
f(DENS) = 0.

f(MOM1) = -nv(1)*tau_xx &                                     ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
          -nv(2)*tau_xy &                                     ! F_euler-mu*(u_y+v_x)
          -nv(3)*tau_xz                                       ! F_euler-mu*(u_z+w_x)
f(MOM2) = -nv(1)*tau_xy &                                     ! F_euler-mu*(u_y+v_x)
          -nv(2)*tau_yy &                                     ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
          -nv(3)*tau_yz                                       ! F_euler-mu*(y_z+w_y)
f(MOM3) = -nv(1)*tau_xz &                                     ! F_euler-mu*(u_z+w_x)
          -nv(2)*tau_yz &                                     ! F_euler-mu*(y_z+w_y)
          -nv(3)*tau_zz                                       ! F_euler-4/3*mu*w_z+2/3*mu*(u_x+v_y)
f(ENER) = -nv(1)*(tau_xx*UPrim(VEL1)+tau_xy*UPrim(VEL2)+tau_xz*UPrim(VEL3)+lambda*gradUx(LIFT_TEMP)) & ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) with q_x=-lambda*T_x
          -nv(2)*(tau_xy*UPrim(VEL1)+tau_yy*UPrim(VEL2)+tau_yz*UPrim(VEL3)+lambda*gradUy(LIFT_TEMP)) & ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) with q_y=-lambda*T_y
          -nv(3)*(tau_xz*UPrim(VEL1)+tau_yz*UPrim(VEL2)+tau_zz*UPrim(VEL3)+lambda*gradUz(LIFT_TEMP))   ! F_euler-(tau_zx*u+tau_zy*v+tau_zz*w-q_z) with q_z=-lambda*T_z
#else
! Precompute entries of shear-stress tensor
tau_xx = mu * ( s43 * gradUx(LIFT_VEL1) - s23 * gradUy(LIFT_VEL2))  ! 4/3*mu*u_x-2/3*mu*v_y -2/3*mu*w*z
tau_yy = mu * (-s23 * gradUx(LIFT_VEL1) + s43 * gradUy(LIFT_VEL2))  !-2/3*mu*u_x+4/3*mu*v_y -2/3*mu*w*z
tau_xy = mu * (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))               ! mu*(u_y+v_x)

! viscous fluxes in x-direction
f(DENS) = 0.
f(MOM1) = -tau_xx                                    ! F_euler-4/3*mu*u_x+2/3*mu*(v_y+w_z)
f(MOM2) = -tau_xy                                    ! F_euler-mu*(u_y+v_x)
f(MOM3) = 0.
f(ENER) = -tau_xx*UPrim(VEL1)-tau_xy*UPrim(VEL2) &   ! F_euler-(tau_xx*u+tau_xy*v+tau_xz*w-q_x) with q_x=-lambda*T_x
          -lambda*gradUx(LIFT_TEMP)
! viscous fluxes in y-direction
g(DENS) = 0.
g(MOM1) = -tau_xy                                    ! F_euler-mu*(u_y+v_x)
g(MOM2) = -tau_yy                                    ! F_euler-4/3*mu*v_y+2/3*mu*(u_x+w_z)
g(MOM3) = 0.
g(ENER) = -tau_xy*UPrim(VEL1)-tau_yy*UPrim(VEL2) &   ! F_euler-(tau_yx*u+tau_yy*v+tau_yz*w-q_y) with q_y=-lambda*T_y
          -lambda*gradUy(LIFT_TEMP)
! viscous fluxes in z-direction
h    = 0.
STOP 'EvalNormalDiffFlux3D not implemented for 2D yet!' ! TODO: Not implemented for 2D yet!
#endif
END SUBROUTINE EvalNormalDiffFlux3D

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single side
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Point(UPrim,gradUx,gradUy,gradUz,f,g,h)
! MODULES
USE MOD_EOS_Vars,ONLY:mu0,cp,Pr
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: f,g,h                !> Physical fluxes in x,y,z directions
#if EDDYVISCOSITY
REAL,DIMENSION(1             ),INTENT(IN)  :: muSGS                !< SGS viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j
REAL                :: mu,lambda
!==================================================================================================================================
mu     = VISCOSITY_PRIM(UPrim)
lambda = THERMAL_CONDUCTIVITY_H(mu)
CALL EvalDiffFlux3D(UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
END SUBROUTINE EvalDiffFlux3D_Point

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single side
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Surface(Nloc,UPrim,gradUx,gradUy,gradUz,f,g,h)
! MODULES
USE MOD_EOS_Vars,ONLY:mu0,cp,Pr
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
REAL                :: mu,lambda
!==================================================================================================================================
DO j=0,ZDIM(Nloc); DO i=0,Nloc
  mu     = VISCOSITY_PRIM(UPrim(:,i,j))
  lambda = THERMAL_CONDUCTIVITY_H(mu)
  CALL EvalDiffFlux3D(UPrim(:,i,j),gradUx(:,i,j),gradUy(:,i,j),gradUz(:,i,j), &
                                        f(:,i,j),     g(:,i,j),     h(:,i,j), &
                                        mu, lambda)
END DO; END DO ! i,j
END SUBROUTINE EvalDiffFlux3D_Surface

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Elem(UPrim,gradUx,gradUy,gradUz,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars,ONLY:mu0,cp,Pr
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim   ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)   :: UPrim                !< Solution vector
REAL,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)   :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar       ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT)  :: f,g,h                !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k
REAL                :: mu,lambda
!==================================================================================================================================
DO k=0,PP_NZ;  DO j=0,PP_N; DO i=0,PP_N
  mu     = VISCOSITY_PRIM(UPrim(:,i,j,k))
  lambda = THERMAL_CONDUCTIVITY_H(mu)
  CALL EvalDiffFlux3D(UPrim(:,i,j,k),gradUx(:,i,j,k),gradUy(:,i,j,k),gradUz(:,i,j,k), &
                                          f(:,i,j,k),     g(:,i,j,k),     h(:,i,j,k), &
                                          mu, lambda)
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalDiffFlux3D_Elem

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single side
!==================================================================================================================================
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE EvalDiffFlux3D_CUDA_Kernel(nDOF,UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)   :: nDOF    !< Polynomial degree of input solution
REAL,VALUE,INTENT(IN)      :: mu      !< viscosity of fluid
REAL,VALUE,INTENT(IN)      :: lambda  !< thermal conductivity of fluid
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,nDOF),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DEVICE,DIMENSION(PP_nVarLifting,nDOF),INTENT(IN)  :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DEVICE,DIMENSION(PP_nVar       ,nDOF),INTENT(OUT) :: f,g,h                !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
!==================================================================================================================================
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF(i.LE.nDOF) CALL EvalDiffFlux3D(UPrim(:,i),gradUx(:,i),gradUy(:,i),gradUz(:,i), &
                                                  f(:,i),     g(:,i),     h(:,i), &
                                                  mu, lambda)
END SUBROUTINE EvalDiffFlux3D_CUDA_Kernel

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single side
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Surface_CUDA(Nloc,UPrim,gradUx,gradUy,gradUz,f,g,h)
! MODULES
USE MOD_EOS_Vars,ONLY:mu0,cp,Pr
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: Nloc  !< Polynomial degree of input solution
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DEVICE,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DEVICE,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: f,g,h                !> Physical fluxes in x,y,z directions
#if EDDYVISCOSITY
REAL,DIMENSION(1             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: muSGS                !< SGS viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: mu,lambda
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
!==================================================================================================================================
!TODO: Implement cleanly!
mu     = VISCOSITY_PRIM(UPrim(:,i,j))
lambda = THERMAL_CONDUCTIVITY_H(mu)

nDOF=(Nloc+1)*(ZDIM(Nloc)+1)
CALL EvalDiffFlux3D_CUDA_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
END SUBROUTINE EvalDiffFlux3D_Surface_CUDA

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single side
!==================================================================================================================================
PPURE SUBROUTINE EvalDiffFlux3D_Sides_CUDA(Nloc,nSides,UPrim,gradUx,gradUy,gradUz,f,g,h)
! MODULES
USE MOD_EOS_Vars,ONLY:mu0,cp,Pr
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: Nloc  !< Polynomial degree of input solution
INTEGER,INTENT(IN)   :: nSides  !< Polynomial degree of input solution
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: UPrim                !< Solution vector
REAL,DEVICE,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DEVICE,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(OUT) :: f,g,h                !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: mu,lambda
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
!==================================================================================================================================
!TODO: Implement cleanly!
mu     = VISCOSITY_PRIM(UPrim(:,i,j,1))
lambda = THERMAL_CONDUCTIVITY_H(mu)

nDOF=(Nloc+1)*(ZDIM(Nloc)+1)*nSides
CALL EvalDiffFlux3D_CUDA_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
END SUBROUTINE EvalDiffFlux3D_Sides_CUDA

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a single volume cell
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Elem_CUDA(UPrim,gradUx,gradUy,gradUz,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars,ONLY:mu0,cp,Pr
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)   :: UPrim                !< Solution vector
REAL,DEVICE,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN)   :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DEVICE,DIMENSION(PP_nVar       ,0:PP_N,0:PP_N,0:PP_NZ),INTENT(OUT)  :: f,g,h                !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: mu,lambda
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
!==================================================================================================================================
!TODO: Implement cleanly!
mu     = VISCOSITY_PRIM(UPrim(:,i,j))
lambda = THERMAL_CONDUCTIVITY_H(mu)

nDOF=(PP_N+1)*(PP_N+1)*(ZDIM(PP_N)+1)
CALL EvalDiffFlux3D_CUDA_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
END SUBROUTINE EvalDiffFlux3D_Elem_CUDA

!==================================================================================================================================
!> Wrapper routine to compute the diffusive part of the Navier-Stokes fluxes for a number of elements
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D_Elems_CUDA(nElems,UPrim,gradUx,gradUy,gradUz,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars,ONLY:mu0,cp,Pr
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: muSGS
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)   :: nElems
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,0:PP_N,0:PP_N,0:PP_NZ,nElems),INTENT(IN)   :: UPrim                !< Solution vector
REAL,DEVICE,DIMENSION(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,nElems),INTENT(IN)   :: gradUx,gradUy,gradUz !> Gradients in x,y,z directions
REAL,DEVICE,DIMENSION(PP_nVar       ,0:PP_N,0:PP_N,0:PP_NZ,nElems),INTENT(OUT)  :: f,g,h                !> Physical fluxes in x,y,z directions
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: mu,lambda
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
!==================================================================================================================================
!TODO: Implement cleanly!
mu     = VISCOSITY_PRIM(UPrim(:,i,j))
lambda = THERMAL_CONDUCTIVITY_H(mu)

nDOF=(PP_N+1)*(PP_N+1)*(ZDIM(PP_N)+1)*nElems
CALL EvalDiffFlux3D_CUDA_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,UPrim,gradUx,gradUy,gradUz,f,g,h,mu,lambda)
END SUBROUTINE EvalDiffFlux3D_Elems_CUDA
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
