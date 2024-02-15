!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
#if FV_ENABLED
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> This module contains only the volume operator of the FV sub-cells method.
!==================================================================================================================================
MODULE MOD_FV_VolInt
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE FV_VolInt
  MODULE PROCEDURE FV_VolInt
  MODULE PROCEDURE FV_VolInt_Conv
END INTERFACE

PUBLIC::FV_VolInt
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Volume operator of the FV sub-cells method.
!> The following steps are performed direction-by-direction (XI/ETA/ZETA) for every inner slice in the respective direction:
!> - reconstruct solution at the sub-cell interfaces
!> - evaluate Riemann solver at the slices
!> - apply fluxes to the left and right sub-cell of the slice
!> are evaluated in the volume integral of the lifting procedure.
!==================================================================================================================================
SUBROUTINE FV_VolInt(UPrim,Ut)
! MODULES
USE MOD_PreProc                         ! all PP_*** variables
USE MOD_Flux         ,ONLY: EvalFlux3D  ! 3D fluxes
USE MOD_FV_Vars
#if VOLINT_VISC
USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
#endif
USE MOD_Riemann      ,ONLY: Riemann
USE MOD_Mesh_Vars    ,ONLY: nElems
USE MOD_EOS          ,ONLY: PrimToCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< solution vector of primitive variables
REAL,INTENT(INOUT) :: Ut(   PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< time derivative of conservative solution vector
                                                                         !< for FV elements
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                       :: i,j,k,p,q,iElem
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_NZ)    :: UPrim_L,UPrim_R
REAL,DIMENSION(PP_nVar    ,0:PP_N,0:PP_NZ)    :: UCons_L,UCons_R
#if VOLINT_VISC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ,0:PP_N) :: diffFlux_x,diffFlux_y
#if PP_dim == 3
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ,0:PP_N) :: diffFlux_z
#endif /*PP_dim == 3*/
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ)        :: Fvisc_FV
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: f,g,h      !< viscous volume fluxes at GP
#endif
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_NZ)        :: F_FV
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: Ut_FV
!==================================================================================================================================

! This routine works as follows:
! The tensor product stucture is used to evaluate the fluxes first for all interfaces/slices in xi-direction, then in eta- and
! at last in zeta-direction.
!   0. The viscous fluxes in all sub-cells are calculated.
! For each direction the following steps are performed.
!   1. copy viscous flux
!   2. reconstruct primitive solution at the cell interfaces
!   3. calculate conservative solution at the interfaces (from the reconstructed primitive solution)
!   4. evaluate Riemann solver to get the advective flux
!   5. calculate viscous flux in normal direction at the interface (mean value of the adjacent viscous fluxes)
!   6. add flux to the left and right sub-cell of the interface

DO iElem=1,nElems
#if FV_ENABLED == 1
  IF (FV_Elems(iElem) .EQ. 0) CYCLE ! DG Element
#elif FV_ENABLED == 2
  IF (FV_alpha(iElem) .LT. EPSILON(0.)) CYCLE
#endif

#if VOLINT_VISC
  ! 0. Eval viscous flux for all sub-cells
  CALL EvalDiffFlux3D(UPrim(:,:,:,:,iElem), gradUx(:,:,:,:,iElem), gradUy(:,:,:,:,iElem), gradUz(:,:,:,:,iElem),f,g,h,iElem)
#endif

  ! === Xi-Direction ============
#if VOLINT_VISC
  ! 1. copy viscous fluxes to temporary array (for performance)
  DO i=0,PP_N
    diffFlux_x(:,:,:,i) = f(:,i,:,:)
    diffFlux_y(:,:,:,i) = g(:,i,:,:)
#if PP_dim == 3
    diffFlux_z(:,:,:,i) = h(:,i,:,:)
#endif
  END DO ! i=0,PP_N
#endif
  ! This is the only nullification we have to do, the rest will be overwritten accordingly. Nullify here to catch only the FV
  ! elements. Might be smarter to do this in a single, long operation?
  Ut_FV(:,0,:,:) = 0.
  ! iterate over all inner slices in xi-direction
  DO i=1,PP_N
    DO q=0,PP_NZ; DO p=0,PP_N
      ! 2. reconstruct solution at left and right side of the interface/slice
#if FV_RECONSTRUCT
      UPrim_L(:,p,q) = UPrim(:,i-1,p,q,iElem) + gradUxi(:,p,q,i-1,iElem) * FV_dx_XI_R(p,q,i-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,i  ,p,q,iElem) - gradUxi(:,p,q,i  ,iElem) * FV_dx_XI_L(p,q,i  ,iElem)
#else
      UPrim_L(:,p,q) = UPrim(:,i-1,p,q,iElem)
      UPrim_R(:,p,q) = UPrim(:,i  ,p,q,iElem)
#endif
    END DO; END DO ! p,q=0,PP_N
    ! 3. convert primitve solution to conservative
    CALL PrimToCons(PP_N,UPrim_L, UCons_L)
    CALL PrimToCons(PP_N,UPrim_R, UCons_R)

    ! 4. calculate advective part of the flux
    CALL Riemann(PP_N,F_FV,UCons_L,UCons_R,UPrim_L,UPrim_R,          &
        FV_NormVecXi (:,:,:,i,iElem), FV_TangVec1Xi(:,:,:,i,iElem), FV_TangVec2Xi(:,:,:,i,iElem))

#if VOLINT_VISC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      Fvisc_FV(:,p,q)=0.5*(FV_NormVecXi(1,p,q,i,iElem)*(diffFlux_x(:,p,q,i-1)+diffFlux_x(:,p,q,i)) &
                          +FV_NormVecXi(2,p,q,i,iElem)*(diffFlux_y(:,p,q,i-1)+diffFlux_y(:,p,q,i)) &
#if PP_dim == 3
                          +FV_NormVecXi(3,p,q,i,iElem)*(diffFlux_z(:,p,q,i-1)+diffFlux_z(:,p,q,i)) &
#endif
                          )
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*VOLINT_VISC*/

    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO k=0,PP_NZ; DO j=0,PP_N
      Ut_FV(:,i-1,j,k) = Ut_FV(:,i-1,j,k) + F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem) * FV_w_inv(i-1)
      ! During our first sweep, the DOF here has never been touched and can thus be overwritten
      Ut_FV(:,i  ,j,k) =               -1.* F_FV(:,j,k) * FV_SurfElemXi_sw(j,k,i,iElem) * FV_w_inv(i)
    END DO; END DO
  END DO ! i

  ! === Eta-Direction ===========
#if VOLINT_VISC
  ! 1. copy fluxes to temporary array (for performance)
  DO j=0,PP_N
    diffFlux_x(:,:,:,j) = f(:,:,j,:)
    diffFlux_y(:,:,:,j) = g(:,:,j,:)
#if PP_dim == 3
    diffFlux_z(:,:,:,j) = h(:,:,j,:)
#endif
  END DO ! j=0,PP_N
#endif

  ! iterate over all inner slices in eta-direction
  DO j=1,PP_N
    DO q=0,PP_NZ; DO p=0,PP_N
      ! 2. reconstruct solution at left and right side of the interface/slice
#if FV_RECONSTRUCT
      UPrim_L(:,p,q) = UPrim(:,p,j-1,q,iElem) + gradUeta(:,p,q,j-1,iElem) * FV_dx_ETA_R(p,q,j-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,j  ,q,iElem) - gradUeta(:,p,q,j  ,iElem) * FV_dx_ETA_L(p,q,j  ,iElem)
#else
      UPrim_L(:,p,q) = UPrim(:,p,j-1,q,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,j  ,q,iElem)
#endif
    END DO; END DO ! p,q=0,PP_N
    ! 3. convert primitve solution to conservative
    CALL PrimToCons(PP_N,UPrim_L, UCons_L)
    CALL PrimToCons(PP_N,UPrim_R, UCons_R)

    ! 4. calculate advective part of the flux
    CALL Riemann(PP_N,F_FV,UCons_L,UCons_R,UPrim_L,UPrim_R,          &
        FV_NormVecEta (:,:,:,j,iElem), FV_TangVec1Eta(:,:,:,j,iElem), FV_TangVec2Eta(:,:,:,j,iElem))
#if VOLINT_VISC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_NZ; DO p=0,PP_N
      Fvisc_FV(:,p,q)=0.5*(FV_NormVecEta(1,p,q,j,iElem)*(diffFlux_x(:,p,q,j-1)+diffFlux_x(:,p,q,j)) &
                          +FV_NormVecEta(2,p,q,j,iElem)*(diffFlux_y(:,p,q,j-1)+diffFlux_y(:,p,q,j)) &
#if PP_dim == 3
                          +FV_NormVecEta(3,p,q,j,iElem)*(diffFlux_z(:,p,q,j-1)+diffFlux_z(:,p,q,j)) &
#endif
                          )
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*VOLINT_VISC*/

    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO k=0,PP_NZ; DO i=0,PP_N
      Ut_FV(:,i,j-1,k) = Ut_FV(:,i,j-1,k) + F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem) * FV_w_inv(j-1)
      Ut_FV(:,i,j  ,k) = Ut_FV(:,i,j  ,k) - F_FV(:,i,k) * FV_SurfElemEta_sw(i,k,j,iElem) * FV_w_inv(j)
    END DO; END DO
  END DO ! j

#if PP_dim == 3
  ! === Zeta-Direction ============
  ! 1. no copy of viscous fluxes to diffFlux_x/y/z required, since f,g,h already have the correct memory layout

  ! iterate over all inner slices in zeta-direction
  DO k=1,PP_N
    DO q=0,PP_N; DO p=0,PP_N
      ! 2. reconstruct solution at left and right side of the interface/slice
#if FV_RECONSTRUCT
      UPrim_L(:,p,q) = UPrim(:,p,q,k-1,iElem) + gradUzeta(:,p,q,k-1,iElem) * FV_dx_ZETA_R(p,q,k-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,q,k  ,iElem) - gradUzeta(:,p,q,k  ,iElem) * FV_dx_ZETA_L(p,q,k  ,iElem)
#else
      UPrim_L(:,p,q) = UPrim(:,p,q,k-1,iElem)
      UPrim_R(:,p,q) = UPrim(:,p,q,k  ,iElem)
#endif
    END DO; END DO ! p,q=0,PP_N
    ! 3. convert primitve solution to conservative
    CALL PrimToCons(PP_N,UPrim_L, UCons_L)
    CALL PrimToCons(PP_N,UPrim_R, UCons_R)

    ! 4. calculate advective part of the flux
    CALL Riemann(PP_N,F_FV,UCons_L,UCons_R,UPrim_L,UPrim_R,          &
        FV_NormVecZeta (:,:,:,k,iElem), FV_TangVec1Zeta(:,:,:,k,iElem), FV_TangVec2Zeta(:,:,:,k,iElem))
#if VOLINT_VISC
    ! 5. compute viscous flux in normal direction of the interface
    DO q=0,PP_N; DO p=0,PP_N
      Fvisc_FV(:,p,q)=0.5*(FV_NormVecZeta(1,p,q,k,iElem)*(f(:,p,q,k-1)+f(:,p,q,k)) &
                          +FV_NormVecZeta(2,p,q,k,iElem)*(g(:,p,q,k-1)+g(:,p,q,k)) &
                          +FV_NormVecZeta(3,p,q,k,iElem)*(h(:,p,q,k-1)+h(:,p,q,k)))
    END DO; END DO
    F_FV = F_FV + Fvisc_FV
#endif /*VOLINT_VISC*/

    ! 6. apply flux to the sub-cells at the left and right side of the interface/slice
    DO j=0,PP_N; DO i=0,PP_N
      Ut_FV(:,i,j,k-1) = Ut_FV(:,i,j,k-1) + F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem) * FV_w_inv(k-1)
      Ut_FV(:,i,j,k  ) = Ut_FV(:,i,j,k  ) - F_FV(:,i,j) * FV_SurfElemZeta_sw(i,j,k,iElem) * FV_w_inv(k)
    END DO; END DO
  END DO ! k
#endif /* PP_dim == 3 */

#if FV_ENABLED == 2
  ! Blend the solutions together
  Ut(:,:,:,:,iElem) = (1 - FV_alpha(iElem)) * Ut(:,:,:,:,iElem) + FV_alpha(iElem)*Ut_FV
#else
  Ut(:,:,:,:,iElem) = Ut_FV
#endif /*FV_BLENDING*/

END DO ! iElem
END SUBROUTINE FV_VolInt

!==================================================================================================================================
!> Volume operator of the FV sub-cells method.
!> The following steps are performed direction-by-direction (XI/ETA/ZETA) for every inner slice in the respective direction:
!> - reconstruct solution at the sub-cell interfaces
!> - evaluate Riemann solver at the slices
!> - apply fluxes to the left and right sub-cell of the slice
!> are evaluated in the volume integral of the lifting procedure.
!==================================================================================================================================
SUBROUTINE FV_VolInt_Conv(d_U,d_UPrim,d_Ut,streamID)
! MODULES
USE CUDAFOR
USE MOD_GPU
USE MOD_PreProc                         ! all PP_*** variables
USE MOD_FV_Vars
USE MOD_DG_Vars      ,ONLY: nDOFElem
USE MOD_Mesh_Vars    ,ONLY: nElems
USE MOD_EOS_Vars     ,ONLY: d_EOS_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DEVICE,INTENT(IN)    :: d_U(    PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< solution vector of primitive variables
REAL,DEVICE,INTENT(IN)    :: d_UPrim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< solution vector of primitive variables
REAL,DEVICE,INTENT(INOUT) :: d_Ut(   PP_nVar    ,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)  !< time derivative of conservative solution vector
INTEGER(KIND=CUDA_STREAM_KIND),OPTIONAL,INTENT(IN) :: streamID
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
INTEGER(KIND=CUDA_STREAM_KIND) :: mystream
!==================================================================================================================================
mystream=DefaultStream
IF (PRESENT(streamID)) mystream=streamID

nDOF=nDOFElem*nElems
CALL FV_VolInt_Conv_GPU<<<nDOF/nThreads+1,nThreads,0,mystream>>>(nElems,PP_N,d_U,d_UPrim,d_Ut &
                                                     ,d_FV_NormVecXi  ,d_FV_TangVec1Xi  ,d_FV_TangVec2Xi   &
                                                     ,d_FV_NormVecEta ,d_FV_TangVec1Eta ,d_FV_TangVec2Eta  &
                                                     ,d_FV_NormVecZeta,d_FV_TangVec1Zeta,d_FV_TangVec2Zeta &
                                                     ,d_FV_SurfElemXi_sw   &
                                                     ,d_FV_SurfElemEta_sw  &
                                                     ,d_FV_SurfElemZeta_sw &
                                                     ,d_FV_w_inv &
                                                     ,d_EOS_Vars &
#if FV_ENABLED==2
                                                     ,d_FV_alpha &
#endif
                                                     )

END SUBROUTINE FV_VolInt_Conv

!==================================================================================================================================
!> Volume operator of the FV sub-cells method.
!> The following steps are performed direction-by-direction (XI/ETA/ZETA) for every inner slice in the respective direction:
!> - reconstruct solution at the sub-cell interfaces
!> - evaluate Riemann solver at the slices
!> - apply fluxes to the left and right sub-cell of the slice
!> are evaluated in the volume integral of the lifting procedure.
!==================================================================================================================================
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE FV_VolInt_Conv_GPU(nElems,Nloc,U,UPrim,Ut &
                                                      ,FV_NormVecXi  ,FV_TangVec1Xi  ,FV_TangVec2Xi   &
                                                      ,FV_NormVecEta ,FV_TangVec1Eta ,FV_TangVec2Eta  &
                                                      ,FV_NormVecZeta,FV_TangVec1Zeta,FV_TangVec2Zeta &
                                                      ,FV_SurfElemXi_sw   &
                                                      ,FV_SurfElemEta_sw  &
                                                      ,FV_SurfElemZeta_sw &
                                                      ,FV_w_inv &
                                                      ,EOS_Vars &
#if FV_ENABLED==2
                                                      ,FV_alpha &
#endif
                                                      )
! MODULES
USE MOD_Riemann      ,ONLY: Riemann
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)  :: nElems
INTEGER,VALUE,INTENT(IN)  :: Nloc
REAL,DEVICE,INTENT(IN)    :: UPrim(PP_nVarPrim ,0:Nloc,0:Nloc ,0:ZDIM(Nloc),1:nElems)  !< solution vector of primitive variables
REAL,DEVICE,INTENT(IN)    :: U(    PP_nVar     ,0:Nloc,0:Nloc ,0:ZDIM(Nloc),1:nElems)  !< solution vector of conservative variables
REAL,DEVICE,INTENT(INOUT) :: Ut(   PP_nVar     ,0:Nloc,0:Nloc ,0:ZDIM(Nloc),1:nElems)  !< time derivative of conservative solution vector
REAL,DEVICE,INTENT(IN)    :: FV_NormVecXi(    3,0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_TangVec1Xi(   3,0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_TangVec2Xi(   3,0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_SurfElemXi_sw(  0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_NormVecEta(   3,0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_TangVec1Eta(  3,0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_TangVec2Eta(  3,0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_SurfElemEta_sw( 0:Nloc,0:ZDIM(Nloc),1:Nloc ,1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_NormVecZeta(  3,0:Nloc,0:Nloc ,ZDIM(1:Nloc),1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_TangVec1Zeta( 3,0:Nloc,0:Nloc ,ZDIM(1:Nloc),1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_TangVec2Zeta( 3,0:Nloc,0:Nloc ,ZDIM(1:Nloc),1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_SurfElemZeta_sw(0:Nloc,0:Nloc ,ZDIM(1:Nloc),1:nElems)
REAL,DEVICE,INTENT(IN)    :: FV_w_inv(0:Nloc)
REAL,DEVICE,INTENT(IN)    :: EOS_Vars(PP_nVarEOS)
#if FV_ENABLED==2
REAL,DEVICE,INTENT(IN)    :: FV_alpha(nElems)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k,iElem,threadID,rest
REAL,DIMENSION(PP_nVar) :: Ut_FV_tmp
REAL,DIMENSION(PP_nVar) :: F_FV
!==================================================================================================================================
! Get thread indices
threadID = (blockidx%x-1) * blockdim%x + threadidx%x
! Get iElem of current thread
iElem   =        (threadID-1)/(Nloc+1)**3+1 ! Elems are 1-indexed
rest    = threadID-(iElem-1)*(Nloc+1)**3
! Get ijk indices of current thread
k       = (rest-1)/(Nloc+1)**2
rest    =  rest- k*(Nloc+1)**2
j       = (rest-1)/(Nloc+1)!**1
rest    =  rest- j*(Nloc+1)!**1
i       = (rest-1)!/(Nloc+1)**0

IF (iElem.LE.nElems) THEN
  ! We have point i,j,k!
  Ut_FV_tmp(:) = 0.

  ! === Xi-Direction ===========
  ! Left
  IF (i.GT.0) THEN ! Flux across surface of DG element is performed in surface integral
    CALL Riemann(F_FV &
                ,U(    :,i-1,j,k,iElem) & ! Left
                ,U(    :,i  ,j,k,iElem) & ! Right
                ,UPrim(:,i-1,j,k,iElem) & ! Left
                ,UPrim(:,i  ,j,k,iElem) & ! Right
                ,FV_NormVecXi( :,j,k,i,iElem) &
                ,FV_TangVec1Xi(:,j,k,i,iElem) &
                ,FV_TangVec2Xi(:,j,k,i,iElem) &
                ,EOS_Vars)
    ! sign contains normal vector in negative xi-direction
    Ut_FV_tmp(:) = Ut_FV_tmp(:) - F_FV(:) * FV_SurfElemXi_sw(j,k,i,iElem) * FV_w_inv(i)
  END IF

  ! Right
  IF (i.LT.Nloc) THEN ! Flux across surface of DG element is performed in surface integral
    CALL Riemann(F_FV &
                ,U(    :,i  ,j,k,iElem) & ! Left
                ,U(    :,i+1,j,k,iElem) & ! Right
                ,UPrim(:,i  ,j,k,iElem) & ! Left
                ,UPrim(:,i+1,j,k,iElem) & ! Right
                ,FV_NormVecXi( :,j,k,i+1,iElem) &
                ,FV_TangVec1Xi(:,j,k,i+1,iElem) &
                ,FV_TangVec2Xi(:,j,k,i+1,iElem) &
                ,EOS_Vars)
    ! sign contains normal vector in positive xi-direction
    Ut_FV_tmp(:) = Ut_FV_tmp(:) + F_FV(:) * FV_SurfElemXi_sw(j,k,i+1,iElem) * FV_w_inv(i)
  END IF

  ! === Eta-Direction ===========
  ! Left
  IF (j.GT.0) THEN ! Flux across surface of DG element is performed in surface integral
    CALL Riemann(F_FV &
                ,U(    :,i,j-1,k,iElem) & ! Left
                ,U(    :,i,j  ,k,iElem) & ! Right
                ,UPrim(:,i,j-1,k,iElem) & ! Left
                ,UPrim(:,i,j  ,k,iElem) & ! Right
                ,FV_NormVecEta( :,i,k,j,iElem) &
                ,FV_TangVec1Eta(:,i,k,j,iElem) &
                ,FV_TangVec2Eta(:,i,k,j,iElem) &
                ,EOS_Vars)
    ! sign contains normal vector in negative eta-direction
    Ut_FV_tmp(:) = Ut_FV_tmp(:) - F_FV(:) * FV_SurfElemEta_sw(i,k,j,iElem) * FV_w_inv(j)
  END IF

  ! Right
  IF (j.LT.Nloc) THEN ! Flux across surface of DG element is performed in surface integral
    CALL Riemann(F_FV &
                ,U(    :,i,j  ,k,iElem) & ! Left
                ,U(    :,i,j+1,k,iElem) & ! Right
                ,UPrim(:,i,j  ,k,iElem) & ! Left
                ,UPrim(:,i,j+1,k,iElem) & ! Right
                ,FV_NormVecEta( :,i,k,j+1,iElem) &
                ,FV_TangVec1Eta(:,i,k,j+1,iElem) &
                ,FV_TangVec2Eta(:,i,k,j+1,iElem) &
                ,EOS_Vars)
    ! sign contains normal vector in positive eta-direction
    Ut_FV_tmp(:) = Ut_FV_tmp(:) + F_FV(:) * FV_SurfElemEta_sw(i,k,j+1,iElem) * FV_w_inv(j)
  END IF

#if PP_dim == 3
  ! === Zeta-Direction ===========
  ! Left
  IF (k.GT.0) THEN ! Flux across surface of DG element is performed in surface integral
    CALL Riemann(F_FV &
                ,U(    :,i,j,k-1,iElem) & ! Left
                ,U(    :,i,j,k  ,iElem) & ! Right
                ,UPrim(:,i,j,k-1,iElem) & ! Left
                ,UPrim(:,i,j,k  ,iElem) & ! Right
                ,FV_NormVecZeta( :,i,j,k,iElem) &
                ,FV_TangVec1Zeta(:,i,j,k,iElem) &
                ,FV_TangVec2Zeta(:,i,j,k,iElem) &
                ,EOS_Vars)
    ! sign contains normal vector in negative xi-direction
    Ut_FV_tmp(:) = Ut_FV_tmp(:) - F_FV(:) * FV_SurfElemZeta_sw(i,j,k,iElem) * FV_w_inv(k)
  END IF

  ! Right
  IF (k.LT.Nloc) THEN ! Flux across surface of DG element is performed in surface integral
    CALL Riemann(F_FV &
                ,U(    :,i,j,k  ,iElem) & ! Left
                ,U(    :,i,j,k+1,iElem) & ! Right
                ,UPrim(:,i,j,k  ,iElem) & ! Left
                ,UPrim(:,i,j,k+1,iElem) & ! Right
                ,FV_NormVecZeta( :,i,j,k,iElem) &
                ,FV_TangVec1Zeta(:,i,j,k,iElem) &
                ,FV_TangVec2Zeta(:,i,j,k,iElem) &
                ,EOS_Vars)
    ! sign contains normal vector in positive xi-direction
    Ut_FV_tmp(:) = Ut_FV_tmp(:) + F_FV(:) * FV_SurfElemZeta_sw(i,j,k+1,iElem) * FV_w_inv(k)
  END IF
#endif

#if FV_ENABLED == 2
  ! Blend the solutions together
  Ut(:,i,j,k,iElem) = (1 - FV_alpha(iElem)) * Ut(:,i,j,k,iElem) + FV_alpha(iElem)*Ut_FV_tmp(:)
#else
  Ut(:,i,j,k,iElem) = Ut_FV_tmp(:)
#endif /*FV_BLENDING*/

END IF
END SUBROUTINE FV_VolInt_Conv_GPU

END MODULE MOD_FV_VolInt
#endif /* FV_ENABLED */
