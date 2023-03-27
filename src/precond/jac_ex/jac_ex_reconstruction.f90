!!=================================================================================================================================
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
!==================================================================================================================================
!> Contains the routines for the derivative of the fv_reconstruction procedure
!==================================================================================================================================
MODULE MOD_Jac_Ex_Reconstruction
#if FV_ENABLED && FV_RECONSTRUCT
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE Fill_ExtendedState
  MODULE PROCEDURE Fill_ExtendedState
END INTERFACE

INTERFACE FV_Reconstruction_Derivative
  MODULE PROCEDURE FV_Reconstruction_Derivative
END INTERFACE

INTERFACE FV_Reconstruction_Derivative_Surf
  MODULE PROCEDURE FV_Reconstruction_Derivative_Surf
END INTERFACE

#if PARABOLIC
INTERFACE JacFVGradients_Vol
  MODULE PROCEDURE JacFVGradients_Vol
END INTERFACE

INTERFACE JacFVGradients_nb
  MODULE PROCEDURE JacFVGradients_nb
END INTERFACE
#endif

PUBLIC::Fill_ExtendedState,FV_Reconstruction_Derivative,FV_Reconstruction_Derivative_Surf
#if PARABOLIC
PUBLIC::JacFVGradients_Vol,JacFVGradients_nb
#endif
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Fills extended arrays containing
!> - volume value per element
!> - values at the nodes of the first inner layer next to the DG element interface 
!> Those arrays are: UPrim_extended (containing primitive solution vector) and
!> FV_sdx_XI/ETA/ZETA_extended containing the inverse of the distance between subcell centers.
!==================================================================================================================================
SUBROUTINE Fill_ExtendedState(t,URec,URec_extended,FV_sdx_XI_extended,FV_sdx_ETA_extended &
#if PP_dim ==3 
                              ,FV_sdx_ZETA_extended &
#endif
                              )
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ProlongToFacePrim  ,ONLY: ProlongToFacePrim
USE MOD_Interpolation_Vars ,ONLY: L_Minus,L_Plus
USE MOD_Mesh_Vars          ,ONLY: nElems,nSides,firstInnerSide,nBCSides,lastMPISide_MINE,firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars          ,ONLY: S2V,SideToElem,NormVec,TangVec1,TangVec2,Face_xGP
USE MOD_GetBoundaryFlux    ,ONLY: GetBoundaryState
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisSurf
#if USE_MPI
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI_Vars           ,ONLY: MPIRequest_Rec_SM,MPIRequest_Rec_MS,nNbProcs,DataSizeSidePrim
#endif
USE MOD_FV_Vars            ,ONLY: FV_Elems_Sum,FV_sdx_Face,FV_Vdm,FV_Elems
USE MOD_FV_Vars            ,ONLY: FV_sdx_XI,FV_sdx_ETA
#if PP_dim == 3
USE MOD_FV_Vars            ,ONLY: FV_sdx_ZETA
#endif
USE MOD_FillMortarPrim     ,ONLY: U_MortarPrim
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)     :: t                                                                 !< physical time
REAL,INTENT(IN)     :: URec(         PP_nVarPrim, 0:PP_N,   0:PP_N,   0:PP_NZ, 1:nElems) !< rec volume solution
REAL,INTENT(OUT)    :: FV_sdx_XI_extended(        0:PP_N,   0:PP_NZ,  0:PP_N+1,1:nElems) !< FV_sdx_XI:storage order   (j,k,i,iElem)
REAL,INTENT(OUT)    :: FV_sdx_ETA_extended(       0:PP_N,   0:PP_NZ,  0:PP_N+1,1:nElems) !< FV_sdx_ETA:storage order  (i,k,j,iElem)
#if PP_dim == 3
REAL,INTENT(OUT)    :: FV_sdx_ZETA_extended(      0:PP_N,   0:PP_N,   0:PP_N+1,1:nElems) !< FV_sdx_ZETA:storage order (i,j,k,iElem)
REAL,INTENT(OUT)    :: URec_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,1:nElems) !< extended rec volume solution
#else
REAL,INTENT(OUT)    :: URec_extended(PP_nVarPrim,-1:PP_N+1,-1:PP_N+1, 0:PP_NZ, 1:nElems) !< extended rec volume solution
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:ZDIM(PP_N),1:nSides)  :: URec_master,URec_slave
INTEGER                                                   :: p,q,l,i,j,k,ijk(3),locSideID,ElemID,SideID,flip,iElem
REAL,DIMENSION(0:PP_N,0:PP_NZ,1:nSides)                   :: FV_sdx_Face_loc
REAL,DIMENSION(0:PP_N,0:PP_NZ)                            :: FV_sdx_Face_tmp
REAL,DIMENSION(1:PP_nVarPrim,0:PP_N,0:PP_NZ)              :: UPrim_tmp
!==================================================================================================================================
! 1. Fill inner data with element volume data
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    URec_extended(:,i,j,k,iElem) = URec(:,i,j,k,iElem)
  END DO; END DO; END DO! i,j,k=0,PP_N
END DO ! iElem

! 2.   Add volume data of first layer of neighbouring cell
! 2.1. Copy first layer to surface and send data from slave to master and vice versa
#if USE_MPI
CALL StartReceiveMPIData(URec_slave ,DataSizeSidePrim,1,nSides,MPIRequest_Rec_SM(:,SEND),SendID=2) ! slave -> master
CALL StartReceiveMPIData(URec_master,DataSizeSidePrim,1,nSides,MPIRequest_Rec_MS(:,SEND),SendID=1) ! master -> slave
CALL ProlongToFacePrim(PP_N,URec,URec_master,URec_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_MortarPrim(URec_master,URec_slave,doMPISides=.TRUE.) ! fill faces at mortar interface with dummy values
CALL StartSendMPIData(   URec_slave ,DataSizeSidePrim,1,nSides,MPIRequest_Rec_SM(:,RECV),SendID=2) ! slave -> master
#endif
CALL ProlongToFacePrim(PP_N,URec,URec_master,URec_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
CALL U_MortarPrim(URec_master,URec_slave,doMPISides=.FALSE.)
#if USE_MPI
! prolongation has to be finished before communication of URec_master as doMPISides=.FALSE. includes MPISidesMINE
CALL StartSendMPIData(   URec_master,DataSizeSidePrim,1,nSides,MPIRequest_Rec_MS(:,RECV),SendID=1) ! master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_SM)
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_MS)
#endif

! 3.   Fill information about distance between solution points at faces (FV_sdx_Face)
! 3.2. Send FV_sdx_Face from master to slave
FV_sdx_Face_loc = 0.
DO SideID=firstInnerSide,lastMPISide_MINE
  IF (FV_Elems_Sum(SideID).GT.0) FV_sdx_Face_loc(:,:,SideID) = FV_sdx_Face(:,:,FV_Elems_Sum(SideID),SideID)
END DO
#if USE_MPI
CALL StartReceiveMPIData(FV_sdx_Face_loc,(PP_N+1)*(PP_NZ+1),1,nSides,MPIRequest_Rec_MS(:,SEND),SendID=1) ! master -> slave
CALL StartSendMPIData(   FV_sdx_Face_loc,(PP_N+1)*(PP_NZ+1),1,nSides,MPIRequest_Rec_MS(:,RECV),SendID=1) ! master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Rec_MS)
#endif
! 3.3. Fill FV_sdx_extended array with inner data
DO iElem=1,nElems
  DO p=0,PP_N; DO q=0,PP_NZ; DO l=1,PP_N
    FV_sdx_XI_extended(  p,q,l,iElem) = FV_sdx_XI(  p,q,l,iElem)
    FV_sdx_ETA_extended( p,q,l,iElem) = FV_sdx_ETA( p,q,l,iElem)
#if PP_dim == 3
    FV_sdx_ZETA_extended(p,q,l,iElem) = FV_sdx_ZETA(p,q,l,iElem)
#endif
  END DO; END DO; END DO
END DO

! 2.2. Switch solution representation of URec_slave/master to FV if corresponding element is a DG element
! 2.3. Add first layer of URec of neighbouring cell to extended array
! 3.4. Add face data to extended FV_sdx array (take care of storage order of FV_sdx!)
DO SideID=1,nSides
  ! Loop over master and slave sides
  DO i = 1, 2
    IF (i.EQ.1) THEN ! master sides
      ElemID = SideToElem(S2E_ELEM_ID,SideID)
      IF(ElemID.EQ.-1) CYCLE ! skip if my proc does not have that element
      IF(FV_Elems(ElemID).EQ.0) CYCLE ! skip for DG elements
      locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
      flip      = 0
      IF (SideID.LE.nBCSides) THEN
        FV_sdx_face_tmp = FV_sdx_face(:,:,3,SideID)
        ! set state with boundary condition
        CALL GetBoundaryState(SideID,t,PP_N,UPrim_tmp,URec_master(:,:,:,SideID),                    &
                              NormVec(:,:,:,1,SideID),TangVec1(:,:,:,1,SideID),TangVec2(:,:,:,1,SideID), &
                              Face_xGP(:,:,:,1,SideID))
      ELSE
        IF ((SideID.LT.firstInnerSide).OR.((SideID.GE.firstMortarMPISide).AND.(SideID.LE.lastMortarMPISide))) THEN ! big mortar side
          FV_sdx_face_tmp = 0. ! derivative of reconstruction of gradients influenced by mortar interface are set to zero.
          UPrim_tmp = URec_master(:,:,:,SideID)
        ELSE
          FV_sdx_face_tmp = FV_sdx_face_loc(:,:,SideID)
          IF(FV_Elems_Sum(SideID).EQ.1) THEN
            CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_Vdm,URec_slave(:,:,:,SideID),UPrim_tmp(:,:,:)) ! switch slave to FV
          ELSE
            UPrim_tmp = URec_slave(:,:,:,SideID)
          END IF
        END IF
      END IF
    ELSE ! slave side
      ElemID = SideToElem(S2E_NB_ELEM_ID,SideID)
      IF(ElemID.EQ.-1) CYCLE ! skip if my proc does not have that element
      IF(FV_Elems(ElemID).EQ.0) CYCLE ! skip for DG elements
      locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
      flip      = SideToElem(S2E_FLIP,SideID)
      FV_sdx_face_tmp = FV_sdx_face_loc(:,:,SideID)
      IF(FV_Elems_Sum(SideID).EQ.2) THEN
        CALL ChangeBasisSurf(PP_nVarPrim,PP_N,PP_N,FV_Vdm,URec_master(:,:,:,SideID),UPrim_tmp(:,:,:)) ! switch master to FV
      ELSE
        UPrim_tmp = URec_master(:,:,:,SideID)
      END IF
    END IF

    DO q=0,PP_NZ; DO p=0,PP_N
      ijk=S2V(:,0,p,q,flip,locSideID)
      SELECT CASE(locSideID)
      CASE(XI_MINUS)
        URec_extended(:,     ijk(1)-1,ijk(2)  ,ijk(3)  ,ElemID) = UPrim_tmp(    :,p,q)
        FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)  ,ElemID) = FV_sdx_face_tmp(p,q)
      CASE(XI_PLUS)
        URec_extended(:,     ijk(1)+1,ijk(2)  ,ijk(3)  ,ElemID) = UPrim_tmp(    :,p,q)
        FV_sdx_XI_extended(  ijk(2)  ,ijk(3)  ,ijk(1)+1,ElemID) = FV_sdx_face_tmp(p,q)
      CASE(ETA_MINUS)
        URec_extended(:,     ijk(1)  ,ijk(2)-1,ijk(3)  ,ElemID) = UPrim_tmp(    :,p,q)
        FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)  ,ElemID) = FV_sdx_face_tmp(p,q)
      CASE(ETA_PLUS)
        URec_extended(:,     ijk(1)  ,ijk(2)+1,ijk(3)  ,ElemID) = UPrim_tmp(    :,p,q)
        FV_sdx_ETA_extended( ijk(1)  ,ijk(3)  ,ijk(2)+1,ElemID) = FV_sdx_face_tmp(p,q)
#if PP_dim == 3
      CASE(ZETA_MINUS)
        URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)-1,ElemID) = UPrim_tmp(    :,p,q)
        FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)  ,ElemID) = FV_sdx_face_tmp(p,q)
      CASE(ZETA_PLUS)
        URec_extended(:,     ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = UPrim_tmp(    :,p,q)
        FV_sdx_ZETA_extended(ijk(1)  ,ijk(2)  ,ijk(3)+1,ElemID) = FV_sdx_face_tmp(p,q)
#endif
      END SELECT
    END DO; END DO
  END DO ! i = 1, 2
END DO
  
END SUBROUTINE Fill_ExtendedState

!===================================================================================================================================
!> Calculates DOF-wise Jacobian of reconstruction procedure
!> dF/dUvol = dF/dU_LR * dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim * dUvol_prim/dUvol
!>               |            |                       |                   |
!>  FD in DGVolIntJac_FV    dCons/dPrim   derivative of reconstruction  dPrim/dCons
!>                          |_____________________________________________________|
!>                                                    |
!>                                              calculated here
!===================================================================================================================================
SUBROUTINE FV_Reconstruction_Derivative(FV_sdx,FV_dx_L,FV_dx_R,UPrim_plus,UPrim_minus, &
                                        URec_extended,dUdUvol_plus,dUdUvol_minus)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_FV_Limiter     ,ONLY: FV_Limiter
USE MOD_FV_Vars        ,ONLY: LimiterType
USE MOD_Jacobian       ,ONLY: dPrimTempdCons,dConsdPrimTemp
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: FV_sdx(0:PP_N+1)                          !< extended inverse of distance between centers of FV subcells
REAL,INTENT(IN)    :: FV_dx_L(0:PP_N)                           !< left distance between face and center of FV subcells
REAL,INTENT(IN)    :: FV_dx_R(0:PP_N)                           !< right distance between face and center of FV subcells
REAL,INTENT(IN)    :: UPrim_plus(PP_nVarPrim,0:PP_N)            !< primitive reconstructed value at plus side
REAL,INTENT(IN)    :: UPrim_minus(PP_nVarPrim,0:PP_N)           !< primitive reconstructed value at minus side
REAL,INTENT(IN)    :: URec_extended(PP_nVarPrim,-1:PP_N+1)      !< extended rec volume solution
REAL,INTENT(OUT)   :: dUdUvol_plus( PP_nVar,PP_nVar,0:PP_N,1:3) !< Jacobian of reconstruction procedure of left  state of interface
REAL,INTENT(OUT)   :: dUdUvol_minus(PP_nVar,PP_nVar,0:PP_N,1:3) !< Jacobian of reconstruction procedure of right state of interface
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                            :: iVar,i,ind
REAL,DIMENSION(PP_nVarPrim,0:PP_N)                 :: s_L,s_R,s_lim
REAL,DIMENSION(PP_nVar,PP_nVarPrim, 0:PP_N  )      :: Jac_ConsPrim_plus,Jac_ConsPrim_minus
REAL,DIMENSION(PP_nVarPrim,PP_nVar,-1:PP_N+1)      :: Jac_PrimCons
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,1:3) :: dUdUvolprim_plus
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,0:PP_N,1:3) :: dUdUvolprim_minus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,1:3) :: matrix_plus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,0:PP_N,1:3) :: matrix_minus
!==================================================================================================================================
! storage order in dUdUvol(1:3):
! dU_LR_i/dU_(i-1), dU_LR_i/dU_i, dU_LR_i/dU_(i+1)
dUdUvol_plus  = 0.
dUdUvol_minus = 0.
dUdUvolprim_plus  = 0.
dUdUvolprim_minus = 0.

! 1. Do derivative of i_th surface data from cons to prim
DO i=0,PP_N
  CALL dConsdPrimTemp(UPrim_plus( :,i),Jac_ConsPrim_plus( :,:,i))
  CALL dConsdPrimTemp(UPrim_minus(:,i),Jac_ConsPrim_minus(:,:,i))
END DO
! 2. Do derivative of i-1:i+1 volume data from prim to cons 
DO i=-1,PP_N+1
  CALL dPrimTempdCons(URec_extended(:,i),Jac_PrimCons(:,:,i))
END DO

! 3. Calculate derivative of reconstruction
DO i=0,PP_N
  ! Calculate left and right unlimited slopes
  s_L(:,i) = (URec_extended(:,i  ) - URec_extended(:,i-1)) * FV_sdx(i  ) 
  s_R(:,i) = (URec_extended(:,i+1) - URec_extended(:,i  )) * FV_sdx(i+1)
  ! Limit slopes
  CALL FV_Limiter(s_L(:,i),s_R(:,i),s_lim(:,i)) ! only null- or minmod-limiter
  SELECT CASE(LimiterType)
  CASE(0) ! NullLimiter
    DO iVar=1,PP_nVarPrim
      dUdUvolprim_minus(iVar,iVar,i,2) = 1.
      dUdUvolprim_plus( iVar,iVar,i,2) = 1.
      ! NullLimiter has no dependencies on neighbouring DOFs!
    END DO !iVar
  CASE(1) ! MinMod
    DO iVar=1,PP_nVarPrim
      IF(s_lim(iVar,i).EQ.0.)THEN ! first order
        dUdUvolprim_plus( iVar,iVar,i,2) = 1.
        dUdUvolprim_minus(iVar,iVar,i,2) = 1.
      ELSEIF(s_lim(iVar,i).EQ.s_L(iVar,i))THEN ! use left slope
        ! only depends on DOF on my left and myself
        dUdUvolprim_plus( iVar,iVar,i,1) = 0. - FV_sdx(i  ) * FV_dx_R(i)
        dUdUvolprim_plus( iVar,iVar,i,2) = 1. + FV_sdx(i  ) * FV_dx_R(i)

        dUdUvolprim_minus(iVar,iVar,i,1) = 0. + FV_sdx(i  ) * FV_dx_L(i)
        dUdUvolprim_minus(iVar,iVar,i,2) = 1. - FV_sdx(i  ) * FV_dx_L(i)
      ELSEIF(s_lim(iVar,i).EQ.s_R(iVar,i))THEN ! use right slope
        ! only depends on DOF on my right and myself
        dUdUvolprim_plus( iVar,iVar,i,2) = 1. - FV_sdx(i+1) * FV_dx_R(i)
        dUdUvolprim_plus( iVar,iVar,i,3) = 0. + FV_sdx(i+1) * FV_dx_R(i)

        dUdUvolprim_minus(iVar,iVar,i,2) = 1. + FV_sdx(i+1) * FV_dx_L(i)
        dUdUvolprim_minus(iVar,iVar,i,3) = 0. - FV_sdx(i+1) * FV_dx_L(i)
      ELSE
        CALL Abort(__STAMP__,'Slopes do not match with minmod in preconditioner!')
      END IF
    END DO !iVar
  CASE(3) ! van Albada
    DO iVar=1,PP_nVarPrim
      IF((s_R(iVar,i)*s_L(iVar,i)).LT.0.)THEN ! first order
        dUdUvolprim_plus( iVar,iVar,i,2) = 1.
        dUdUvolprim_minus(iVar,iVar,i,2) = 1.
      ELSE
        dUdUvolprim_plus( iVar,iVar,i,1) = 0. + FV_dx_R(i) * FV_sdx(i) * (&
                       -(s_R(iVar,i)*(s_L(iVar,i)+s_R(iVar,i))+s_L(iVar,i)*s_R(iVar,i))/MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13) &
                       +(2.*s_L(iVar,i)*(s_L(iVar,i)**2*s_R(iVar,i)+s_R(iVar,i)**2*s_L(iVar,i)))/(MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13)**2))
    
        dUdUvolprim_plus( iVar,iVar,i,2) = 1. + FV_dx_R(i) * (&
                       -(s_L(iVar,i)**2*FV_sdx(i+1)-s_R(iVar,i)**2*FV_sdx(i) & 
                         +2.*s_L(iVar,i)*s_R(iVar,i)*FV_sdx(i+1)-2.*s_L(iVar,i)*s_R(iVar,i)*FV_sdx(i))/MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13) &
                       -(2.*(s_L(iVar,i)*s_R(iVar,i)**2+s_L(iVar,i)**2*s_R(iVar,i)) &
                         *(s_L(iVar,i)*FV_sdx(i)-s_R(iVar,i)*FV_sdx(i+1)))/(MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13)**2))
    
        dUdUvolprim_plus( iVar,iVar,i,3) = 0. + FV_dx_R(i) * FV_sdx(i+1) * (&
                       -(2*s_L(iVar,i)*s_R(iVar,i)**2*(s_L(iVar,i)+s_R(iVar,i)))/(MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13)**2) &
                       +(s_L(iVar,i)*(2.*s_R(iVar,i)+s_L(iVar,i)))/MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13))
    
        dUdUvolprim_minus(iVar,iVar,i,1) = 0. - FV_dx_L(i) * FV_sdx(i) * (&
                       -(s_R(iVar,i)*(s_L(iVar,i)+s_R(iVar,i))+s_L(iVar,i)*s_R(iVar,i))/MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13) &
                       +(2.*s_L(iVar,i)*(s_L(iVar,i)**2*s_R(iVar,i)+s_R(iVar,i)**2*s_L(iVar,i)))/(MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13)**2))
    
        dUdUvolprim_minus(iVar,iVar,i,2) = 1. - FV_dx_L(i) * (&
                       -(s_L(iVar,i)**2*FV_sdx(i+1)-s_R(iVar,i)**2*FV_sdx(i) & 
                         +2.*s_L(iVar,i)*s_R(iVar,i)*FV_sdx(i+1)-2.*s_L(iVar,i)*s_R(iVar,i)*FV_sdx(i))/MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13) &
                       -(2.*(s_L(iVar,i)*s_R(iVar,i)**2+s_L(iVar,i)**2*s_R(iVar,i)) &
                         *(s_L(iVar,i)*FV_sdx(i)-s_R(iVar,i)*FV_sdx(i+1)))/(MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13)**2))
    
        dUdUvolprim_minus(iVar,iVar,i,3) = 0. - FV_dx_L(i) * FV_sdx(i+1) * (&
                       -(2*s_L(iVar,i)*s_R(iVar,i)**2*(s_L(iVar,i)+s_R(iVar,i)))/(MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13)**2) &
                       +(s_L(iVar,i)*(2.*s_R(iVar,i)+s_L(iVar,i)))/MAX(s_L(iVar,i)**2+s_R(iVar,i)**2,1e-13))
      END IF
    END DO !iVar
  CASE(9) ! Central
    ! Central: 0.5*( s_L(U_(i-1),U_i) + s_R(U_(i+1),U_i) )
    DO iVar=1,PP_nVarPrim
      dUdUvolprim_plus( iVar,iVar,i,1) = 0. + FV_dx_R(i) * 0.5*(-FV_sdx(i))
      dUdUvolprim_plus( iVar,iVar,i,2) = 1. + FV_dx_R(i) * 0.5*( FV_sdx(i)-FV_sdx(i+1))
      dUdUvolprim_plus( iVar,iVar,i,3) = 0. + FV_dx_R(i) * 0.5*( FV_sdx(i+1))

      dUdUvolprim_minus(iVar,iVar,i,1) = 0. - FV_dx_L(i) * 0.5*(-FV_sdx(i))
      dUdUvolprim_minus(iVar,iVar,i,2) = 1. - FV_dx_L(i) * 0.5*( FV_sdx(i) - FV_sdx(i+1))
      dUdUvolprim_minus(iVar,iVar,i,3) = 0. - FV_dx_L(i) * 0.5*( FV_sdx(i+1))
    END DO !iVar
  CASE DEFAULT 
    CALL Abort(__STAMP__,'No preconditioner for chosen limiter implemented!')
  END SELECT

  ! multiply: dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim
  DO ind=1,3
    matrix_plus( :,:,i,ind) = MATMUL(Jac_ConsPrim_plus( :,:,i),dUdUvolprim_plus( :,:,i,ind))
    matrix_minus(:,:,i,ind) = MATMUL(Jac_ConsPrim_minus(:,:,i),dUdUvolprim_minus(:,:,i,ind))
  END DO
  ! multiply: (dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim) * dUvol_prim/dUvol
  dUdUvol_plus( :,:,i,1) = MATMUL(matrix_plus(:,:,i,1),Jac_PrimCons(:,:,i-1))
  dUdUvol_plus( :,:,i,2) = MATMUL(matrix_plus(:,:,i,2),Jac_PrimCons(:,:,i))
  dUdUvol_plus( :,:,i,3) = MATMUL(matrix_plus(:,:,i,3),Jac_PrimCons(:,:,i+1))

  dUdUvol_minus(:,:,i,1) = MATMUL(matrix_minus(:,:,i,1),Jac_PrimCons(:,:,i-1))
  dUdUvol_minus(:,:,i,2) = MATMUL(matrix_minus(:,:,i,2),Jac_PrimCons(:,:,i))
  dUdUvol_minus(:,:,i,3) = MATMUL(matrix_minus(:,:,i,3),Jac_PrimCons(:,:,i+1))
END DO !i

END SUBROUTINE FV_Reconstruction_Derivative

!===================================================================================================================================
!> Calculates Jacobian of reconstruction procedure of surface DOFs
!> dF/dUvol = dF/dU_LR * dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim * dUvol_prim/dUvol
!>               |            |                       |                   |
!>    FD in JacSurfInt    dCons/dPrim   derivative of reconstruction  dPrim/dCons
!>                        |_____________________________________________________|
!>                                                  |
!>                                           calculated here
!> If a side is a big mortar side, set the derivative to zero.
!===================================================================================================================================
SUBROUTINE FV_Reconstruction_Derivative_Surf(FV_sdx,FV_dx_L,FV_dx_R,FV_dx_L_nb,FV_dx_R_nb,UPrim_plus,UPrim_minus, &
                                             UPrim_plus_nb,UPrim_minus_nb,URec_extended,Mortar_minus,Mortar_plus, &
                                             dUdUvol_plus,dUdUvol_minus)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_FV_Limiter     ,ONLY: FV_Limiter
USE MOD_FV_Vars        ,ONLY: LimiterType
USE MOD_Jacobian       ,ONLY: dPrimTempdCons,dConsdPrimTemp
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)    :: FV_sdx(0:PP_N+1)                               !< extended inverse of distance between centers of FV subcells
REAL,INTENT(IN)    :: FV_dx_L                                        !< left distance between face and center of FV subcells
REAL,INTENT(IN)    :: FV_dx_R                                        !< right distance between face and center of FV subcells
REAL,INTENT(IN)    :: FV_dx_L_nb                                     !< left dist. between face and center of FV subcells (nbElem)
REAL,INTENT(IN)    :: FV_dx_R_nb                                     !< right dist. between face and center of FV subcells (nbElem)
REAL,INTENT(IN)    :: UPrim_plus(PP_nVarPrim)                        !< primitive reconstructed value at plus side
REAL,INTENT(IN)    :: UPrim_minus(PP_nVarPrim)                       !< primitive reconstructed value at minus side
REAL,INTENT(IN)    :: UPrim_plus_nb(PP_nVarPrim)                     !< primitive reconstructed value at plus side (nbElem)
REAL,INTENT(IN)    :: UPrim_minus_nb(PP_nVarPrim)                    !< primitive reconstructed value at minus side (nbElem)
LOGICAL,INTENT(IN) :: Mortar_minus                                   !< True if minus side is a big mortar side
LOGICAL,INTENT(IN) :: Mortar_plus                                    !< True if plus side is a big mortar side
REAL,INTENT(IN)    :: URec_extended(PP_nVarPrim,-1:PP_N+1)           !< extended rec volume solution
REAL,INTENT(OUT)   :: dUdUvol_plus( PP_nVar,PP_nVar,PP_N:PP_N+1,1:2) !< Jacobian of reconstr. procedure at plus side of element
REAL,INTENT(OUT)   :: dUdUvol_minus(PP_nVar,PP_nVar,-1:0,2:3)        !< Jacobian of reconstr. procedure at minus side of element
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                                :: iVar,ind,i
REAL,DIMENSION(PP_nVarPrim)                            :: s_L_minus,s_R_minus,s_L_plus,s_R_plus,s_lim_minus,s_lim_plus!,s_L_nb,s_R_nb
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,PP_N:PP_N+1)    :: Jac_ConsPrim_plus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,-1:0)           :: Jac_ConsPrim_minus
REAL,DIMENSION(PP_nVarPrim,PP_nVar    , 0:1)           :: Jac_PrimCons_minus
REAL,DIMENSION(PP_nVarPrim,PP_nVar    ,PP_N-1:PP_N)    :: Jac_PrimCons_plus
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,PP_N:PP_N+1,1:2):: dUdUvolprim_plus
REAL,DIMENSION(PP_nVarPrim,PP_nVarPrim,-1:0,2:3)       :: dUdUvolprim_minus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,PP_N:PP_N+1,1:2):: matrix_plus
REAL,DIMENSION(PP_nVar    ,PP_nVarPrim,-1:0,2:3)       :: matrix_minus
!==================================================================================================================================
! storage order in dUdUvol(1:3):
! dU_LR_i/dU_(i-1), dU_LR_i/dU_i, dU_LR_i/dU_(i+1)
dUdUvol_minus = 0.
dUdUvolprim_minus = 0.
dUdUvol_plus  = 0.
dUdUvolprim_plus  = 0.

! 1. Do derivative of surface data from cons to prim
IF (.NOT.Mortar_minus) THEN
  CALL dConsdPrimTemp(UPrim_minus_nb,Jac_ConsPrim_minus(:,:,-1    ))
  CALL dConsdPrimTemp(UPrim_minus   ,Jac_ConsPrim_minus(:,:, 0    ))
END IF
IF (.NOT.Mortar_plus) THEN
CALL dConsdPrimTemp(UPrim_plus    ,Jac_ConsPrim_plus( :,:,PP_N  ))
CALL dConsdPrimTemp(UPrim_plus_nb ,Jac_ConsPrim_plus( :,:,PP_N+1))
END IF

! 2. Do derivative of interface and neighbouring volume data from prim to cons 
IF (.NOT.Mortar_minus) THEN
  DO i=0,1
    CALL dPrimTempdCons(URec_extended(:,i),Jac_PrimCons_minus(:,:,i))
  END DO
END IF
IF (.NOT.Mortar_plus) THEN
  DO i=PP_N-1,PP_N
    CALL dPrimTempdCons(URec_extended(:,i),Jac_PrimCons_plus( :,:,i))
  END DO
END IF

! 3. Calculate derivative of reconstruction
! For NullLimiter and central limiter this is a constant for all variables. This is not the case for the minmod limiter. Hence,
! the Jacobian of the reconstruction is a diagonal matrix with variable entries on the main diagonal.
! Calculate left and right gradients
IF (.NOT.Mortar_minus) THEN
  s_L_minus(:) = (URec_extended(:,0     ) - URec_extended(:,-1    )) * FV_sdx(0     ) 
  s_R_minus(:) = (URec_extended(:,1     ) - URec_extended(:, 0    )) * FV_sdx(1     )
  ! Limit gradients
  CALL FV_Limiter(s_L_minus(:),s_R_minus(:),s_lim_minus(:)) ! only null-, minmod- or central-limiter
END IF
IF (.NOT.Mortar_plus) THEN
  s_L_plus( :) = (URec_extended(:,PP_N  ) - URec_extended(:,PP_N-1)) * FV_sdx(PP_N  ) 
  s_R_plus( :) = (URec_extended(:,PP_N+1) - URec_extended(:,PP_N  )) * FV_sdx(PP_N+1)
  ! Limit gradients
  CALL FV_Limiter(s_L_plus( :),s_R_plus( :),s_lim_plus( :)) ! only null-, minmod- or central-limiter
END IF
SELECT CASE(LimiterType)
CASE(0) ! NullLimiter
  DO iVar=1,PP_nVarPrim
  IF (.NOT.Mortar_minus) dUdUvolprim_minus(iVar,iVar,0   ,2) = 1.
  IF (.NOT.Mortar_plus)  dUdUvolprim_plus( iVar,iVar,PP_N,2) = 1.
  END DO !iVar
CASE(1) ! MinMod
  IF (.NOT.Mortar_minus) THEN
    DO iVar=1,PP_nVarPrim
      ! minus side of element
      ! derivatives of minus value at minus interface
      IF(s_lim_minus(iVar).EQ.0.)THEN ! first order
        dUdUvolprim_minus(iVar,iVar,0,2) = 1.
      ELSEIF(s_lim_minus(iVar).EQ.s_L_minus(iVar))THEN ! use left gradient
        dUdUvolprim_minus(iVar,iVar,0,2) = 1. - FV_sdx(0) * FV_dx_L
      ELSEIF(s_lim_minus(iVar).EQ.s_R_minus(iVar))THEN ! use right gradient
        dUdUvolprim_minus(iVar,iVar,0,2) = 1. + FV_sdx(1) * FV_dx_L
        dUdUvolprim_minus(iVar,iVar,0,3) = 0. - FV_sdx(1) * FV_dx_L
      ELSE
        CALL Abort(__STAMP__,'Slopes do not match with minmod in preconditioner!')
      END IF
      ! The reconstructed value from the left neighbour depends on my own value IF the reconstruction in the left neighbour cell has
      ! been done with the slope between my cell and that neighbour (s_L_minus). Check if that was the case and set the dependency
      ! accordingly. 
      !        reconstr. value      cell mean value        recon. with s_L_minus
      IF(ABS(UPrim_minus_nb(iVar)-(URec_extended(iVar,-1)+s_L_minus(iVar)*FV_dx_R_nb)).LE.1E-12)THEN 
        dUdUvolprim_minus(iVar,iVar,-1,3) = 0. + FV_dx_R_nb * ( FV_sdx(0))
      END IF
    END DO !iVar
  END IF
  IF (.NOT.Mortar_plus) THEN
    DO iVar=1,PP_nVarPrim
      ! plus side of element
      ! derivatives of plus values at plus interface
      IF(s_lim_plus(iVar).EQ.0.)THEN ! first order
        dUdUvolprim_plus(iVar,iVar,PP_N,2) = 1.
      ELSEIF(s_lim_plus(iVar).EQ.s_L_plus(iVar))THEN ! use left gradient
        dUdUvolprim_plus(iVar,iVar,PP_N,1) = 0. - FV_sdx(PP_N  ) * FV_dx_R
        dUdUvolprim_plus(iVar,iVar,PP_N,2) = 1. + FV_sdx(PP_N  ) * FV_dx_R
      ELSEIF(s_lim_plus(iVar).EQ.s_R_plus(iVar))THEN ! use right gradient
        dUdUvolprim_plus(iVar,iVar,PP_N,2) = 1. - FV_sdx(PP_N+1) * FV_dx_R
      ELSE
        CALL Abort(__STAMP__,'Slopes do not match with minmod in preconditioner!')
      END IF
      ! The reconstructed value from the right neighbour depends on my own value IF the reconstruction in the right neighbour cell has
      ! been done with the slope between my cell and that neighbour (s_R_plus). Check if that was the case and set the dependency
      ! accordingly. 
      !        reconstr. value      cell mean value        recon. with s_R_plus
      IF(ABS(UPrim_plus_nb(iVar)-(URec_extended(iVar,PP_N+1)-s_R_plus(iVar)*FV_dx_L_nb)).LE.1E-12)THEN 
        dUdUvolprim_plus( iVar,iVar,PP_N+1,1) = 0. - FV_dx_L_nb * (-FV_sdx(PP_N+1))
      END IF
    END DO !iVar
  END IF
CASE(3) ! van Albada
  IF (.NOT.Mortar_minus) THEN
    DO iVar=1,PP_nVarPrim
      ! minus side of interface
      IF((s_L_minus(iVar)*s_R_minus(iVar)).LT.0)THEN ! zero slope
        IF (.NOT.Mortar_minus) dUdUvolprim_minus(iVar,iVar,0   ,2) = 1.
      ELSE
        ! derivatives of minus value at minus interface
        dUdUvolprim_minus(iVar,iVar,0     ,2) = 1. - FV_dx_L * (&
                         -(s_L_minus(iVar)**2*FV_sdx(1)-s_R_minus(iVar)**2*FV_sdx(0)+2.*s_L_minus(iVar)*s_R_minus(iVar)*FV_sdx(1) &
                          -2.*s_L_minus(iVar)*s_R_minus(iVar)*FV_sdx(0))/MAX(s_L_minus(iVar)**2+s_R_minus(iVar)**2,1e-13) &
                         -(2.*(s_L_minus(iVar)*s_R_minus(iVar)**2+s_L_minus(iVar)**2*s_R_minus(iVar))* &
                           (s_L_minus(iVar)*FV_sdx(0)-s_R_minus(iVar)*FV_sdx(1)))/(MAX(s_L_minus(iVar)**2+s_R_minus(iVar)**2,1e-13)**2))
        dUdUvolprim_minus(iVar,iVar,0     ,3) = 0. - FV_dx_L * FV_sdx(1) * (&
                         -(2*s_L_minus(iVar)*s_R_minus(iVar)**2*(s_L_minus(iVar)+s_R_minus(iVar)))/(MAX(s_L_minus(iVar)**2+s_R_minus(iVar)**2,1e-13)**2) &
                         +(s_L_minus(iVar)*(2.*s_R_minus(iVar)+s_L_minus(iVar)))/MAX(s_L_minus(iVar)**2+s_R_minus(iVar)**2,1e-13))
        ! derivatives of plus value at minus interface
        ! ATTENTION: s_L and s_R are not the same as above!!! s_L_nb is not known as this is from the neighbor!
        !s_L_nb(iVar) = 0.
        !dUdUvolprim_minus(iVar,iVar,-1    ,3) = 0. + FV_dx_R_nb *FV_sdx(0) *(&
                         !-(2*s_L_nb(iVar)*s_L_minus(iVar)**2*(s_L_nb(iVar)+s_L_minus(iVar)))/(MAX(s_L_nb(iVar)**2+s_L_minus(iVar)**2,1e-13)**2) &
                         !+(s_L_nb(iVar)*(2.*s_L_minus(iVar)+s_L_nb(iVar)))/MAX(s_L_nb(iVar)**2+s_L_minus(iVar)**2,1e-13))
      END IF !zero slope
    END DO !iVar
  END IF

  IF (.NOT.Mortar_plus) THEN
    DO iVar=1,PP_nVarPrim
      ! plus side of interface
      IF((s_L_plus(iVar)*s_R_plus(iVar)).LT.0)THEN ! zero slope
        IF (.NOT.Mortar_plus)  dUdUvolprim_plus( iVar,iVar,PP_N,2) = 1.
      ELSE
      ! derivatives of plus values at plus interface
        dUdUvolprim_plus( iVar,iVar,PP_N  ,1) = 0. + FV_dx_R * FV_sdx(PP_N) *(&
                         -(s_R_plus(iVar)*(s_L_plus(iVar)+s_R_plus(iVar))+s_L_plus(iVar)*s_R_plus(iVar))/MAX(s_L_plus(iVar)**2+s_R_plus(iVar)**2,1e-13) &
                         +(2.*s_L_plus(iVar)*(s_L_plus(iVar)**2*s_R_plus(iVar)+s_R_plus(iVar)**2*s_L_plus(iVar)))/(MAX(s_L_plus(iVar)**2+s_R_plus(iVar)**2,1e-13)**2))
        dUdUvolprim_plus( iVar,iVar,PP_N  ,2) = 1. + FV_dx_R * (&
                         -(s_L_plus(iVar)**2*FV_sdx(PP_N+1)-s_R_plus(iVar)**2*FV_sdx(PP_N) & 
                           +2.*s_L_plus(iVar)*s_R_plus(iVar)*FV_sdx(PP_N+1)-2.*s_L_plus(iVar)*s_R_plus(iVar)*FV_sdx(PP_N))/MAX(s_L_plus(iVar)**2+s_R_plus(iVar)**2,1e-13) &
                         -(2.*(s_L_plus(iVar)*s_R_plus(iVar)**2+s_L_plus(iVar)**2*s_R_plus(iVar)) &
                           *(s_L_plus(iVar)*FV_sdx(PP_N)-s_R_plus(iVar)*FV_sdx(PP_N+1)))/(MAX(s_L_plus(iVar)**2+s_R_plus(iVar)**2,1e-13)**2))
        ! derivatives of minus values at plus interface
        ! ATTENTION: s_L and s_R are not the same as above!!! s_R_nb is not known as this is from the neighbor!
        !s_R_nb(iVar) = 0.
        !dUdUvolprim_plus( iVar,iVar,PP_N+1,1) = 0. - FV_dx_L_nb * FV_sdx(PP_N+1) * (&
                         !-(s_R_nb(iVar)*(s_R_plus(iVar)+s_R_nb(iVar))+s_R_plus(iVar)*s_R_nb(iVar))/MAX(s_R_plus(iVar)**2+s_R_nb(iVar)**2,1e-13) &
                         !+(2.*s_R_plus(iVar)*(s_R_plus(iVar)**2*s_R_nb(iVar)+s_R_nb(iVar)**2*s_R_plus(iVar)))/(MAX(s_R_plus(iVar)**2+s_R_nb(iVar)**2,1e-13)**2))
      END IF !zero slope
    END DO !iVar
  END IF
CASE(9) ! Central
  IF (.NOT.Mortar_minus) THEN
    DO iVar=1,PP_nVarPrim
      ! minus side of interface
      ! derivatives of minus value at minus interface
      dUdUvolprim_minus(iVar,iVar,0     ,2) = 1. - FV_dx_L    * 0.5*( FV_sdx(0) - FV_sdx(1))
      dUdUvolprim_minus(iVar,iVar,0     ,3) = 0. - FV_dx_L    * 0.5*( FV_sdx(1))
      ! derivatives of plus value at minus interface
      dUdUvolprim_minus(iVar,iVar,-1    ,3) = 0. + FV_dx_R_nb * 0.5*( FV_sdx(0))
    END DO !iVar
  END IF

  IF (.NOT.Mortar_plus) THEN
    DO iVar=1,PP_nVarPrim
      ! plus side of interface
      ! derivatives of plus values at plus interface
      dUdUvolprim_plus( iVar,iVar,PP_N  ,1) = 0. + FV_dx_R    * 0.5*(-FV_sdx(PP_N))
      dUdUvolprim_plus( iVar,iVar,PP_N  ,2) = 1. + FV_dx_R    * 0.5*( FV_sdx(PP_N)-FV_sdx(PP_N+1))
      ! derivatives of minus values at plus interface
      dUdUvolprim_plus( iVar,iVar,PP_N+1,1) = 0. - FV_dx_L_nb * 0.5*(-FV_sdx(PP_N+1))
    END DO !iVar
  END IF
CASE DEFAULT 
  CALL Abort(__STAMP__,'No preconditioner for chosen limiter implemented!')
END SELECT

! multiply: dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim
IF (.NOT.Mortar_minus) THEN
  DO ind=2,3
    matrix_minus(:,:,0,ind) = MATMUL(Jac_ConsPrim_minus(:,:,0),dUdUvolprim_minus(:,:,0,ind))
  END DO
  matrix_minus(:,:,-1,3) = MATMUL(Jac_ConsPrim_minus(:,:,-1),dUdUvolprim_minus(:,:,-1,3))
  ! multiply: (dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim) * dUvol_prim/dUvol
  dUdUvol_minus(:,:,-1,3) = MATMUL(matrix_minus(:,:,-1,3),Jac_PrimCons_minus(:,:,0))
  dUdUvol_minus(:,:,0,2)  = MATMUL(matrix_minus(:,:,0,2),Jac_PrimCons_minus(:,:,0))
  dUdUvol_minus(:,:,0,3)  = MATMUL(matrix_minus(:,:,0,3),Jac_PrimCons_minus(:,:,1))
END IF
IF (.NOT.Mortar_plus) THEN
  DO ind=1,2
    matrix_plus( :,:,PP_N,ind) = MATMUL(Jac_ConsPrim_plus( :,:,PP_N),dUdUvolprim_plus( :,:,PP_N,ind))
  END DO
  matrix_plus( :,:,PP_N+1,1) = MATMUL(Jac_ConsPrim_plus( :,:,PP_N+1),dUdUvolprim_plus( :,:,PP_N+1,1))
  ! multiply: (dU_LR/dU_LR_prim * dU_LR_prim/dUvol_prim) * dUvol_prim/dUvol
  dUdUvol_plus( :,:,PP_N,1)   = MATMUL(matrix_plus(:,:,PP_N,1),Jac_PrimCons_plus(:,:,PP_N-1))
  dUdUvol_plus( :,:,PP_N,2)   = MATMUL(matrix_plus(:,:,PP_N,2),Jac_PrimCons_plus(:,:,PP_N))
  dUdUvol_plus( :,:,PP_N+1,1) = MATMUL(matrix_plus(:,:,PP_N+1,1),Jac_PrimCons_plus(:,:,PP_N))
END IF

END SUBROUTINE FV_Reconstruction_Derivative_Surf

#if PARABOLIC
!===================================================================================================================================
!> Computes the Volume gradient Jacobian of the reconstruction procedure dQprim/dUprim (Q= Grad U) only w.r.t. the volume DOFs of
!> the current element.
!> As the reconstruction procedure is the same for all primitive variables and is done independently for each variable the Jacobian
!> has diagonal shape with the same coefficient on the diagonal. Therefore, we calculate only this scalar value and give it back.
!> The reconstruction procedure for the gradients takes the central limiter. Therefore the gradient at a DOF (i,j,k) is depending on
!> the left and right neighbour and the DOF itself.
!> If side is a big mortar side set the Jacobian of the evolved gradient to zero.
!===================================================================================================================================
SUBROUTINE JacFVGradients_Vol(dir,iElem,Jac_reconstruct)
! MODULES
USE MOD_PreProc
USE MOD_FV_Vars            ,ONLY: FV_sdx_XI,FV_sdx_ETA,FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ
USE MOD_Jac_Ex_Vars        ,ONLY: FV_sdx_XI_extended,FV_sdx_ETA_extended
#if PP_dim == 3     
USE MOD_FV_Vars            ,ONLY: FV_sdx_ZETA,FV_Metrics_hTilde_sJ
USE MOD_Jac_Ex_Vars        ,ONLY: FV_sdx_ZETA_extended
#endif
USE MOD_Mesh_Vars          ,ONLY: ElemToSide,MortarType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem ! < current element index
INTEGER,INTENT(IN) :: dir   ! < current physical direction
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Jac_reconstruct(0:PP_N,0:PP_N,0:PP_NZ,0:PP_N,PP_dim) !< Jacobian of volume gradients in direction dir
                                                                           !> w.r.t. primitive volume solution
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,ll,Mortar_minus,Mortar_plus
!===================================================================================================================================
Jac_reconstruct = 0.
DO k=0,PP_NZ
  DO j=0,PP_N
    DO i=0,PP_N
      ! Contribution by gradient reconstruction procedure with central gradient (i,j,k) -> gradient, (ll) -> state 
      ! --------- XI direction ------------ !
      Mortar_minus = MortarType(1,ElemToSide(E2S_SIDE_ID,XI_MINUS,iElem))
      Mortar_plus  = MortarType(1,ElemToSide(E2S_SIDE_ID,XI_PLUS ,iElem))
      IF(i.GT.0)THEN
        Jac_reconstruct(i,j,k,i-1,1) = 0.5*(-FV_sdx_XI(j,k,i,iElem))
      END IF
      IF(.NOT.(((i.EQ.0).AND.(Mortar_minus.GT.0)).OR.((i.EQ.PP_N).AND.(Mortar_plus.GT.0))))THEN
        Jac_reconstruct(i,j,k,i  ,1) = 0.5*FV_sdx_XI_extended(j,k,i,iElem) - 0.5*FV_sdx_XI_extended(j,k,i+1,iElem)
      END IF
      IF(i.LT.PP_N)THEN          
        Jac_reconstruct(i,j,k,i+1,1) = 0.5*( FV_sdx_XI(j,k,i+1,iElem))
      END IF
      ! --------- ETA direction ------------ !
      Mortar_minus = MortarType(1,ElemToSide(E2S_SIDE_ID,ETA_MINUS,iElem))
      Mortar_plus  = MortarType(1,ElemToSide(E2S_SIDE_ID,ETA_PLUS ,iElem))
      IF(j.GT.0)THEN
        Jac_reconstruct(i,j,k,j-1,2) = 0.5*(-FV_sdx_ETA(i,k,j,iElem))
      END IF                     
      IF(.NOT.(((j.EQ.0).AND.(Mortar_minus.GT.0)).OR.((j.EQ.PP_N).AND.(Mortar_plus.GT.0))))THEN
        Jac_reconstruct(i,j,k,j  ,2) = 0.5*FV_sdx_ETA_extended(i,k,j,iElem) - 0.5*FV_sdx_ETA_extended(i,k,j+1,iElem)
      END IF
      IF(j.LT.PP_N)THEN          
        Jac_reconstruct(i,j,k,j+1,2) = 0.5*( FV_sdx_ETA(i,k,j+1,iElem))
      END IF
#if PP_dim==3
      ! --------- ZETA direction ------------ !
      Mortar_minus = MortarType(1,ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem))
      Mortar_plus  = MortarType(1,ElemToSide(E2S_SIDE_ID,ZETA_PLUS ,iElem))
      IF(k.GT.0)THEN
        Jac_reconstruct(i,j,k,k-1,3) = 0.5*(-FV_sdx_ZETA(i,j,k,iElem))
      END IF                     
      IF(.NOT.(((k.EQ.0).AND.(Mortar_minus.GT.0)).OR.((k.EQ.PP_N).AND.(Mortar_plus.GT.0))))THEN
        Jac_reconstruct(i,j,k,k  ,3) = 0.5*FV_sdx_ZETA_extended(i,j,k,iElem) - 0.5*FV_sdx_ZETA_extended(i,j,k+1,iElem)
      END IF
      IF(k.LT.PP_N)THEN          
        Jac_reconstruct(i,j,k,k+1,3) = 0.5*( FV_sdx_ZETA(i,j,k+1,iElem))
      END IF
#endif
      ! Contribution by lifting volume integral (only apply metrics for FV elements)
      DO ll=0,PP_N
        Jac_reconstruct(i,j,k,ll,1) = Jac_reconstruct(i,j,k,ll,1)*FV_Metrics_fTilde_sJ(dir,ll,j,k,iElem)
        Jac_reconstruct(i,j,k,ll,2) = Jac_reconstruct(i,j,k,ll,2)*FV_Metrics_gTilde_sJ(dir,i,ll,k,iElem)
#if PP_dim==3
        Jac_reconstruct(i,j,k,ll,3) = Jac_reconstruct(i,j,k,ll,3)*FV_Metrics_hTilde_sJ(dir,i,j,ll,iElem)
#endif
      END DO ! ll
    END DO !i
  END DO !j
END DO !k
END SUBROUTINE JacFVGradients_Vol

!===================================================================================================================================
!> Computes the dependency of FV gradient of neighbouring element at the interface on the volume DOFs of the considered element.
!> Again, as the reconstruction of the gradients is independent of the variable (iVar) itself, the output of this routine is a
!> scalar.
!> For the second order reconstruction with central limiter the gradient of the neighbouring element only depends on the outer layer
!> of the current element.
!> |...:...:...: g | x :...:...: x | g :...:...:...| => gradients g only depend on the outer layer DOF x of the current element
!> |...:...:...:...|...:...:...:...|...:...:...:...|
!>    nbElem_minus       iElem        nbElem_plus
!===================================================================================================================================
SUBROUTINE JacFVGradients_nb(dir,iElem,dQ_dUVolOuter)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars                 ,ONLY: ElemToSide,nBCSides
#if FV_ENABLED
USE MOD_Jac_Ex_Vars               ,ONLY: FV_sdx_XI_extended,FV_sdx_ETA_extended
USE MOD_FV_Vars                   ,ONLY: FV_Metrics_fTilde_sJ,FV_Metrics_gTilde_sJ
#if PP_dim == 3
USE MOD_Jac_Ex_Vars               ,ONLY: FV_sdx_ZETA_extended
USE MOD_FV_Vars                   ,ONLY: FV_Metrics_hTilde_sJ
#endif
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: dir   !< considered direction (1,2,3)
INTEGER,INTENT(IN) :: iElem !< considered element ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
#if PP_dim == 3
REAL,INTENT(OUT) :: dQ_dUVolOuter(0:PP_N,0:PP_NZ,6,0:PP_N)!< Jacobian of surface gradients of the neighbour element w.r.t. primitive
                                                          !> solution of current element
#else
REAL,INTENT(OUT) :: dQ_dUVolOuter(0:PP_N,0:PP_NZ,2:5,0:PP_N)
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iLocSide,SideID,i,j,k
!===================================================================================================================================
dQ_dUVolOuter=0.

#if PP_dim == 3
DO iLocSide=1,6
#else    
DO iLocSide=2,5
#endif    
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
  IF (SideID.LE.nBCSides) CYCLE  !for boundary conditions, dQ_dUVolOuter=0.

  SELECT CASE(iLocSide)
  CASE(XI_MINUS,XI_PLUS)
    DO k=0,PP_NZ
      DO j=0,PP_N
        IF(iLocSide.EQ.XI_MINUS)THEN
          dQ_dUVolOuter(j,k,XI_MINUS,   0) = 0.5*FV_sdx_XI_extended(j,k,     0,iElem) * FV_Metrics_fTilde_sJ(dir,   0,j,k,iElem)
        ELSE
          dQ_dUVolOuter(j,k,XI_PLUS ,PP_N) = 0.5*FV_sdx_XI_extended(j,k,PP_N+1,iElem) * FV_Metrics_fTilde_sJ(dir,PP_N,j,k,iElem)
        END IF
      END DO !j
    END DO !k
  CASE(ETA_MINUS,ETA_PLUS)
    DO k=0,PP_NZ
      DO i=0,PP_N
        IF(iLocSide.EQ.ETA_MINUS)THEN
          dQ_dUVolOuter(i,k,ETA_MINUS,   0) = 0.5*FV_sdx_ETA_extended(i,k,     0,iElem) * FV_Metrics_gTilde_sJ(dir,i,   0,k,iElem)
        ELSE
          dQ_dUVolOuter(i,k,ETA_PLUS ,PP_N) = 0.5*FV_sdx_ETA_extended(i,k,PP_N+1,iElem) * FV_Metrics_gTilde_sJ(dir,i,PP_N,k,iElem)
        END IF
      END DO !i
    END DO !k
#if PP_dim==3
  CASE(ZETA_MINUS,ZETA_PLUS)
    DO j=0,PP_N
      DO i=0,PP_N
        IF(iLocSide.EQ.ZETA_MINUS)THEN
          dQ_dUVolOuter(i,j,ZETA_MINUS,   0) = 0.5*FV_sdx_ZETA_extended(i,j,     0,iElem) * FV_Metrics_hTilde_sJ(dir,i,j,   0,iElem)
        ELSE
          dQ_dUVolOuter(i,j,ZETA_PLUS ,PP_N) = 0.5*FV_sdx_ZETA_extended(i,j,PP_N+1,iElem) * FV_Metrics_hTilde_sJ(dir,i,j,PP_N,iElem)
        END IF
      END DO !i
    END DO !j
#endif
  END SELECT
END DO !iLocSide
END SUBROUTINE JacFVGradients_nb
#endif /*PARABOLIC*/
#endif /*FV_ENABLED && FV_RECONSTRUCT*/

END MODULE MOD_Jac_Ex_Reconstruction
