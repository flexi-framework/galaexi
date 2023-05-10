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

!==================================================================================================================================
!> \brief Prepares the integrand for the flux integrals, i.e. computes the common inteface fluxes.
!> This module prepares the computation of the fluxes over the sides by filling the two parts (advective and viscous)
!> of the common interface flux by calling the associated numerical flux functions.
!> We distinguish between inner sides, sides with boundary conditions and mpi sides.
!==================================================================================================================================
MODULE MOD_FillFlux
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------


INTERFACE FillFlux
  MODULE PROCEDURE FillFlux
END INTERFACE


PUBLIC::FillFlux
!==================================================================================================================================



CONTAINS



!==================================================================================================================================
!> Computes the fluxes for inner sides, MPI sides where the local proc is "master"  and boundary conditions.
!> The flux computation is performed separately for advection and diffusion fluxes in case
!> parabolic terms are considered.
!==================================================================================================================================
SUBROUTINE FillFlux(t,d_Flux_master,d_Flux_slave,d_U_master,d_U_slave,d_UPrim_master,d_UPrim_slave,doMPISides)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_DG_Vars,         ONLY: nDOFFace
USE MOD_Mesh_Vars,       ONLY: d_NormVec, d_TangVec1, d_TangVec2, d_SurfElem, Face_xGP
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars,       ONLY: nSides,firstBCSide
USE MOD_ChangeBasisByDim,ONLY: ChangeBasisSurf
USE MOD_Riemann,         ONLY: Riemann,Riemann_CPU
USE MOD_GetBoundaryFlux, ONLY: GetBoundaryFlux
USE MOD_EOS,             ONLY: ConsToPrim
USE MOD_Mesh_Vars,       ONLY: nBCSides
USE MOD_EOS,             ONLY: PrimToCons
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN) :: doMPISides  !< = .TRUE. only MINE (where the proc is master)  MPISides are filled, =.FALSE. InnerSides
REAL,INTENT(IN)    :: t           !< physical time required for BC state evaluation in case of time dependent BCs
REAL,INTENT(OUT)   :: d_Flux_master(1:PP_nVar,0:PP_N,0:PP_NZ,1:nSides)      !< sum of advection and diffusion fluxes across the boundary
REAL,INTENT(OUT)   :: d_Flux_slave (1:PP_nVar,0:PP_N,0:PP_NZ,1:nSides)      !< sum of advection and diffusion fluxes across the boundary
REAL,INTENT(INOUT) :: d_U_master(    PP_nVar,0:PP_N, 0:PP_NZ,1:nSides)      !< solution on master sides
REAL,INTENT(INOUT) :: d_U_slave(     PP_nVar,0:PP_N, 0:PP_NZ,1:nSides)      !< solution on slave sides
REAL,INTENT(IN)    :: d_UPrim_master(PP_nVarPrim,0:PP_N, 0:PP_NZ, 1:nSides) !< primitive solution on master sides
REAL,INTENT(IN)    :: d_UPrim_slave( PP_nVarPrim,0:PP_N, 0:PP_NZ, 1:nSides) !< primitive solution on slave sides
!@cuf ATTRIBUTES(DEVICE) :: d_U_master,d_U_slave,d_UPrim_master,d_UPrim_slave,d_Flux_master,d_Flux_slave
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER  :: nBlockSides=128
INTEGER            :: firstBlockSide,lastBlockSide,nMyBlockSides
INTEGER :: SideID,p,q,firstSideID_wo_BC,firstSideID ,lastSideID,FVEM
INTEGER :: FV_Elems_Max(1:nSides) ! 0 if both sides DG, 1 else
!==================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID using Riemann solver for advection and viscous terms
! Set the side range according to MPI or no MPI
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides (where the local proc is master)
  firstSideID_wo_BC = firstMPISide_MINE
  firstSideID = firstMPISide_MINE
   lastSideID =  lastMPISide_MINE
ELSE
  ! fill only InnerSides that do not need communication
  firstSideID_wo_BC = firstInnerSide ! for fluxes
  firstSideID = firstBCSide    ! include BCs for master sides
   lastSideID = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  FV_Elems_Max(SideID) = MAX(FV_Elems_master(SideID),FV_Elems_slave(SideID))
END DO

! =============================
! Workflow:
!
!  1.  compute flux for non-BC sides
!  1.1) advective flux
!  1.2) viscous flux
!  1.3) add up viscous flux to Flux_master
!  2.  compute flux for BC sides
!  3.  multiply by SurfElem
!  4.  copy flux from Flux_master to Flux_slave
!  5.  convert FV flux to DG flux at mixed interfaces
!==============================

! 1. compute flux for non-BC sides
DO firstBlockSide=firstSideID_wo_BC,lastSideID,nBlockSides
  lastBlockSide = MIN(lastSideID,firstBlockSide+nBlockSides-1)
  nMyBlockSides = lastBlockSide-firstBlockSide+1
  ! 1.1) advective part of flux
  !CALL Riemann_CPU(PP_N,Flux_master(:,:,:,SideID),&
  !    U_master    (:,:,:,SideID),U_slave    (:,:,:,SideID),       &
  !    UPrim_master(:,:,:,SideID),UPrim_slave(:,:,:,SideID),       &
  !    NormVec (:,:,:,FV_Elems_Max(SideID),SideID), &
  !    TangVec1(:,:,:,FV_Elems_Max(SideID),SideID), &
  !    TangVec2(:,:,:,FV_Elems_Max(SideID),SideID),doBC=.FALSE.)
  CALL Riemann<<<(nDOFFace*nMyBlockSides/256+1),256>>>(   &
      (nDOFFace*nMyBlockSides),                           &
      d_Flux_master( :,:,:,firstBlockSide:lastBlockSide), &
      d_U_master(    :,:,:,firstBlockSide:lastBlockSide), &
      d_U_slave(     :,:,:,firstBlockSide:lastBlockSide), &
      d_UPrim_master(:,:,:,firstBlockSide:lastBlockSide), &
      d_UPrim_slave( :,:,:,firstBlockSide:lastBlockSide), &
      d_NormVec (  :,:,:,:,firstBlockSide:lastBlockSide), &
      d_TangVec1(  :,:,:,:,firstBlockSide:lastBlockSide), &
      d_TangVec2(  :,:,:,:,firstBlockSide:lastBlockSide))!,doBC=.FALSE.)
END DO ! SideID

!Flux_master(:,:,:,firstSideID_wo_BC:lastSideID) = d_Flux_master(:,:,:,firstSideID_wo_BC:lastSideID)

! 2. Compute the fluxes at the boundary conditions: 1..nBCSides
IF(.NOT.doMPISides)THEN
  DO SideID=1,nBCSides
    CALL ABORT(__STAMP__,"No boundaries currently supported!")
!    FVEM = FV_Elems_master(SideID)
!    CALL GetBoundaryFlux(SideID,t,PP_N,&
!       Flux_master(  :,:,:,     SideID),&
!       UPrim_master( :,:,:,     SideID),&
!#if PARABOLIC
!       gradUx_master(:,:,:,     SideID),&
!       gradUy_master(:,:,:,     SideID),&
!       gradUz_master(:,:,:,     SideID),&
!#endif
!       NormVec(      :,:,:,FVEM,SideID),&
!       TangVec1(     :,:,:,FVEM,SideID),&
!       TangVec2(     :,:,:,FVEM,SideID),&
!       Face_xGP(     :,:,:,FVEM,SideID))
  END DO
END IF ! .NOT. MPISIDES


! 3. multiply by SurfElem
!$cuf kernel do(3) <<< *, 256 >>>
DO SideID=firstSideID,lastSideID
  ! multiply with SurfElem
  DO q=0,PP_NZ; DO p=0,PP_N
    d_Flux_master(:,p,q,SideID) = d_Flux_master(:,p,q,SideID) * d_SurfElem(p,q,0,SideID)
    d_Flux_slave( :,p,q,SideID) = d_Flux_master(:,p,q,SideID)
  END DO; END DO
END DO ! SideID

!! 4. copy flux from master side to slave side
!d_Flux_slave(:,:,:,firstSideID:lastSideID) = d_Flux_master(:,:,:,firstSideID:lastSideID)


END SUBROUTINE FillFlux



END MODULE MOD_FillFlux
