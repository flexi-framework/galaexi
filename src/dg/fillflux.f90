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
USE MOD_DG_Vars,         ONLY: nDOFFace,Flux_master
USE MOD_Mesh_Vars,       ONLY: d_NormVec, d_TangVec1, d_TangVec2, d_SurfElem, Face_xGP
USE MOD_Mesh_Vars,       ONLY: firstInnerSide,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE
USE MOD_Mesh_Vars,       ONLY: nSides,firstBCSide
USE MOD_Riemann,         ONLY: Riemann
USE MOD_GetBoundaryFlux, ONLY: GetBoundaryFlux
USE MOD_Mesh_Vars,       ONLY: nBCSides
#if PARABOLIC
USE MOD_Riemann,         ONLY: ViscousFlux
USE MOD_Lifting_Vars,    ONLY: d_gradUx_master,d_gradUy_master,d_gradUz_master
USE MOD_Lifting_Vars,    ONLY: d_gradUx_slave ,d_gradUy_slave ,d_gradUz_slave
#endif
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
INTEGER :: firstSideID_wo_BC,firstSideID,lastSideID
INTEGER :: SideID,p,q
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

! =============================
! Workflow:
!
!  1.  compute flux for non-BC sides
!  1.1) advective flux
!  1.2) viscous flux
!  1.3) add up viscous flux to Flux_master
!  2.  compute flux for BC sides
!  3.  multiply by SurfElem and copy flux from master side to slave side
!==============================
! 1.  compute flux for non-BC sides
IF (firstSideID_wo_BC.LE.lastSideID) THEN
  ! 1.1) advective flux
  CALL Riemann(PP_N, &
               lastSideID-firstSideID_wo_BC+1, & ! Number of sides in array
               d_Flux_master (:,:,:,firstSideID_wo_BC:lastSideID), &
               d_U_master    (:,:,:,firstSideID_wo_BC:lastSideID), &
               d_U_slave     (:,:,:,firstSideID_wo_BC:lastSideID), &
               d_UPrim_master(:,:,:,firstSideID_wo_BC:lastSideID), &
               d_UPrim_slave (:,:,:,firstSideID_wo_BC:lastSideID), &
               d_NormVec   (:,:,:,:,firstSideID_wo_BC:lastSideID), &
               d_TangVec1  (:,:,:,:,firstSideID_wo_BC:lastSideID), &
               d_TangVec2  (:,:,:,:,firstSideID_wo_BC:lastSideID))
#if PARABOLIC
  ! 1.2) viscous flux
  CALL ViscousFlux(PP_N, &
                   lastSideID-firstSideID_wo_BC+1, & ! Number of sides in array
                   d_Flux_master  (:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_UPrim_master (:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_UPrim_slave  (:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_gradUx_master(:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_gradUx_slave (:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_gradUy_master(:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_gradUy_slave (:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_gradUz_master(:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_gradUz_slave (:,:,:,firstSideID_wo_BC:lastSideID), &
                   d_NormVec    (:,:,:,:,firstSideID_wo_BC:lastSideID))
#endif
END IF

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


! 3. multiply by SurfElem and copy flux from master side to slave side
!$cuf kernel do(3) <<< *, 256, 0 >>>
DO SideID=firstSideID,lastSideID
  DO q=0,PP_NZ; DO p=0,PP_N
    ! multiply with SurfElem
    d_Flux_master(:,p,q,SideID) = d_Flux_master(:,p,q,SideID) * d_SurfElem(p,q,0,SideID)
    ! copy master to slave
    d_Flux_slave( :,p,q,SideID) = d_Flux_master(:,p,q,SideID)
  END DO; END DO
END DO ! SideID

END SUBROUTINE FillFlux

END MODULE MOD_FillFlux
