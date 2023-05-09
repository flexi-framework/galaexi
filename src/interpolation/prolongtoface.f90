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
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Contains routines to interpolate the interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFace
IMPLICIT NONE
PRIVATE
#define WITHnVar 1

INTERFACE ProlongToFace
  MODULE PROCEDURE ProlongToFace
END INTERFACE

INTERFACE EvalElemFace
  MODULE PROCEDURE EvalElemFaceG
  MODULE PROCEDURE EvalElemFaceGL
END INTERFACE

PUBLIC::ProlongToFace,EvalElemFace

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFace

!==================================================================================================================================
!> Contains routines to interpolate the conservative interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFaceCons
IMPLICIT NONE
PRIVATE
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE ProlongToFaceCons
  MODULE PROCEDURE ProlongToFace
END INTERFACE

INTERFACE ProlongToFaceCons_GPU
  MODULE PROCEDURE ProlongToFace_GPU
END INTERFACE

PUBLIC::ProlongToFaceCons,ProlongToFaceCons_GPU

CONTAINS
#include "prolongtoface.t90"

! Avoid compiling it for all TP_nVar to improve compilation speed

!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
SUBROUTINE ProlongToFace_GPU(&
#ifdef WITHnVar
    TP_nVar,&
#endif
    Nloc,Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,doMPISides &
)
! MODULES
USE MOD_Mesh_Vars,          ONLY: nElems,nSides
!@cuf USE MOD_Mesh_Vars,             ONLY: d_SideToElem,d_S2V2,d_ElemToSide
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
#ifdef WITHnVar
INTEGER,INTENT(IN)              :: TP_nVar
#endif
INTEGER,INTENT(IN)              :: Nloc
REAL,DEVICE,INTENT(IN)          :: Uvol(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)
REAL,DEVICE,INTENT(INOUT)       :: Uface_master(TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,DEVICE,INTENT(INOUT)       :: Uface_slave( TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,DEVICE,INTENT(IN)          :: L_Minus(0:Nloc),L_Plus(0:Nloc)
LOGICAL,INTENT(IN)              :: doMPISides  != .TRUE. only YOUR MPISides are filled, =.FALSE. BCSides +InnerSides +MPISides MINE
#if FV_ENABLED
LOGICAL,INTENT(IN),OPTIONAL     :: pureDG      != .TRUE. prolongates all elements as DG elements
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DEVICE                     :: Uface_work(TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides,2)
INTEGER,PARAMETER               :: nThreads=12
!==================================================================================================================================

! CALL ProlongToFace_Kernel<<<nSides/nThreads+1,nThreads>>>(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,d_SideToElem,d_S2V2)
! CALL ProlongToFace_Kernel_Elem<<<nElems/nThreads+1,nThreads>>>(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,d_ElemToSide,d_S2V2)
! CALL ProlongToFace_Kernel_Elem_locSide<<<(nElems*6)/nThreads+1,nThreads>>>(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,d_ElemToSide,d_S2V2)
CALL ProlongToFace_Kernel_Elem_locSide_PreAlloc<<<(nElems*6)/nThreads+1,nThreads>>>(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,Uface_work,L_Minus,L_Plus,d_ElemToSide,d_S2V2)
END SUBROUTINE ProlongToFace_GPU



!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE ProlongToFace_Kernel(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,SideToElem,S2V2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
#ifdef WITHnVar
INTEGER,VALUE,INTENT(IN)        :: TP_nVar
#endif
INTEGER,VALUE,INTENT(IN)        :: Nloc
INTEGER,VALUE,INTENT(IN)        :: nSides
INTEGER,VALUE,INTENT(IN)        :: nElems
REAL,INTENT(IN)                 :: Uvol(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)
REAL,INTENT(INOUT)              :: Uface_master(TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave( TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: p,q
INTEGER                         :: ElemID,nbElemID,locSide,nblocSide,SideID,flip
REAL                            :: Uface(TP_nVar,0:Nloc,0:ZDIM(Nloc))
!==================================================================================================================================
SideID = (blockidx%x-1) * blockdim%x + threadidx%x
IF (SideID.LE.nSides) THEN
  ElemID    = SideToElem(S2E_ELEM_ID,SideID)
  nbElemID  = SideToElem(S2E_NB_ELEM_ID,SideID)

  !master sides
  IF(ElemID.GT.0)THEN
    locSide = SideToElem(S2E_LOC_SIDE_ID,SideID)
    flip    = 0

    CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)

    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      Uface_master(:,p,q,SideID)=Uface(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
    END DO; END DO
  END IF

  !slave side (ElemID,locSide and flip =-1 if not existing)
  IF(nbElemID.GT.0)THEN
   nblocSide = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
   flip      = SideToElem(S2E_FLIP,SideID)

    CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,nbElemID),Uface,L_Minus,L_Plus,nblocSide)

    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      Uface_slave( :,p,q,SideID)=Uface(:,S2V2(1,p,q,flip,nblocSide),S2V2(2,p,q,flip,nblocSide))
    END DO; END DO
  END IF
END IF

END SUBROUTINE ProlongToFace_Kernel


!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE ProlongToFace_Kernel_Elem(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,L_Minus,L_Plus,ElemToSide,S2V2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
#ifdef WITHnVar
INTEGER,VALUE,INTENT(IN)        :: TP_nVar
#endif
INTEGER,VALUE,INTENT(IN)        :: Nloc
INTEGER,VALUE,INTENT(IN)        :: nSides
INTEGER,VALUE,INTENT(IN)        :: nElems
REAL,INTENT(IN)                 :: Uvol(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)
REAL,INTENT(INOUT)              :: Uface_master(TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave( TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
INTEGER,INTENT(IN)              :: ElemToSide(3,6,nElems)
!INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: p,q
INTEGER                         :: ElemID,nbElemID,locSide,SideID,flip
REAL                            :: Uface(TP_nVar,0:Nloc,0:ZDIM(Nloc))
LOGICAL                         :: isMaster
!==================================================================================================================================
ElemID = (blockidx%x-1) * blockdim%x + threadidx%x
IF (ElemID.LE.nElems) THEN
  DO locSide=1,6
    SideID   = ElemToSide(E2S_SIDE_ID  ,locSide,ElemID)
    flip     = ElemToSide(E2S_FLIP     ,locSide,ElemID)
    isMaster = ElemToSide(E2S_IS_MASTER,locSide,ElemID).EQ.1 ! master side for current elem


    IF(isMaster) THEN
      CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        Uface_master(:,p,q,SideID)=Uface(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
      END DO; END DO
    ELSE ! slave
      CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        Uface_slave( :,p,q,SideID)=Uface(:,S2V2(1,p,q,flip,locSide),S2V2(2,p,q,flip,locSide))
      END DO; END DO
    END IF
  END DO !locSideID
END IF

END SUBROUTINE ProlongToFace_Kernel_Elem

!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE ProlongToFace_Kernel_Elem_locSide(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,Uface,L_Minus,L_Plus,ElemToSide,S2V2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
#ifdef WITHnVar
INTEGER,VALUE,INTENT(IN)        :: TP_nVar
#endif
INTEGER,VALUE,INTENT(IN)        :: Nloc
INTEGER,VALUE,INTENT(IN)        :: nSides
INTEGER,VALUE,INTENT(IN)        :: nElems
REAL,INTENT(IN)                 :: Uvol(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)
REAL,INTENT(INOUT)              :: Uface_master(TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave( TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
INTEGER,INTENT(IN)              :: ElemToSide(3,6,nElems)
!INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: p,q
INTEGER                         :: ElemID,nbElemID,locSide,SideID,flip,threadID
REAL                            :: Uface(TP_nVar,0:Nloc,0:ZDIM(Nloc))
LOGICAL                         :: isMaster
!==================================================================================================================================
threadID = (blockidx%x-1) * blockdim%x + threadidx%x
ElemID  = (threadID-1)/6 + 1
locSide = threadID-(ElemID-1)*6
IF (ElemID.LE.nElems) THEN
  !DO locSide=1,6
    SideID   = ElemToSide(E2S_SIDE_ID  ,locSide,ElemID)
    flip     = ElemToSide(E2S_FLIP     ,locSide,ElemID)
    isMaster = ElemToSide(E2S_IS_MASTER,locSide,ElemID).EQ.1 ! master side for current elem

    CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)

    IF(isMaster) THEN
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        Uface_master(:,p,q,SideID)=Uface(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
      END DO; END DO
    ELSE ! slave
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        Uface_slave( :,p,q,SideID)=Uface(:,S2V2(1,p,q,flip,locSide),S2V2(2,p,q,flip,locSide))
      END DO; END DO
    END IF
  !END DO !locSideID
END IF

END SUBROUTINE ProlongToFace_Kernel_Elem_locSide


!==================================================================================================================================
!> Interpolates the interior volume data (stored at the Gauss or Gauss-Lobatto points) to the surface
!> integration points, using fast 1D Interpolation and store in global side structure
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE ProlongToFace_Kernel_Elem_locSide_PreAlloc(Nloc,nSides,nElems,Uvol,Uface_master,Uface_slave,Uface,L_Minus,L_Plus,ElemToSide,S2V2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
#ifdef WITHnVar
INTEGER,VALUE,INTENT(IN)        :: TP_nVar
#endif
INTEGER,VALUE,INTENT(IN)        :: Nloc
INTEGER,VALUE,INTENT(IN)        :: nSides
INTEGER,VALUE,INTENT(IN)        :: nElems
REAL,INTENT(IN)                 :: Uvol(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)
REAL,INTENT(INOUT)              :: Uface_master(TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(INOUT)              :: Uface_slave( TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides)
REAL,INTENT(OUT)                :: Uface(TP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides,2)
REAL,INTENT(IN)                 :: L_Minus(0:Nloc),L_Plus(0:Nloc)
INTEGER,INTENT(IN)              :: ElemToSide(3,6,nElems)
!INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: p,q
INTEGER                         :: ElemID,nbElemID,locSide,SideID,flip,threadID
!REAL                            :: Uface(TP_nVar,0:Nloc,0:ZDIM(Nloc))
LOGICAL                         :: isMaster
!==================================================================================================================================
threadID = (blockidx%x-1) * blockdim%x + threadidx%x
ElemID  = (threadID-1)/6 + 1
locSide = threadID-(ElemID-1)*6
IF (ElemID.LE.nElems) THEN
  !DO locSide=1,6
    SideID   = ElemToSide(E2S_SIDE_ID  ,locSide,ElemID)
    flip     = ElemToSide(E2S_FLIP     ,locSide,ElemID)
    isMaster = ElemToSide(E2S_IS_MASTER,locSide,ElemID).EQ.1 ! master side for current elem

    !CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,ElemID),Uface,L_Minus,L_Plus,locSide)

    IF(isMaster) THEN
      CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,ElemID),Uface(:,:,:,SideID,1),L_Minus,L_Plus,locSide)
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        !Uface_master(:,p,q,SideID)=Uface(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide))
        Uface_master(:,p,q,SideID)=Uface(:,S2V2(1,p,q,0,locSide),S2V2(2,p,q,0,locSide),SideID,1)
      END DO; END DO
    ELSE ! slave
      CALL EvalElemFaceG_GPU(Nloc,UVol(:,:,:,:,ElemID),Uface(:,:,:,SideID,2),L_Minus,L_Plus,locSide)
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        !Uface_slave( :,p,q,SideID)=Uface(:,S2V2(1,p,q,flip,locSide),S2V2(2,p,q,flip,locSide))
        Uface_slave( :,p,q,SideID)=Uface(:,S2V2(1,p,q,flip,locSide),S2V2(2,p,q,flip,locSide),SideID,2)
      END DO; END DO
    END IF
  !END DO !locSideID
END IF

END SUBROUTINE ProlongToFace_Kernel_Elem_locSide_PreAlloc


END MODULE MOD_ProlongToFaceCons

!==================================================================================================================================
!> Contains routines to interpolate the primitive interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFacePrim
IMPLICIT NONE
PRIVATE
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE ProlongToFacePrim
  MODULE PROCEDURE ProlongToFace
END INTERFACE

PUBLIC::ProlongToFacePrim

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFacePrim

!==================================================================================================================================
!> Contains routines to interpolate the primitive interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFaceLifting
IMPLICIT NONE
PRIVATE
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarLifting

INTERFACE ProlongToFaceLifting
  MODULE PROCEDURE ProlongToFace
END INTERFACE

PUBLIC::ProlongToFaceLifting

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFaceLifting

!==================================================================================================================================
!> Contains routines to interpolate a scalar interior solution to the boundary
!==================================================================================================================================
MODULE MOD_ProlongToFace1
! MODULES
IMPLICIT NONE
PRIVATE
#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = 1

INTERFACE ProlongToFace1
  MODULE PROCEDURE ProlongToFace
END INTERFACE

PUBLIC::ProlongToFace1

CONTAINS
#include "prolongtoface.t90"
END MODULE MOD_ProlongToFace1

!==================================================================================================================================
