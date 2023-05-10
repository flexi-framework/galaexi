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
!> Contains the different Surface integral formulations
!> Computes the Surface integral for all faces using U and updates Ut
!> Computes only inner surface integrals!
!> Surface integrals are separated for each direction
!==================================================================================================================================
MODULE MOD_SurfInt
IMPLICIT NONE
PRIVATE

#define WITHnVars 1

INTERFACE SurfInt
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfInt
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfInt,DoSurfInt

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfInt

!==================================================================================================================================
!> Contains the surface integral for conservative quantities
!==================================================================================================================================
MODULE MOD_SurfIntCons
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE SurfIntCons
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE SurfIntCons_GPU
  MODULE PROCEDURE SurfIntCons_GPU
END INTERFACE

INTERFACE DoSurfIntCons
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfIntCons,DoSurfIntCons,SurfIntCons_GPU

CONTAINS
#include "surfint.t90"

!==================================================================================================================================
!> In this routine, the surface integral will be computed
!==================================================================================================================================
SUBROUTINE SurfIntCons_GPU(Nloc,Flux_master,Flux_slave,Ut,doMPISides,L_HatMinus,L_HatPlus)
! MODULES
USE MOD_DG_Vars   ,ONLY: nDOFElem
USE MOD_Mesh_Vars ,ONLY: SideToElem,nSides,ElemToSide
USE MOD_Mesh_Vars ,ONLY: S2V2,nElems
!@cuf USE MOD_Mesh_Vars ,ONLY: d_SideToElem,d_ElemToSide,d_S2V2,d_S2V2_inv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc        !< (IN) Polynomial degree
LOGICAL,INTENT(IN) :: doMPISides  !<= .TRUE. only MPISides_YOUR+MPIMortar are filled
                                  !<=.FALSE. BCSides+(Mortar-)InnerSides+MPISides_MINE
REAL,DEVICE,INTENT(IN)    :: Flux_master(1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on master side
REAL,DEVICE,INTENT(IN)    :: Flux_slave (1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on slave side
!> (IN) Lagrange polynomials evaluated at \f$\xi=+1\f$ and \f$\xi=-1\f$ and premultiplied by mass matrix
REAL,DEVICE,INTENT(IN)    :: L_HatPlus(0:Nloc),L_HatMinus(0:Nloc)
REAL,DEVICE,INTENT(INOUT) :: Ut(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< (INOUT) Time derivative of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER     :: nThreads=128
!==================================================================================================================================

!CALL SurfIntCons_Kernel<<<nElems/nThreads+1,nThreads>>>(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,d_L_HatMinus,d_L_HatPlus,d_ElemToSide,d_S2V2)
!CALL SurfIntCons_Kernel_Point<<<nElems*nDOFElem/nThreads+1,nThreads>>>(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,L_HatMinus,L_HatPlus,d_ElemToSide,d_S2V2_inv)
CALL SurfIntCons_Kernel_Point_Contract<<<nElems*nDOFElem/nThreads+1,nThreads>>>(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,L_HatMinus,L_HatPlus,d_ElemToSide,d_S2V2_inv)
END SUBROUTINE SurfIntCons_GPU

!==================================================================================================================================
!> In this routine, the surface integral will be computed
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE SurfIntCons_Kernel(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,L_HatMinus,L_HatPlus,ElemToSide,S2V2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN) :: Nloc  !< (IN) Polynomial degree
INTEGER,VALUE,INTENT(IN) :: nSides
INTEGER,VALUE,INTENT(IN) :: nElems
REAL,INTENT(IN)    :: Flux_master(1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on master side
REAL,INTENT(IN)    :: Flux_slave (1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on slave side
REAL,INTENT(IN)    :: L_HatPlus(0:Nloc),L_HatMinus(0:Nloc)
REAL,INTENT(INOUT) :: Ut(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< (INOUT) Time derivative of the solution
INTEGER,INTENT(IN)              :: ElemToSide(3,6,nElems)
!INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,nbElemID,locSideID,nblocSideID,SideID,p,q,flip
INTEGER            :: firstSideID,lastSideID
REAL               :: FluxTmp(1:TP_nVar,0:Nloc,0:ZDIM(Nloc))
LOGICAL            :: isMaster
!==================================================================================================================================
ElemID = (blockidx%x-1) * blockdim%x + threadidx%x
IF (ElemID.LE.nElems) THEN
  DO locSideID=1,6
    SideID   = ElemToSide(E2S_SIDE_ID  ,locSideID,ElemID)
    flip     = ElemToSide(E2S_FLIP     ,locSideID,ElemID)
    isMaster = ElemToSide(E2S_IS_MASTER,locSideID,ElemID).EQ.1 ! master side for current elem

    IF(isMaster) THEN
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        ! note: for master sides, the mapping S2V2 should be a unit matrix
        FluxTmp(:,S2V2(1,p,q,flip,locSideID),S2V2(2,p,q,flip,locSideID)) = Flux_master(:,p,q,SideID)
      END DO; END DO ! p,q
    ELSE ! slave
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        ! note: for master sides, the mapping S2V2 should be a unit matrix
        FluxTmp(:,S2V2(1,p,q,flip,locSideID),S2V2(2,p,q,flip,locSideID)) =-Flux_slave(:,p,q,SideID)
      END DO; END DO ! p,q
    END IF
    CALL DoSurfInt_GPU( Nloc,FluxTmp,L_HatMinus, L_HatPlus,locSideID,Ut(:,:,:,:,ElemID))
  END DO !locSideID
END IF
END SUBROUTINE SurfIntCons_Kernel

!==================================================================================================================================
!> In this routine, the surface integral will be computed
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE SurfIntCons_Kernel_Point(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,L_HatMinus,L_HatPlus,ElemToSide,S2V2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN) :: Nloc  !< (IN) Polynomial degree
INTEGER,VALUE,INTENT(IN) :: nSides
INTEGER,VALUE,INTENT(IN) :: nElems
REAL,INTENT(IN)    :: Flux_master(1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on master side
REAL,INTENT(IN)    :: Flux_slave (1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on slave side
REAL,INTENT(IN)    :: L_HatPlus(0:Nloc),L_HatMinus(0:Nloc)
REAL,INTENT(INOUT) :: Ut(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< (INOUT) Time derivative of the solution
INTEGER,INTENT(IN)              :: ElemToSide(3,6,nElems)
!INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,locSideID,SideID,p,q,flip,rest
INTEGER            :: i,j,k,threadID
INTEGER            :: firstSideID,lastSideID
!REAL               :: FluxTmp(1:TP_nVar,0:Nloc,0:ZDIM(Nloc))
REAL               :: FluxTmp(1:TP_nVar)
LOGICAL            :: isMaster
!==================================================================================================================================
! Get thread indices
threadID = (blockidx%x-1) * blockdim%x + threadidx%x
! Get ElemID of current thread
ElemID   =        (threadID-1)/(Nloc+1)**3+1 ! Elems are 1-indexed
rest     = threadID-(ElemID-1)*(Nloc+1)**3
! Get ijk indices of current thread
k        = (rest-1)/(Nloc+1)**2
rest     =  rest- k*(Nloc+1)**2
j        = (rest-1)/(Nloc+1)!**1
rest     =  rest- j*(Nloc+1)!**1
i        = (rest-1)!/(Nloc+1)**0
rest     =  rest- j!*(Nloc+1)**0

IF (ElemID.LE.nElems) THEN
  DO locSideID=1,6
    SideID   = ElemToSide(E2S_SIDE_ID  ,locSideID,ElemID)
    flip     = ElemToSide(E2S_FLIP     ,locSideID,ElemID)
    isMaster = ElemToSide(E2S_IS_MASTER,locSideID,ElemID).EQ.1 ! master side for current elem

    ! Get indices of pointwise flux on locside
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      p=S2V2(1,j,k,flip,locSideID)
      q=S2V2(2,j,k,flip,locSideID)
    CASE(ETA_MINUS)
      p=S2V2(1,i,k,flip,locSideID)
      q=S2V2(2,i,k,flip,locSideID)
    CASE(ZETA_MINUS)
      p=S2V2(1,i,j,flip,locSideID)
      q=S2V2(2,i,j,flip,locSideID)
    CASE(XI_PLUS)
      p=S2V2(1,j,k,flip,locSideID)
      q=S2V2(2,j,k,flip,locSideID)
    CASE(ETA_PLUS)
      p=S2V2(1,i,k,flip,locSideID)
      q=S2V2(2,i,k,flip,locSideID)
    CASE(ZETA_PLUS)
      p=S2V2(1,i,j,flip,locSideID)
      q=S2V2(2,i,j,flip,locSideID)
    END SELECT !locSideID

    ! Extract flux either from master or slave array
    IF(isMaster) THEN
      FluxTmp(:) = Flux_master(:,p,q,SideID)
    ELSE ! slave
      FluxTmp(:) =-Flux_slave(:,p,q,SideID)
    END IF

    ! Add contribution from this locSide to my own DOF
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +FluxTmp(:)*L_hatMinus(i)
    CASE(ETA_MINUS)
      Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +FluxTmp(:)*L_hatMinus(j)
    CASE(ZETA_MINUS)
      Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +FluxTmp(:)*L_hatMinus(k)
    CASE(XI_PLUS)
      Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +FluxTmp(:)*L_hatPlus(i)
    CASE(ETA_PLUS)
      Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +FluxTmp(:)*L_hatPlus(j)
    CASE(ZETA_PLUS)
      Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +FluxTmp(:)*L_hatPlus(k)
    END SELECT !locSideID
  END DO
ENDIF
END SUBROUTINE SurfIntCons_Kernel_Point

!==================================================================================================================================
!> In this routine, the surface integral will be computed
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE SurfIntCons_Kernel_Point_Contract(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,L_HatMinus,L_HatPlus,ElemToSide,S2V2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN) :: Nloc  !< (IN) Polynomial degree
INTEGER,VALUE,INTENT(IN) :: nSides
INTEGER,VALUE,INTENT(IN) :: nElems
REAL,INTENT(IN)    :: Flux_master(1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on master side
REAL,INTENT(IN)    :: Flux_slave (1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on slave side
REAL,INTENT(IN)    :: L_HatPlus(0:Nloc),L_HatMinus(0:Nloc)
REAL,INTENT(INOUT) :: Ut(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< (INOUT) Time derivative of the solution
INTEGER,INTENT(IN)              :: ElemToSide(3,6,nElems)
!INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,locSideID,SideID,p,q,flip,rest
INTEGER            :: i,j,k,threadID
INTEGER            :: firstSideID,lastSideID
LOGICAL            :: isMaster
!==================================================================================================================================
! Get thread indices
threadID = (blockidx%x-1) * blockdim%x + threadidx%x
! Get ElemID of current thread
ElemID   =        (threadID-1)/(Nloc+1)**3+1 ! Elems are 1-indexed
rest     = threadID-(ElemID-1)*(Nloc+1)**3
! Get ijk indices of current thread
k        = (rest-1)/(Nloc+1)**2
rest     =  rest- k*(Nloc+1)**2
j        = (rest-1)/(Nloc+1)!**1
rest     =  rest- j*(Nloc+1)!**1
i        = (rest-1)!/(Nloc+1)**0
rest     =  rest- j!*(Nloc+1)**0

IF (ElemID.LE.nElems) THEN
  DO locSideID=1,6
    SideID   = ElemToSide(E2S_SIDE_ID  ,locSideID,ElemID)
    flip     = ElemToSide(E2S_FLIP     ,locSideID,ElemID)
    isMaster = ElemToSide(E2S_IS_MASTER,locSideID,ElemID).EQ.1 ! master side for current elem

    ! Get indices of pointwise flux on locside
    SELECT CASE(locSideID)
    CASE(XI_MINUS)
      p=S2V2(1,j,k,flip,locSideID)
      q=S2V2(2,j,k,flip,locSideID)
      IF(isMaster) THEN
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +Flux_master(:,p,q,SideID)*L_hatMinus(i)
      ELSE ! slave
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   -Flux_slave( :,p,q,SideID)*L_hatMinus(i)
      END IF
    CASE(ETA_MINUS)
      p=S2V2(1,i,k,flip,locSideID)
      q=S2V2(2,i,k,flip,locSideID)
      IF(isMaster) THEN
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +Flux_master(:,p,q,SideID)*L_hatMinus(j)
      ELSE ! slave
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   -Flux_slave( :,p,q,SideID)*L_hatMinus(j)
      END IF
    CASE(ZETA_MINUS)
      p=S2V2(1,i,j,flip,locSideID)
      q=S2V2(2,i,j,flip,locSideID)
      IF(isMaster) THEN
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +Flux_master(:,p,q,SideID)*L_hatMinus(k)
      ELSE ! slave
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   -Flux_slave( :,p,q,SideID)*L_hatMinus(k)
      END IF
    CASE(XI_PLUS)
      p=S2V2(1,j,k,flip,locSideID)
      q=S2V2(2,j,k,flip,locSideID)
      IF(isMaster) THEN
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +Flux_master(:,p,q,SideID)*L_hatPlus(i)
      ELSE ! slave
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   -Flux_slave( :,p,q,SideID)*L_hatPlus(i)
      END IF
    CASE(ETA_PLUS)
      p=S2V2(1,i,k,flip,locSideID)
      q=S2V2(2,i,k,flip,locSideID)
      IF(isMaster) THEN
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +Flux_master(:,p,q,SideID)*L_hatPlus(j)
      ELSE ! slave
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   -Flux_slave( :,p,q,SideID)*L_hatPlus(j)
      END IF
    CASE(ZETA_PLUS)
      p=S2V2(1,i,j,flip,locSideID)
      q=S2V2(2,i,j,flip,locSideID)
      IF(isMaster) THEN
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   +Flux_master(:,p,q,SideID)*L_hatPlus(k)
      ELSE ! slave
        Ut(:,i,j,k,ElemID) =Ut(:,i,j,k,ElemID)   -Flux_slave( :,p,q,SideID)*L_hatPlus(k)
      END IF
    END SELECT !locSideID
  END DO
ENDIF
END SUBROUTINE SurfIntCons_Kernel_Point_Contract

END MODULE MOD_SurfIntCons

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_SurfIntLifting
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVarLifting

INTERFACE SurfIntLifting
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfIntLifting
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfIntLifting,DoSurfIntLifting

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntLifting

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_SurfIntPrim
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE SurfIntPrim
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE DoSurfIntPrim
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::SurfIntPrim,DoSurfIntPrim

CONTAINS
#include "surfint.t90"
END MODULE MOD_SurfIntPrim
