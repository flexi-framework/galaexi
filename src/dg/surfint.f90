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
USE MOD_Mesh_Vars ,ONLY: SideToElem,nSides
USE MOD_Mesh_Vars ,ONLY: S2V2,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc        !< (IN) Polynomial degree
LOGICAL,INTENT(IN) :: doMPISides  !<= .TRUE. only MPISides_YOUR+MPIMortar are filled
                                  !<=.FALSE. BCSides+(Mortar-)InnerSides+MPISides_MINE
REAL,DEVICE,INTENT(IN)    :: Flux_master(1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on master side
REAL,DEVICE,INTENT(IN)    :: Flux_slave (1:TP_nVar,0:Nloc,0:ZDIM(Nloc),nSides) !< (IN) Flux on slave side
!> (IN) Lagrange polynomials evaluated at \f$\xi=+1\f$ and \f$\xi=-1\f$ and premultiplied by mass matrix
REAL,INTENT(IN)    :: L_HatPlus(0:Nloc),L_HatMinus(0:Nloc)
REAL,DEVICE,INTENT(INOUT) :: Ut(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< (INOUT) Time derivative of the solution
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DEVICE                     :: d_L_HatMinus(0:Nloc),d_L_HatPlus(0:Nloc)
INTEGER,DEVICE                  :: d_SideToElem(5,nSides)
INTEGER,DEVICE                  :: d_S2V2(2,0:Nloc,0:Nloc,0:4,6)
REAL,DEVICE :: d_Ut_surf(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< (INOUT) Time derivative of the solution
REAL        :: Ut_surf(TP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems)   !< (INOUT) Time derivative of the solution
INTEGER     :: iElem,i,j,k
!==================================================================================================================================
d_L_HatMinus   = L_HatMinus
d_L_HatPlus    = L_HatPlus
d_SideToElem   = SideToElem
d_S2V2         = S2V2
Ut_surf = 0.
d_Ut_surf = Ut_surf
d_Ut_surf = 0.

!CALL SurfIntCons_Kernel<<<nSides/32+1,32>>>(Nloc,nSides,nElems,Flux_master,Flux_slave,d_Ut_Surf,d_L_HatMinus,d_L_HatPlus,d_SideToElem,d_S2V2)
CALL SurfIntCons_Kernel<<<nSides,1>>>(Nloc,nSides,nElems,Flux_master,Flux_slave,d_Ut_Surf,d_L_HatMinus,d_L_HatPlus,d_SideToElem,d_S2V2)
!$cuf kernel do(4) <<< *, 256 >>>
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+d_Ut_surf(:,i,j,k,iElem)
  END DO; END DO; END DO
END DO
!CALL SurfIntCons_Kernel<<<nSides/32+1,32>>>(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,d_L_HatMinus,d_L_HatPlus,d_SideToElem,d_S2V2)
END SUBROUTINE SurfIntCons_GPU

!==================================================================================================================================
!> In this routine, the surface integral will be computed
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE SurfIntCons_Kernel(Nloc,nSides,nElems,Flux_master,Flux_slave,Ut,L_HatMinus,L_HatPlus,SideToElem,S2V2)
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
INTEGER,INTENT(IN)              :: SideToElem(5,nSides)
INTEGER,INTENT(IN)              :: S2V2(2,0:Nloc,0:Nloc,0:4,6)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,nbElemID,locSideID,nblocSideID,SideID,p,q,flip
INTEGER            :: firstSideID,lastSideID
REAL               :: FluxTmp(1:TP_nVar,0:Nloc,0:ZDIM(Nloc))
!==================================================================================================================================
SideID = (blockidx%x-1) * blockdim%x + threadidx%x
IF (SideID.LE.nSides) THEN
  ElemID      = SideToElem(S2E_ELEM_ID,   SideID)
  nbElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)

  ! master sides
  IF(ElemID.GT.0)THEN
    locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
    flip      = 0
    ! orient flux to fit flip and locSide to element local system
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! note: for master sides, the mapping S2V2 should be a unit matrix
      FluxTmp(:,S2V2(1,p,q,flip,locSideID),S2V2(2,p,q,flip,locSideID)) = Flux_master(:,p,q,SideID)
    END DO; END DO ! p,q
    CALL DoSurfInt_GPU( Nloc,FluxTmp,L_HatMinus, L_HatPlus,locSideID,Ut(:,:,:,:,ElemID))
  END IF

  ! slave sides
  IF(nbElemID.GT.0)THEN
    nblocSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
    flip        = SideToElem(S2E_FLIP,SideID)
    ! orient flux to fit flip and locSide to element local system
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! p,q are in the master RHS system, they need to be transformed to the slave volume system using S2V2 mapping
      FluxTmp(:,S2V2(1,p,q,flip,nblocSideID),S2V2(2,p,q,flip,nblocSideID))=-Flux_slave(:,p,q,SideID)
    END DO; END DO ! p,q
    CALL DoSurfInt_GPU(Nloc,FluxTmp,L_HatMinus,L_HatPlus,nblocSideID,Ut(:,:,:,:,nbElemID))
  END IF
END IF
END SUBROUTINE SurfIntCons_Kernel


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
