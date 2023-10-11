!=================================================================================================================================
! Copyright (c) 2010-2023  Prof. Claus-Dieter Munz
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
!> Multiples volume data within an DG element linewise with a metric
!==================================================================================================================================
MODULE MOD_ApplyDMatrix
IMPLICIT NONE
PRIVATE

#define WITHnVars 1

INTERFACE ApplyDMatrix
  MODULE PROCEDURE ApplyDMatrix
  MODULE PROCEDURE ApplyDMatrix_CUDA
END INTERFACE

PUBLIC::ApplyDMatrix

CONTAINS
#include "applydmatrix.t90"
END MODULE MOD_ApplyDMatrix

!==================================================================================================================================
!> Contains the surface integral for conservative quantities
!==================================================================================================================================
MODULE MOD_ApplyDMatrixCons
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE ApplyDMatrixCons
  MODULE PROCEDURE ApplyDMatrix
  MODULE PROCEDURE ApplyDMatrix_CUDA
END INTERFACE

PUBLIC::ApplyDMatrixCons

CONTAINS
#include "applydmatrix.t90"
END MODULE MOD_ApplyDMatrixCons

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_ApplyDMatrixLifting
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVarLifting

INTERFACE ApplyDMatrixLifting
  MODULE PROCEDURE ApplyDMatrix
  MODULE PROCEDURE ApplyDMatrix_CUDA
END INTERFACE

PUBLIC::ApplyDMatrixLifting

CONTAINS
#include "applydmatrix.t90"
END MODULE MOD_ApplyDMatrixLifting

!==================================================================================================================================
!> Contains the surface integral for primitive quantities
!==================================================================================================================================
MODULE MOD_ApplyDMatrixPrim
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE ApplyDMatrixPrim
  MODULE PROCEDURE ApplyDMatrix
  MODULE PROCEDURE ApplyDMatrix_CUDA
END INTERFACE

PUBLIC::ApplyDMatrixPrim

CONTAINS
#include "applydmatrix.t90"
END MODULE MOD_ApplyDMatrixPrim
