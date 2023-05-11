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

MODULE MOD_GPU
! MODULES
USE CUDAFOR
IMPLICIT NONE

PRIVATE

INTEGER(KIND=cuda_stream_kind) :: stream1
INTEGER(KIND=cuda_stream_kind) :: stream2

INTERFACE InitGPU
  MODULE PROCEDURE InitGPU
END INTERFACE

INTERFACE FinalizeGPU
  MODULE PROCEDURE FinalizeGPU
END INTERFACE

PUBLIC::InitGPU,FinalizeGPU
PUBLIC::stream1,stream2
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Initialize indicators
!==================================================================================================================================
SUBROUTINE InitGPU()
! MODULES
USE MOD_Preproc
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               ::   istat
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT GPU...'

SWRITE(UNIT_stdOut,'(A)') ' DUMMY ' 


istat = cudaStreamCreate(stream1)
istat = cudaStreamCreate(stream2)

SWRITE(UNIT_stdOut,'(A)')' INIT GPU DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitGPU

!==================================================================================================================================
!> Initialize indicators
!==================================================================================================================================
SUBROUTINE FinalizeGPU()
! MODULES
USE MOD_Preproc
USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               ::   istat
!==================================================================================================================================

istat = cudaStreamDestroy(stream1)
istat = cudaStreamDestroy(stream2)

END SUBROUTINE FinalizeGPU

END MODULE MOD_GPU
