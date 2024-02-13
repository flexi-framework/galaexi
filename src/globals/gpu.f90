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

INTEGER(KIND=CUDA_STREAM_KIND) :: DefaultStream
INTEGER(KIND=CUDA_STREAM_KIND) :: stream1
INTEGER(KIND=CUDA_STREAM_KIND) :: stream2
INTEGER(KIND=CUDA_STREAM_KIND) :: stream3
INTEGER(KIND=CUDA_STREAM_KIND) :: stream4
INTEGER(KIND=CUDA_STREAM_KIND) :: stream5
INTEGER(KIND=CUDA_STREAM_KIND) :: stream6

INTERFACE DefineParametersGPU
  MODULE PROCEDURE DefineParametersGPU
END INTERFACE

INTERFACE InitGPU
  MODULE PROCEDURE InitGPU
END INTERFACE

INTERFACE FinalizeGPU
  MODULE PROCEDURE FinalizeGPU
END INTERFACE

PUBLIC::DefineParametersGPU,InitGPU,FinalizeGPU
PUBLIC::stream1,stream2,stream3,DefaultStream
PUBLIC::stream4,stream5,stream6
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> DefineLifting
!==================================================================================================================================
SUBROUTINE DefineParametersGPU()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("GPU")
CALL prms%CreateLogicalOption('useStreams', "Set true to use concurrent streaming on device.", '.TRUE.')
END SUBROUTINE DefineParametersGPU

!==================================================================================================================================
!> Initialize indicators
!==================================================================================================================================
SUBROUTINE InitGPU()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_ReadInTools,ONLY: GETLOGICAL
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL               :: useStreams
INTEGER               :: prio_min,prio_max,prio_mid
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT GPU...'

useStreams = GETLOGICAL('useStreams')

DefaultStream=cudaforGetDefaultStream()
IF (useStreams) THEN
  ! Get priority range of current GPU
  iError = cudaDeviceGetStreamPriorityRange(prio_min,prio_max)
  prio_mid = MIN(prio_max+1,prio_min) ! Lower values mean higher priority

  ! Create individual streams with corresponding priority:
  ! Highest: MPI Sides
  ! Lowest:  Inner Volume work
  ! Mid:     All the rest of local work
  iError = cudaStreamCreateWithPriority(stream1,cudaStreamDefault,prio_min) ! Volume
  iError = cudaStreamCreateWithPriority(stream2,cudaStreamDefault,prio_mid) ! Inner Sides
  iError = cudaStreamCreateWithPriority(stream3,cudaStreamDefault,prio_max) ! MPI Sides
  iError = cudaStreamCreateWithPriority(stream4,cudaStreamDefault,prio_mid) ! gradUx
  iError = cudaStreamCreateWithPriority(stream5,cudaStreamDefault,prio_mid) ! gradUy
  iError = cudaStreamCreateWithPriority(stream6,cudaStreamDefault,prio_mid) ! gradUz
ELSE
  stream1 = DefaultStream
  stream2 = DefaultStream
  stream3 = DefaultStream
  stream4 = DefaultStream
  stream5 = DefaultStream
  stream6 = DefaultStream
END IF

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
istat = cudaStreamDestroy(stream3)
istat = cudaStreamDestroy(stream4)
istat = cudaStreamDestroy(stream5)
istat = cudaStreamDestroy(stream6)

END SUBROUTINE FinalizeGPU

END MODULE MOD_GPU
