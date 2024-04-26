!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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

!===================================================================================================================================
! Contains tools for array related operations.
!===================================================================================================================================
MODULE MOD_Array_Operations
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ChangeSizeArray
  MODULE PROCEDURE ChangeSizeArrayLOG1
  MODULE PROCEDURE ChangeSizeArrayLOG2
  MODULE PROCEDURE ChangeSizeArrayINT1
  MODULE PROCEDURE ChangeSizeArrayINT1_KIND8
  MODULE PROCEDURE ChangeSizeArrayINT2
  MODULE PROCEDURE ChangeSizeArrayREAL1
  MODULE PROCEDURE ChangeSizeArrayREAL2
  MODULE PROCEDURE ChangeSizeArrayREAL3
END INTERFACE

PUBLIC :: ChangeSizeArray
!===================================================================================================================================

CONTAINS

SUBROUTINE ChangeSizeArrayLOG1(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Logical 1D
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT) :: Vec(:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
LOGICAL,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL,ALLOCATABLE               :: TempVec(:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF (NewSize.GT.OldSize) THEN
  TempVec(1:OldSize) = Vec
  IF (PRESENT(Default)) TempVec(OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayLOG1


SUBROUTINE ChangeSizeArrayLOG2(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Logical 2D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,ALLOCATABLE,INTENT(INOUT) :: Vec(:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
LOGICAL,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL,ALLOCATABLE               :: TempVec(:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF (NewSize.GT.OldSize) THEN
  TempVec(:,1:OldSize) = Vec
  IF (PRESENT(Default)) TempVec(:,OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayLOG2


SUBROUTINE ChangeSizeArrayINT1(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Integer 1D
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,ALLOCATABLE,INTENT(INOUT) :: Vec(:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
INTEGER,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE               :: TempVec(:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF (NewSize.GT.OldSize) THEN
  TempVec(1:OldSize) = Vec
  IF (PRESENT(Default)) TempVec(OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayINT1


SUBROUTINE ChangeSizeArrayINT1_KIND8(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Integer 1D
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),ALLOCATABLE,INTENT(INOUT) :: Vec(:)
INTEGER,INTENT(IN)                        :: OldSize, NewSize
INTEGER,INTENT(IN),OPTIONAL               :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8),ALLOCATABLE               :: TempVec(:)
INTEGER                                   :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF (NewSize.GT.OldSize) THEN
  TempVec(1:OldSize) = Vec
  IF (PRESENT(Default)) TempVec(OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayINT1_KIND8


SUBROUTINE ChangeSizeArrayINT2(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Integer 2D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,ALLOCATABLE,INTENT(INOUT) :: Vec(:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
INTEGER,INTENT(IN),OPTIONAL       :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE               :: TempVec(:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF (NewSize.GT.OldSize) THEN
  TempVec(:,1:OldSize) = Vec
  IF(PRESENT(Default)) TempVec(:,OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayINT2


SUBROUTINE ChangeSizeArrayREAL1(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Real 1D
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,ALLOCATABLE,INTENT(INOUT)    :: Vec(:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
REAL,INTENT(IN),OPTIONAL          :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                  :: TempVec(:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF(NewSize.GT.OldSize) THEN
  TempVec(1:OldSize) = Vec
  IF (PRESENT(Default)) TempVec(OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayREAL1


SUBROUTINE ChangeSizeArrayREAL2(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Real 2D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,ALLOCATABLE,INTENT(INOUT)    :: Vec(:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
REAL,INTENT(IN),OPTIONAL          :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                  :: TempVec(:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF (NewSize.GT.OldSize) THEN
  TempVec(:,1:OldSize) = Vec
  IF (PRESENT(Default)) TempVec(:,OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayREAL2


SUBROUTINE ChangeSizeArrayREAL3(Vec,OldSize,NewSize,Default)
!===================================================================================================================================
! Change size of a vector, Real 3D, last dimension
! If NewSize>OldSize and PRESENT(Default): New entries are initialized with default
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,ALLOCATABLE,INTENT(INOUT)    :: Vec(:,:,:)
INTEGER,INTENT(IN)                :: OldSize, NewSize
REAL,INTENT(IN),OPTIONAL          :: Default
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE                  :: TempVec(:,:,:)
INTEGER                           :: ALLOCSTAT
!===================================================================================================================================
! Allocate new memory space
ALLOCATE(TempVec(SIZE(Vec,1),SIZE(Vec,2),NewSize),STAT=ALLOCSTAT)
IF (ALLOCSTAT.NE.0) CALL Abort(__STAMP__,'Cannot allocate new Array in ChangeSizeArray')

! Write old data to new memory space
IF (NewSize.GT.OldSize) THEN
  TempVec(:,:,1:OldSize) = Vec
  IF (PRESENT(Default)) TempVec(:,:,OldSize+1:NewSize) = Default
ELSE
  TempVec = Vec(:,:,1:NewSize)
END IF

! Switch array pointer
CALL MOVE_ALLOC(TempVec,Vec)

END SUBROUTINE ChangeSizeArrayREAL3

END MODULE MOD_Array_Operations
