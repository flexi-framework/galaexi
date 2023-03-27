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
!> Routines providing general math functions
!==================================================================================================================================
MODULE MOD_Mathtools
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE INVERSE
   MODULE PROCEDURE INVERSE
END INTERFACE

INTERFACE CROSS
  MODULE PROCEDURE CROSS
END INTERFACE CROSS

PUBLIC::INVERSE
PUBLIC::CROSS
PUBLIC::GlobalVectorDotProduct
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Computes matrix inverse using LAPACK
!> Input matrix should be a square matrix
!==================================================================================================================================
FUNCTION INVERSE(A) RESULT(AINV)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: A(:,:)                      !< input matrix
REAL             :: AINV(SIZE(A,1),SIZE(A,2))   !< result: inverse of A
!----------------------------------------------------------------------------------------------------------------------------------
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI
! LOCAL VARIABLES
REAL    :: sdet
REAL    :: work(SIZE(A,1))  ! work array for LAPACK
INTEGER :: ipiv(SIZE(A,1))  ! pivot indices
INTEGER :: n,info
!==================================================================================================================================
! Store A in Ainv to prevent it from being overwritten by LAPACK
Ainv = A
n = size(A,1)

IF (n.EQ.2) THEN
  sdet=1./(A(1,1) * A(2,2) - A(1,2)*A(2,1))
  AINV(1,1) = (  A(2,2) ) * sdet
  AINV(1,2) = (- A(1,2) ) * sdet
  AINV(2,1) = (- A(2,1) ) * sdet
  AINV(2,2) = (  A(1,1) ) * sdet
ELSE IF (n.EQ.3) THEN
  sdet = 1./ (  ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) * A(3,3) &
              + ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) * A(3,1) &
              + ( A(1,3) * A(2,1) - A(1,1) * A(2,3) ) * A(3,2))
  AINV(1,1) = ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) * sdet
  AINV(1,2) = ( A(1,3) * A(3,2) - A(1,2) * A(3,3) ) * sdet
  AINV(1,3) = ( A(1,2) * A(2,3) - A(1,3) * A(2,2) ) * sdet
  AINV(2,1) = ( A(2,3) * A(3,1) - A(2,1) * A(3,3) ) * sdet
  AINV(2,2) = ( A(1,1) * A(3,3) - A(1,3) * A(3,1) ) * sdet
  AINV(2,3) = ( A(1,3) * A(2,1) - A(1,1) * A(2,3) ) * sdet
  AINV(3,1) = ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) * sdet
  AINV(3,2) = ( A(1,2) * A(3,1) - A(1,1) * A(3,2) ) * sdet
  AINV(3,3) = ( A(1,1) * A(2,2) - A(1,2) * A(2,1) ) * sdet
ELSE
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  CALL DGETRF(n, n, Ainv, n, ipiv, info)

  IF(info.NE.0)THEN
    STOP 'Matrix is numerically singular!'
  END IF

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  CALL DGETRI(n, Ainv, n, ipiv, work, n, info)

  IF(info.NE.0)THEN
    STOP 'Matrix inversion failed!'
  END IF
END IF
END FUNCTION INVERSE


!!===================================================================================================================================
!!> Bruce Dawson quote:
!!> "There is no silver bullet. You have to choose wisely."
!!>    * "If you are comparing against zero, then relative epsilons and ULPs based comparisons are usually meaningless.
!!>      You’ll need to use an absolute epsilon, whose value might be some small multiple of FLT_EPSILON and the inputs
!!>      to your calculation. Maybe."
!!>    * "If you are comparing against a non-zero number then relative epsilons or ULPs based comparisons are probably what you want.
!!>      You’ll probably want some small multiple of FLT_EPSILON for your relative epsilon, or some small number of ULPs.
!!>      An absolute epsilon could be used if you knew exactly what number you were comparing against."
!!>    * "If you are comparing two arbitrary numbers that could be zero or non-zero then you need the kitchen sink.
!!>      Good luck and God speed."
!!>
!!>      NOTE: The functions below are implemented as preprocessor macros, which are by definition
!!>            inlined and are thus beneficial in terms of performance and accuracy as they do not
!!>            depend on the data type.
!!>
!!===================================================================================================================================
!PURE FUNCTION ALMOSTEQUALRELATIVE(x,y,tol)
!! MODULES
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!REAL,INTENT(IN) :: x                !< (IN)  first scalar to be compared
!REAL,INTENT(IN) :: y                !< (IN)  second scalar to be compared
!REAL,INTENT(IN) :: tol              !< (IN) relative epsilon value as input
!LOGICAL         :: ALMOSTEQUALRELATIVE
!!==================================================================================================================================
!ALMOSTEQUALRELATIVE=(ABS(x-y).LE.MAX(ABS(x),ABS(y))*tol)
!END FUNCTION ALMOSTEQUALRELATIVE
!PURE FUNCTION ALMOSTEQUALABSOLUTE(x,y,tol)
!! MODULES
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!REAL,INTENT(IN) :: x                !< (IN)  first scalar to be compared
!REAL,INTENT(IN) :: y                !< (IN)  second scalar to be compared
!REAL,INTENT(IN) :: tol              !< (IN) relative epsilon value as input
!LOGICAL         :: ALMOSTEQUALABSOLUTE
!!==================================================================================================================================
!ALMOSTEQUALRELATIVE=(ABS(x-y).LE.tol)
!
!END FUNCTION ALMOSTEQUALABSOLUTE


!==================================================================================================================================
!> computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
!==================================================================================================================================
PPURE FUNCTION CROSS(v1,v2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    !< input vector 1
REAL,INTENT(IN) :: v2(3)    !< input vector 2
REAL            :: CROSS(3) !< cross product of vectors
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CROSS=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS

!===================================================================================================================================
!> Computes the global dot Product for vectors a and b: resu=a.b. Each processor computes the local dot product, then an allreduce
!> is called.
!===================================================================================================================================
SUBROUTINE GlobalVectorDotProduct(A,B,nDOFProc,Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: A(nDOFProc)
REAL,INTENT(IN)    :: B(nDOFProc)
INTEGER,INTENT(IN) :: nDOFProc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
!===================================================================================================================================
Resu=0.
DO i=1,nDOFProc
  Resu=Resu + A(i)*B(i)
END DO

#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_FLEXI,iError)
#endif

END SUBROUTINE GlobalVectorDotProduct

END MODULE MOD_Mathtools
