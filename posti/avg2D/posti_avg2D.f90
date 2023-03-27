!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
!> Tool which takes a HDF5 file (no matter what kind) that belongs to a ijk-sorted 3D mesh and averages all datasets in the file in
!> the third dimension. Works for either pointwise or elementwise data sets.
!> The original file will be copied, so it keeps all attributes, user block contents etc.
!>
!> Usage: posti_avg2D statefile.h5
!===================================================================================================================================
PROGRAM avg2D
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_IO_HDF5
USE MOD_HDF5_Input
USE MOD_HDF5_Output
USE MOD_MPI,                     ONLY: InitMPI
#if USE_MPI
USE MOD_MPI,                     ONLY: FinalizeMPI
#endif
USE MOD_Interpolation,           ONLY: GetNodesAndWeights
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255),ALLOCATABLE       :: tmpDatasetNames(:)
INTEGER,ALLOCATABLE                  :: Elem_IJK(:,:)
INTEGER                              :: nElems_IJK(3)
INTEGER                              :: nElems,iElem,iDataset,l
REAL,ALLOCATABLE                     :: RealArray(:,:,:,:,:),RealAvg(:,:,:,:,:)
REAL,ALLOCATABLE                     :: RealElemArray(:,:),RealElemAvg(:,:,:)
INTEGER                              :: nVar,N
INTEGER                              :: p,q,i,j
CHARACTER(LEN=255)                   :: MeshFile,NodeType,NewFileName
REAL,ALLOCATABLE                     :: xGP(:),wGP(:)
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()
! Check for correct number of input files
IF (nArgs.NE.1) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti_avg2D statefile.h5')
END IF

! Read the mesh file and node type from the statefile
CALL OpenDataFile(TRIM(Args(1)),create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
CALL ReadAttribute(File_ID,'NodeType',1,StrScalar=NodeType)
CALL CloseDataFile()


! Read in the ijk sorting of the mesh
CALL OpenDataFile(TRIM(MeshFile),create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
CALL GetDataSize(File_ID,'Elem_IJK',nDims,HSize)
nElems = INT(HSIZE(nDims))
ALLOCATE(Elem_IJK(3,nElems))
CALL ReadArray('Elem_IJK',2,(/3,nElems/),0,2,IntArray=Elem_IJK)
CALL ReadArray('nElems_IJK',1,(/3/),0,1,IntArray=nElems_IJK)
CALL CloseDataFile()


! Get the names of the data sets in the file, and information about array size etc.
CALL OpenDataFile(TRIM(Args(1)),create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
CALL GetDatasetNamesInGroup("/",tmpDatasetNames)
CALL CloseDataFile()

! Copy the current file, so we keep all the attributes etc.
NewFileName = Args(1)(:LEN(TRIM(Args(1)))-3)//'_avg2D.h5'
IF (MPIRoot) CALL EXECUTE_COMMAND_LINE("cp -f "//TRIM(Args(1))//" "//TRIM(NewFileName))

! Loop over all the datasets
DO iDataset = 1, SIZE(tmpDatasetNames)
  ! Read in the elementwise or pointwise arrays
  WRITE(*,*) ''
  WRITE(*,*) 'Read dataset ',TRIM(tmpDatasetNames(iDataset))
  CALL OpenDataFile(TRIM(Args(1)),create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
  CALL GetDataSize(File_ID,TRIM(tmpDatasetNames(iDataset)),nDims,HSize)
  IF (nDims.EQ.2) THEN
    ! Elementwise data set
    ALLOCATE(RealElemArray(INT(HSize(1)),INT(HSize(2))))
    CALL ReadArray(TRIM(tmpDatasetNames(iDataset)),2,INT(HSize),0,2,RealArray=RealElemArray)
    nVar = INT(HSize(1))
  ELSE IF (nDims.EQ.5) THEN
    ! Pointwise data set
    ALLOCATE(RealArray(INT(HSize(1)),0:INT(HSize(2))-1,0:INT(HSize(3))-1,0:INT(HSize(4))-1,INT(HSize(5))))
    CALL ReadArray(TRIM(tmpDatasetNames(iDataset)),5,INT(HSize),0,5,RealArray=RealArray)
    nVar = INT(HSize(1))
    N = INT(HSize(2)-1)
    ! Prepare integration weights for averaging
    ALLOCATE(xGP(0:N))
    ALLOCATE(wGP(0:N))
    CALL GetNodesAndWeights(N,TRIM(NodeType),xGP,wGP)
  ELSE
    CALL CloseDataFile()
    WRITE(*,*) 'Skip dataset ',TRIM(tmpDatasetNames(iDataset)), ' (wrong dimension)'
    CYCLE
  END IF
  CALL CloseDataFile()


  ! Compute the averages
  IF (nDims.EQ.2) THEN
    ! Elementwise data set
    ALLOCATE(RealElemAvg(nVar,nElems_IJK(1),nElems_IJK(2)))
    RealElemAvg = 0.
    DO iElem=1,nElems
      ! Indixes in the IJK sorting
      i = Elem_IJK(1,iElem)
      j = Elem_IJK(2,iElem)
      ! Build sum
      RealElemAvg(:,i,j) = RealElemAvg(:,i,j) + RealElemArray(:,iElem)
    END DO ! iElem
    ! Finish sum
    RealElemAvg = RealElemAvg / nElems_IJK(3)
    ! Write back to the array itself
    DO iElem=1,nElems
      ! Indixes in the IJK sorting
      i = Elem_IJK(1,iElem)
      j = Elem_IJK(2,iElem)
      RealElemArray(:,iElem) = RealElemAvg(:,i,j)
    END DO ! iElem
  ELSE IF (nDims.EQ.5) THEN
    ! Pointwise data set
    ALLOCATE(RealAvg(nVar,0:N,0:N,nElems_IJK(1),nElems_IJK(2)))
    RealAvg = 0.
    DO iElem=1,nElems
      ! Indixes in the IJK sorting
      i = Elem_IJK(1,iElem)
      j = Elem_IJK(2,iElem)
      ! Build sum
      DO q=0,N; DO p=0,N
        DO l = 0,N
          RealAvg(:,p,q,i,j) = RealAvg(:,p,q,i,j) + RealArray(:,p,q,l,iElem)*wGP(l)/2.
        END DO ! l = 0,N
      END DO; END DO ! p,q=0,PP_N
    END DO ! iElem
    ! Finish sum
    RealAvg = RealAvg / (nElems_IJK(3))
    ! Write back to the array itself
    DO iElem=1,nElems
      ! Indixes in the IJK sorting
      i = Elem_IJK(1,iElem)
      j = Elem_IJK(2,iElem)
      ! Build sum
      DO q=0,N; DO p=0,N
        DO l = 0,N
          RealArray(:,p,q,l,iElem) = RealAvg(:,p,q,i,j)
        END DO ! l = 0,N
      END DO; END DO ! p,q=0,PP_N
    END DO ! iElem
  END IF


  ! Open new file and write the array
  WRITE(*,*) 'Write dataset ',TRIM(tmpDatasetNames(iDataset))
  CALL OpenDataFile(TRIM(NewFileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
  IF (nDims.EQ.2) THEN
    CALL WriteArray(TRIM(tmpDatasetNames(iDataset)),2,&
                    INT(HSize),&
                    INT(HSize),&
                    (/0,0/),.FALSE.,RealArray=RealElemArray)
  ELSE IF (nDims.EQ.5) THEN
    CALL WriteArray(TRIM(tmpDatasetNames(iDataset)),5,&
                    INT(HSize),&
                    INT(HSize),&
                    (/0,0,0,0,0/),.FALSE.,RealArray=RealArray)
  END IF
  CALL CloseDataFile()

  SDEALLOCATE(RealArray)
  SDEALLOCATE(RealAvg)
  SDEALLOCATE(RealElemArray)
  SDEALLOCATE(RealElemAvg)
  SDEALLOCATE(wGP)
  SDEALLOCATE(xGP)

END DO ! iDataset = 1, SIZE(tmpDatasetNames)
SDEALLOCATE(tmpDatasetNames)
SDEALLOCATE(Elem_IJK)


#if USE_MPI
CALL FinalizeMPI()
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' AVG2D TOOL FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM avg2D
