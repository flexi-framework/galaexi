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
!> Module that provides functions for computing the solutions time history at a defined set of points ("recordpoints")
!==================================================================================================================================
MODULE MOD_RecordPoints
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersRecordPoints
  MODULE PROCEDURE DefineParametersRecordPoints
END INTERFACE

INTERFACE InitRecordPoints
  MODULE PROCEDURE InitRecordPoints
END INTERFACE

INTERFACE RecordPoints
  MODULE PROCEDURE RecordPoints
END INTERFACE

INTERFACE WriteRP
  MODULE PROCEDURE WriteRP
END INTERFACE

INTERFACE FinalizeRecordPoints
  MODULE PROCEDURE FinalizeRecordPoints
END INTERFACE

PUBLIC :: DefineParametersRecordPoints
PUBLIC :: InitRecordPoints
PUBLIC :: RecordPoints
PUBLIC :: WriteRP
PUBLIC :: FinalizeRecordPoints
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersRecordPoints()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("RecordPoints")
CALL prms%CreateLogicalOption('RP_inUse',          "Set true to compute solution history at points defined in recordpoints file.",&
                                                   '.FALSE.')
CALL prms%CreateStringOption( 'RP_DefFile',        "File containing element-local parametric recordpoint coordinates and structure.")
CALL prms%CreateIntOption(    'RP_MaxMemory',      "Maximum memory in MiB to be used for storing recordpoint state history. "//&
                                                   "If memory is exceeded before regular IO level states are written to file.",&
                                                   '100')
CALL prms%CreateIntOption(    'RP_SamplingOffset', "Multiple of timestep at which recordpoints are evaluated.",&
                                                   '1')
END SUBROUTINE DefineParametersRecordPoints


!==================================================================================================================================
!> Read RP parameters from ini file and RP definitions from HDF5
!==================================================================================================================================
SUBROUTINE InitRecordPoints(nVar)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools         ,ONLY: GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Interpolation_Vars  ,ONLY: InterpolationInitIsDone
USE MOD_RecordPoints_Vars
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL :: nVar
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: RP_maxMemory
INTEGER                     :: maxRP
INTEGER                     :: nVar_loc
!==================================================================================================================================
! check if recordpoints are activated
RP_inUse = GETLOGICAL('RP_inUse')
IF(.NOT.RP_inUse) RETURN

IF((.NOT.InterpolationInitIsDone) .OR. RecordPointsInitIsDone) &
   CALL Abort(__STAMP__,"InitRecordPoints not ready to be called or already called.")

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT RECORDPOINTS...'

RPDefFile = GETSTR('RP_DefFile')                ! Filename with RP coords
CALL ReadRPList(RPDefFile)                      ! RP_inUse is set to FALSE by ReadRPList if no RP is on proc.
maxRP = nGlobalRP
#if USE_MPI
CALL InitRPCommunicator()
#endif /*USE_MPI*/

IF (PRESENT(nVar)) THEN; nVar_loc = nVar
ELSE                   ; nVar_loc = PP_nVar
END IF

RP_maxMemory      = GETINT('RP_MaxMemory')      ! Max buffer (100MB)
RP_SamplingOffset = GETINT('RP_SamplingOffset') ! Sampling offset (iteration)
IF (RP_onProc) THEN
  maxRP = nGlobalRP
#if USE_MPI
  CALL MPI_ALLREDUCE(nRP,maxRP,1,MPI_INTEGER,MPI_MAX,RP_COMM,iError)
#endif /*USE_MPI*/
  RP_MaxBufferSize = RP_MaxMemory*131072/(maxRP*(nVar_loc+1)) ! = size in bytes/(real*maxRP*nVar)
  iSample          = 0
  RP_BufferSize    = 1
  ALLOCATE(RP_Data(   1:nVar_loc+1,nRP,RP_BufferSize))
!@cuf ALLOCATE(d_U_RP(nVar_loc,nRP))
  ALLOCATE(lastSample(1:nVar_loc+1,nRP))
  lastSample       = 0.
END IF

RecordPointsInitIsDone = .TRUE.

SWRITE(UNIT_stdOut,'(A)')' INIT RECORDPOINTS DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitRecordPoints


#if USE_MPI
!==================================================================================================================================
!> Read RP parameters from ini file and RP definitions from HDF5
!==================================================================================================================================
SUBROUTINE InitRPCommunicator()
! MODULES
USE MOD_Globals
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc,myRPrank,RP_COMM,nRP_Procs
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: color
!==================================================================================================================================

!--- Split communicator from MPI_COMM_FLEXI
color = MERGE(2,MPI_UNDEFINED,RP_onProc)

! create new RP communicator for RP output. Pass MPI_INFO_NULL as rank to follow the original ordering
CALL MPI_COMM_SPLIT(MPI_COMM_FLEXI, color, MPI_INFO_NULL, RP_COMM, iError)

! return if proc not on RP_COMM
IF (.NOT. RP_onProc) RETURN

! Find my rank on the RP communicator, comm size and proc name
CALL MPI_COMM_RANK(RP_COMM, myRPrank , iError)
CALL MPI_COMM_SIZE(RP_COMM, nRP_Procs, iError)

IF (myRPrank.EQ.0) WRITE(UNIT_stdOut,'(A,I0,A)') ' | RP COMM: ',nRP_Procs,' procs'

END SUBROUTINE InitRPCommunicator
#endif /*USE_MPI*/


!==================================================================================================================================
!> Read Recordpoint coordinates from HDF5 file
!==================================================================================================================================
SUBROUTINE ReadRPList(FileString)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_HDF5_Input
USE MOD_Mesh_Vars             ,ONLY:MeshFile,nGlobalElems
USE MOD_Mesh_Vars             ,ONLY:OffsetElem
USE MOD_Mesh_Vars             ,ONLY:nElems
USE MOD_RecordPoints_Vars     ,ONLY:RP_onProc,L_xi_RP,L_eta_RP,L_zeta_RP
USE MOD_RecordPoints_Vars     ,ONLY:d_L_xi_RP,d_L_eta_RP,d_L_zeta_RP,d_RP_ElemID
USE MOD_RecordPoints_Vars     ,ONLY:offsetRP,RP_ElemID,nRP,nGlobalRP
#if FV_ENABLED
USE MOD_RecordPoints_Vars     ,ONLY:FV_RP_ijk
#endif /*FV_ENABLED*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: FileString !< name of hdf5 file for readin of recordpoints
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: MeshFile_RPList
INTEGER                       :: nGlobalElems_RPList
INTEGER                       :: iElem,iRP,iRP_glob
INTEGER                       :: OffsetRPArray(2,nElems)
REAL,ALLOCATABLE              :: xi_RP(:,:)
!==================================================================================================================================

IF(MPIRoot)THEN
  IF(.NOT.FILEEXISTS(FileString))  CALL ABORT(__STAMP__, &
          'RPList from data file "'//TRIM(FileString)//'" does not exist')
END IF

SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' Read recordpoint definitions from data file "'//TRIM(FileString)//'" ...'

! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! compare mesh file names
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile_RPList)
IF(TRIM(MeshFile_RPList).NE.TRIM(MeshFile)) THEN
  ! Print empty line to break the ADVANCE=NO
  SWRITE(UNIT_stdOut,'(/,A,A,A)') ' WARNING: MeshFileName ',TRIM(MeshFile_RPList), &
                                  ' from RPList differs from Mesh File specified in parameterfile!'
END IF

! Readin OffsetRP
CALL GetDataSize(File_ID,'OffsetRP',nDims,HSize)
CHECKSAFEINT(HSize(2),4)
nGlobalElems_RPList=INT(HSize(2),4) !global number of elements
DEALLOCATE(HSize)
IF(nGlobalElems_RPList.NE.nGlobalElems) CALL ABORT(__STAMP__, &
          'nGlobalElems from RPList differs from nGlobalElems from Mesh File!')

CALL ReadArray('OffsetRP',2,(/2,nElems/),OffsetElem,2,IntArray=OffsetRPArray)

! Check if local domain contains any record points
! OffsetRP: first index: 1: offset in RP list for first RP on elem,
!                        2: offset in RP list for last RP on elem
! If these offsets are equal, no RP on elem.
nRP      = OffsetRPArray(2,nElems)-OffsetRPArray(1,1)
offsetRP = OffsetRPArray(1,1)

! Read in RP reference coordinates
CALL GetDataSize(File_ID,'xi_RP',nDims,HSize)
CHECKSAFEINT(HSize(2),4)
nGlobalRP=INT(HSize(2),4) !global number of RecordPoints
DEALLOCATE(HSize)

ALLOCATE(xi_RP(3,nRP))
CALL ReadArray('xi_RP',2,(/3,nRP/),offsetRP,2,RealArray=xi_RP)

IF(nRP.LT.1) THEN
  RP_onProc = .FALSE.
ELSE
  RP_onProc = .TRUE.
  ! create mapping to elements
  ALLOCATE(RP_ElemID(nRP))
!@cuf ALLOCATE(d_RP_ElemID(nRP))
  DO iRP=1,nRP
    iRP_glob = offsetRP+iRP
    DO iElem=1,nElems
      IF(iRP_glob.LE.OffsetRPArray(2,iElem) .AND. iRP_glob.GT.OffsetRPArray(1,iElem)) &
        RP_ElemID(iRP) = iElem
    END DO
  END DO
  d_RP_ElemID = RP_ElemID
END IF

CALL CloseDataFile()

IF (RP_onProc) THEN
  ALLOCATE(L_xi_RP  (0:PP_N,nRP))
  ALLOCATE(L_eta_RP (0:PP_N,nRP))
  ALLOCATE(L_zeta_RP(0:PP_N,nRP))
  CALL InitRPBasis(nRP,xi_RP,L_xi_RP,L_eta_RP,L_zeta_RP)

!@cuf ALLOCATE(d_L_xi_RP  (0:PP_N,nRP))
!@cuf ALLOCATE(d_L_eta_RP (0:PP_N,nRP))
!@cuf ALLOCATE(d_L_zeta_RP(0:PP_N,nRP))

  d_L_xi_RP   = L_xi_RP
  d_L_eta_RP  = L_eta_RP
  d_L_zeta_RP = L_zeta_RP

#if FV_ENABLED
  ALLOCATE(FV_RP_ijk(3,nRP))
  !=====================================================================================
  ! Two variants possible for FV:
  ! 1. The RP state is the nearest (reference-space) FV cells cell average:
  !    + Most simple solution
  !    - Spatial accuracy
  !    - Parameter space estimate can be wrong in case of strongly deformed meshes
  ! 2. The RP state is obtained by tri-linear interpolation from the 8 nearest
  !    (physical space) FV cell averages:
  !    + Probably gives the best quality
  !    - Requires geometry info in physical space
  !    - General implementation requires interpolation across macro-cell boundaries.
  !      Very difficult especially in an MPI setting.
  !
  ! We implement the first variant, the second may be an option if higher accuracy
  ! is desired, possibly with only element local interpolation.
  !=====================================================================================
  FV_RP_ijk = INT((xi_RP+1.)*0.5*(PP_N+1))
  FV_RP_ijk = MAX(FV_RP_ijk,0)
  FV_RP_ijk = MIN(FV_RP_ijk,PP_N)

#if PP_dim==2
  FV_RP_ijk(3,:) = 0
#endif /*PP_dim==2*/

#endif /*FV_ENABLED*/
END IF

DEALLOCATE(xi_RP)

SWRITE(UNIT_stdOut,'(A)',ADVANCE='YES')' DONE.'

END SUBROUTINE ReadRPList


!==================================================================================================================================
!> Precompute Lagrange basis function values at recordpoints
!==================================================================================================================================
SUBROUTINE InitRPBasis(nRP,xi_RP,L_xi_RP,L_eta_RP,L_zeta_RP)
! MODULES
USE MOD_PreProc
USE MOD_Interpolation_Vars    ,ONLY: xGP,wBary
USE MOD_Basis                 ,ONLY: LagrangeInterpolationPolys
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nRP                      !< size of recordpointarray
REAL,INTENT(IN)               :: xi_RP(3,nRP)             !< coordinates of recordpoints in reference space
REAL,INTENT(OUT)              :: L_xi_RP(  0:PP_N,nRP)    !< Lagrange basis evaluated at recordpoints (\f$\xi\f$-direction)
REAL,INTENT(OUT)              :: L_eta_RP( 0:PP_N,nRP)    !< Lagrange basis evaluated at recordpoints (\f$\eta\f$-direction)
REAL,INTENT(OUT)              :: L_zeta_RP(0:PP_N,nRP)    !< Lagrange basis evaluated at recordpoints (\f$\zeta\f$-direction)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: iRP
!==================================================================================================================================
! build local basis for Recordpoints
DO iRP=1,nRP
  CALL LagrangeInterpolationPolys(xi_RP(1,iRP),PP_N,xGP,wBary,L_xi_RP(  :,iRP))
  CALL LagrangeInterpolationPolys(xi_RP(2,iRP),PP_N,xGP,wBary,L_eta_RP( :,iRP))
#if PP_dim == 3
  CALL LagrangeInterpolationPolys(xi_RP(3,iRP),PP_N,xGP,wBary,L_zeta_RP(:,iRP))
#endif /*PP_dim==3*/
END DO

END SUBROUTINE InitRPBasis


!==================================================================================================================================
!> Evaluate solution at current time t at recordpoint positions and fill output buffer
!==================================================================================================================================
SUBROUTINE RecordPoints(nVar,StrVarNames,iter,t,forceSampling,streamID)
! MODULES
USE MOD_Globals
USE MOD_GPU
USE CUDAFOR
USE MOD_Preproc
USE MOD_Array_Operations, ONLY: ChangeSizeArray
USE MOD_Analyze_Vars,     ONLY: WriteData_dt,tWriteData
USE MOD_DG_Vars          ,ONLY: d_U
USE MOD_RecordPoints_Vars,ONLY: RP_Data,d_RP_ElemID
USE MOD_RecordPoints_Vars,ONLY: RP_BufferSize,RP_MaxBufferSize,RP_SamplingOffset,iSample
USE MOD_RecordPoints_Vars,ONLY: d_l_xi_RP,d_l_eta_RP,nRP
USE MOD_RecordPoints_Vars,ONLY: d_U_RP
USE MOD_Timedisc_Vars,    ONLY: dt
#if PP_dim==3
USE MOD_RecordPoints_Vars,ONLY: d_l_zeta_RP
#endif /*PP_dim==3*/
#if FV_ENABLED
USE MOD_FV_Vars          ,ONLY: FV_Elems
USE MOD_RecordPoints_Vars,ONLY: FV_RP_ijk
#endif /*FV_ENABLED*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVar                    !< Number of variables in U array
CHARACTER(LEN=255),INTENT(IN)  :: StrVarNames(nVar)       !< String with the names of the variables
INTEGER(KIND=8),INTENT(IN)     :: iter                    !< current number of timesteps
REAL,INTENT(IN)                :: t                       !< current time t
LOGICAL,INTENT(IN)             :: forceSampling           !< force sampling (e.g. at first/last timestep of computation)
INTEGER(KIND=CUDA_STREAM_KIND),OPTIONAL,INTENT(IN) :: streamID
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k,iRP
INTEGER                        :: NewSize
INTEGER(KIND=CUDA_STREAM_KIND) :: mystream
INTEGER                        :: N_loc
!----------------------------------------------------------------------------------------------------------------------------------

IF(MOD(iter,INT(RP_SamplingOffset,KIND=8)).NE.0 .AND. .NOT. forceSampling) RETURN

! U needs to be filled
iError = cudaDeviceSynchronize()

IF (iSample+1.GT.RP_BufferSize) THEN
  ! Grow the array if possible
  IF (RP_Buffersize.LT.RP_MaxBuffersize) THEN
    ! Compute required buffersize from timestep and add 20% tolerance
    ! +1 is added to ensure a minimum buffersize of 2
    NewSize = MIN(MAX(RP_BufferSize+1,FLOOR(REAL(RP_BufferSize)*1.2)),RP_MaxBuffersize)
    CALL ChangeSizeArray(RP_Data,RP_BufferSize,NewSize)
    RP_BufferSize = NewSize
    ! Force early writeout
  ELSE
    CALL WriteRP(nVar,StrVarNames,t)
  END IF
END IF ! iSample.GT.RP_BufferSize
iSample = iSample + 1

mystream=DefaultStream
IF (PRESENT(streamID)) mystream=streamID

! d_U_RP = 0.
iError = cudaMemsetAsync(d_U_RP,0.,nVar*nRP,mystream)
N_loc  = PP_N

!$cuf kernel do <<< *, *, 0, mystream >>>
DO iRP=1,nRP
  DO k=0,N_loc; DO j=0,N_loc; DO i=0,N_loc
    d_U_RP(:,iRP) = d_U_RP(:,iRP) + d_U(:,i,j,k,d_RP_ElemID(iRP))*d_l_xi_RP(i,iRP)*d_l_eta_RP(j,iRP)*d_l_zeta_RP(k,iRP)
  END DO; END DO; END DO ! k,j,i
END DO ! iRP

! Explicit wait
! iError = cudaStreamSynchronize(mystream)

! RP_Data(2:nVar+1,:,iSample) = d_U_RP
iError = cudaMemcpyAsync(RP_Data(2:nVar+1,:,iSample),d_U_RP,nVar*nRP,cudaMemcpyDeviceToHost,mystream)
RP_Data(1,       :,iSample) = t

END SUBROUTINE RecordPoints


!==================================================================================================================================
!> Writes the time history of the solution at the recordpoints to an HDF5 file
!==================================================================================================================================
SUBROUTINE WriteRP(nVar,StrVarNames,OutputTime,streamID)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_GPU
USE CUDAFOR
USE MOD_Array_Operations  ,ONLY: ChangeSizeArray
USE MOD_HDF5_Output       ,ONLY: WriteAttribute,GatheredWriteArray,MarkWriteSuccessfull
USE MOD_IO_HDF5           ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_Mesh_Vars         ,ONLY: MeshFile
USE MOD_Output_Vars       ,ONLY: ProjectName
USE MOD_Recordpoints_Vars ,ONLY: lastSample
USE MOD_Recordpoints_Vars ,ONLY: RPDefFile,RP_Data,iSample
USE MOD_Recordpoints_Vars ,ONLY: offsetRP,nRP,nGlobalRP
USE MOD_Recordpoints_Vars ,ONLY: RP_BufferSize,RP_Maxbuffersize
#if USE_MPI
USE MOD_Recordpoints_Vars ,ONLY: RP_COMM
USE MOD_RecordPoints_Vars ,ONLY: myRPrank
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: nVar                  !< Number of variables to write
CHARACTER(LEN=255),INTENT(IN)  :: StrVarNames(nVar)     !< String with the names of the variables
REAL,INTENT(IN)                :: OutputTime            !< time
INTEGER(KIND=CUDA_STREAM_KIND),OPTIONAL,INTENT(IN) :: streamID
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255)             :: tmp255
INTEGER(KIND=CUDA_STREAM_KIND) :: mystream
REAL                           :: startT,endT
!==================================================================================================================================

IF (iSample.LE.0) RETURN

FileName = TRIM(TIMESTAMP(TRIM(ProjectName)//'_RP',OutputTime))//'.h5'

#if USE_MPI
IF(myRPrank.EQ.0)THEN
#endif /*USE_MPI*/
  WRITE(UNIT_stdOut,'(A,I0,A,I0,A,I0,A)') ' RP Buffer  : ',iSample,', ',RP_BufferSize,', ', RP_MaxBuffersize, ' [used,alloc,avail]'
  WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')   ' WRITE RECORDPOINT DATA TO HDF5 FILE...'
  GETTIME(startT)

  CALL OpenDataFile(FileName,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

  ! Create dataset attributes
  CALL WriteAttribute(File_ID,'File_Type'  ,1,StrScalar=(/CHARACTER(LEN=255)::'RecordPoints_Data'/))
  tmp255=TRIM(MeshFile)
  CALL WriteAttribute(File_ID,'MeshFile'   ,1,StrScalar=(/tmp255/))
  tmp255=TRIM(ProjectName)
  CALL WriteAttribute(File_ID,'ProjectName',1,StrScalar=(/tmp255/))
  tmp255=TRIM(RPDefFile)
  CALL WriteAttribute(File_ID,'RPDefFile'  ,1,StrScalar=(/tmp255/))
  CALL WriteAttribute(File_ID,'VarNames'   ,nVar,StrArray=StrVarNames)
  CALL WriteAttribute(File_ID,'Time'       ,1,RealScalar=OutputTime)

  CALL CloseDataFile()
#if USE_MPI
END IF
#endif /*USE_MPI*/

! Ensure RP_Data is up-to-date
mystream=DefaultStream
IF (PRESENT(streamID)) mystream=streamID
iError = cudaStreamSynchronize(mystream)

CALL GatheredWriteArray(FileName                                       ,&
                        create      = .FALSE.                          ,&
                        DataSetName = 'RP_Data'                        ,&
                        rank        = 3                                ,&
                        nValGlobal  = (/nVar+1 ,nGlobalRP,iSample   /) ,&
                        nVal        = (/nVar+1 ,nRP      ,iSample   /) ,&
                        offset      = (/0      ,offsetRP ,0         /) ,&
                        collective  = .TRUE.                           ,&
                        RealArray   = RP_Data(:,:,1:iSample))

! Store the last sample for the next RP file
lastSample = RP_Data(:,:,iSample)

! Reset buffer
RP_Data = 0.
iSample = 0

! Shrink the array
CALL ChangeSizeArray(RP_Data,RP_BufferSize,1)
RP_BufferSize = 1

#if USE_MPI
IF(myRPrank.EQ.0)THEN
#endif /*USE_MPI*/
  CALL MarkWriteSuccessfull(FileName)
  GETTIME(EndT)
  WRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' DONE  [',EndT-StartT,'s]'
#if USE_MPI
END IF
#endif /*USE_MPI*/

END SUBROUTINE WriteRP


!==================================================================================================================================
!> Deallocate recordpoint arrays
!==================================================================================================================================
SUBROUTINE FinalizeRecordPoints()
! MODULES
USE MOD_Globals,                 ONLY: iError
USE MOD_RecordPoints_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================

SDEALLOCATE(RP_Data)
!@cuf SDEALLOCATE(d_U_RP)
SDEALLOCATE(RP_ElemID)
!@cuf SDEALLOCATE(d_RP_ElemID)
SDEALLOCATE(L_xi_RP)
SDEALLOCATE(L_eta_RP)
SDEALLOCATE(L_zeta_RP)
!@cuf SDEALLOCATE(d_L_xi_RP)
!@cuf SDEALLOCATE(d_L_eta_RP)
!@cuf SDEALLOCATE(d_L_zeta_RP)
SDEALLOCATE(lastSample)
#if FV_ENABLED
SDEALLOCATE(FV_RP_ijk)
#endif /*FV_ENABLED*/

#if USE_MPI
! Free MPI communicator
IF(RP_COMM.NE.MPI_COMM_NULL) CALL MPI_COMM_FREE(RP_COMM, iError)
#endif /*USE_MPI*/

RecordPointsInitIsDone = .FALSE.

END SUBROUTINE FinalizeRecordPoints

END MODULE MOD_RecordPoints
