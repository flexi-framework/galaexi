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
!> Module containing the main procedures for the visu tool: visu_requestInformation is called by ParaView to create a
!> list of available variables and visu is the main routine which is either called by ParaView to get the data it visualizes
!> or by the standalone tool.
!===================================================================================================================================
MODULE MOD_Visu
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE visu_getVarNamesAndFileType
  MODULE PROCEDURE visu_getVarNamesAndFileType
END INTERFACE

INTERFACE Visu_InitFile
  MODULE PROCEDURE Visu_InitFile
END INTERFACE

INTERFACE visu
  MODULE PROCEDURE visu
END INTERFACE

INTERFACE FinalizeVisu
  MODULE PROCEDURE FinalizeVisu
END INTERFACE

PUBLIC:: visu_getVarNamesAndFileType
PUBLIC:: visu_InitFile
PUBLIC:: visu
PUBLIC:: FinalizeVisu

CONTAINS

!===================================================================================================================================
!> Create a list of available variables for ParaView. This list contains the conservative, primitive and derived quantities
!> that are available in the current equation system as well as the additional variables read from the state file.
!> The additional variables are stored in the datasets 'ElemData' (elementwise data) and 'FieldData' (pointwise data).
!> Also a list of all available boundary names is created for surface visualization.
!===================================================================================================================================
SUBROUTINE visu_getVarNamesAndFileType(statefile,meshfile,varnames_loc, bcnames_loc)
USE MOD_Globals
USE MOD_EOS_Posti_Vars ,ONLY: DepNames,nVarDepEOS
USE MOD_IO_HDF5        ,ONLY: GetDatasetNamesInGroup,File_ID
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,GetDataSize,GetVarNames,ISVALIDMESHFILE,ISVALIDHDF5FILE,ReadAttribute
USE MOD_HDF5_Input     ,ONLY: DatasetExists,HSize,nDims,ReadArray
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_Restart        ,ONLY: InitRestartFile
USE MOD_Restart_Vars   ,ONLY: RestartMode
USE MOD_Visu_Vars      ,ONLY: FileType,VarNamesHDF5,nBCNamesAll,nVarIni,nVar_State,IJK_exists
! IMPLICIT VARIABLE HANDLINGs
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)                       :: statefile
CHARACTER(LEN=*)  ,INTENT(IN)                       :: meshfile
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: varnames_loc(:)
CHARACTER(LEN=255),INTENT(INOUT),ALLOCATABLE,TARGET :: bcnames_loc(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                             :: i,j,nVar,dims
LOGICAL                                             :: varnames_found,readDGsolutionVars,sameVars,VarNamesExist,file_exists
CHARACTER(LEN=255),ALLOCATABLE                      :: datasetNames(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: varnames_tmp(:)
CHARACTER(LEN=255),ALLOCATABLE                      :: tmp(:)
CHARACTER(LEN=255)                                  :: MeshFile_loc
INTEGER                                             :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================

IF (ISVALIDMESHFILE(statefile)) THEN      ! MESH
  SDEALLOCATE(varnames_loc)

  ! IJK-sorted mesh
  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL DatasetExists(File_ID,'Elem_IJK',IJK_exists)
  IF (IJK_exists) THEN
    ALLOCATE(varnames_loc(5))
  ELSE
    ALLOCATE(varnames_loc(2))
  END IF

  varnames_loc(1) = 'ScaledJacobian'
  varnames_loc(2) = 'ScaledJacobianElem'
  IF (IJK_exists) THEN
    varnames_loc(3) = 'Elem_I'
    varnames_loc(4) = 'Elem_J'
    varnames_loc(5) = 'Elem_K'
  END IF

  FileType='Mesh'

ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! other file
  SDEALLOCATE(varnames_loc)
  CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL ReadAttribute(File_ID,'File_Type',   1,StrScalar =FileType)

  SELECT CASE(TRIM(FileType))
    ! check if variables in state file are the same as in the EQNSYS and set FileType to 'Generic' if not
    CASE('State')
      SDEALLOCATE(VarNamesHDF5)
      CALL GetVarNames("VarNames",VarNamesHDF5,VarNamesExist)

      sameVars = .FALSE.
      IF (VarNamesExist .AND. PP_nVar.EQ.SIZE(VarNamesHDF5)) THEN
        sameVars = .TRUE.
        DO i = 1,SIZE(VarNamesHDF5)
          sameVars = sameVars.AND.(STRICMP(VarNamesHDF5(i),DepNames(i)))
        END DO
      END IF
      IF (.NOT.sameVars) FileType = 'Generic'

    CASE('TimeAvg')
      IF (nVarIni.EQ.0) THEN
        FileType = 'Generic'
      ELSE
        CALL CloseDataFile()
        ! This routine requires the file to be closed
        CALL InitRestartFile(statefile)
        CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
        SELECT CASE(RestartMode)
          CASE(0)
            SDEALLOCATE(VarNamesHDF5)
            CALL GetVarNames("VarNames_Mean",VarNamesHDF5,VarNamesExist)
            FileType = 'Generic'
          CASE(2,3)
            SDEALLOCATE(VarNamesHDF5)
            CALL GetVarNames("VarNames_Mean",VarNamesHDF5,VarNamesExist)
            FileType = 'State'
            ! When restarting from a time-averaged file, we convert to U array to PP_nVar
            nVar_State = PP_nVar
          CASE DEFAULT
            FileType = 'Generic'
        END SELECT
      END IF
  END SELECT

  IF (STRICMP(FileType,'State')) THEN
    nVar = nVarDepEOS
    ALLOCATE(varnames_loc(nVar))
    varnames_loc(1:nVar) = DepNames
    readDGsolutionVars   = .FALSE.
  ELSE
    nVar=0
    readDGsolutionVars   = .TRUE.
  END IF

  CALL GetDatasetNamesInGroup("/",datasetNames)

  DO i = 1,SIZE(datasetNames)
    SDEALLOCATE(varnames_tmp)
    VarNamesExist=.FALSE.
    CALL DatasetExists(File_ID,"VarNames_"//TRIM(datasetNames(i)),varnames_found,attrib=.TRUE.)
    IF (varnames_found) THEN
      CALL GetVarNames("VarNames_"//TRIM(datasetNames(i)),varnames_tmp,VarNamesExist)
    ELSE
      IF (STRICMP(datasetNames(i), "DG_Solution")) THEN
        IF (readDGsolutionVars) THEN
          CALL GetVarNames("VarNames",varnames_tmp,VarNamesExist)
        END IF
      ! ELSEIF (RestartMode.GT.1 .AND. STRICMP(datasetNames(i), "Mean")) THEN
      !   IF (readDGsolutionVars) THEN
      !     CALL GetVarNames("VarNames_Mean",varnames_tmp,VarNamesExist)
      !   END IF
      ELSE IF(STRICMP(datasetNames(i), "ElemData")) THEN
        CALL GetVarNames("VarNamesAdd",varnames_tmp,VarNamesExist)
      ELSE IF(STRICMP(datasetNames(i), "FieldData")) THEN
        CALL GetVarNames("VarNamesAddField",varnames_tmp,VarNamesExist)
      ELSE
        CALL GetDataSize(File_ID,TRIM(datasetNames(i)),dims,HSize)
        IF ((dims.NE.5).AND.(dims.NE.2)) CYCLE ! Do not add datasets to the list that can not contain elementwise or field data
        ALLOCATE(varnames_tmp(INT(HSize(1))))
        DO j=1,INT(HSize(1))
          WRITE(varnames_tmp(j),'(I0)') j
        END DO
        VarNamesExist=.TRUE.
      END IF
    END IF
    IF (.NOT.VarNamesExist) CYCLE

    ! increase array 'varnames_loc'
    IF (nVar.GT.0) THEN
      ALLOCATE(tmp(nVar))
      tmp = varnames_loc
    END IF
    SDEALLOCATE(varnames_loc)
    ALLOCATE(varnames_loc(nVar+SIZE(varnames_tmp)))
    IF (nVar.GT.0) varnames_loc(1:nVar) = tmp(1:nVar)
    SDEALLOCATE(tmp)

    ! copy new varnames from varnames_tmp to varnames_loc
    DO j=1,SIZE(varnames_tmp)
      varnames_loc(nVar+j) = TRIM(datasetNames(i))//":"//TRIM(varnames_tmp(j))
    END DO
    nVar = nVar + SIZE(varnames_tmp)
  END DO ! i = 1,SIZE(datasetNames)

  IF (LEN_TRIM(meshfile).EQ.0) THEN
    ! Save mesh file to get boundary names later
    CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar =MeshFile_loc)
  ELSE
    MeshFile_loc = meshfile
  END IF

  CALL CloseDataFile()

  INQUIRE(FILE=TRIM(MeshFile_loc), EXIST=file_exists)
  IF (file_exists) THEN
    ! Open the mesh file and read all boundary names for surface visualization
    CALL OpenDataFile(MeshFile_loc,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
    CHECKSAFEINT(HSize(1),4)
    nBCNamesAll=INT(HSize(1),4)
    DEALLOCATE(HSize)
    SDEALLOCATE(bcnames_loc)
    ALLOCATE(bcnames_loc(nBCNamesAll))
    CALL ReadArray('BCNames',1,(/nBCNamesAll/),Offset,1,StrArray=bcnames_loc)
    CALL CloseDataFile()
  END IF

  SDEALLOCATE(datasetNames)
END IF

END SUBROUTINE visu_getVarNamesAndFileType

!===================================================================================================================================
!> This routine is used to prepare everything we need to visualize data from a statefile.
!> This includes:
!> * Get the mesh file
!> * Read the desired visualization polynomial degree, the visualization dimennsion, the node type we want to visualize on and the
!>   Dg only option
!> * Decide whether the state file, the mesh file, the visualization polynomial degree or the dg only option changed. This is
!>   needed to decide what parts of the visualization routines should be called.
!> * Call routines that build the distribution between FV and DG elements and the mappings needed to calculate and visualize the
!>   desired variables.
!===================================================================================================================================
SUBROUTINE visu_InitFile(statefile,postifile)
! MODULES
USE HDF5
USE MOD_Preproc
USE MOD_Globals
USE MOD_Visu_Vars
USE MOD_EOS_Posti_Vars
USE MOD_MPI                ,ONLY: InitMPI
USE MOD_HDF5_Input         ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,GetArrayAndName
USE MOD_HDF5_Input         ,ONLY: ReadAttribute,File_ID,OpenDataFile,GetDataProps,CloseDataFile,ReadArray,DatasetExists
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_Output_Vars        ,ONLY: ProjectName
USE MOD_StringTools        ,ONLY: STRICMP,INTTOSTR
USE MOD_ReadInTools        ,ONLY: prms,GETINT,GETLOGICAL,addStrListEntry,GETSTR,FinalizeParameters,CountOption
USE MOD_Posti_Mappings     ,ONLY: Build_FV_DG_distribution,Build_mapDepToCalc_mapAllVarsToVisuVars
USE MOD_Visu_Avg2D         ,ONLY: InitAverage2D,BuildVandermonds_Avg2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)    :: statefile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL
LOGICAL                          :: RestartMean
CHARACTER(LEN=255)               :: NodeType_State, cwd
INTEGER                          :: nElems_State
#if PP_N!=N
INTEGER                          :: N_State
#endif
!===================================================================================================================================
IF (STRICMP(fileType,'Mesh')) THEN
    CALL CollectiveStop(__STAMP__, &
        "FileType==Mesh, but we try to initialize a state file!")
END IF

! open state file to be able to read attributes
CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! read the meshfile attribute from statefile
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar =MeshFile_state)

! read options from posti parameter file
CALL prms%read_options(postifile)

! Get number of variables to be visualized
nVarIni = CountOption("VarName")

! get properties
#if EQNSYSNR != 1
! Check if the file is a time-averaged file
CALL DatasetExists(File_ID,'Mean',RestartMean)
! Read in attributes
IF (.NOT.RestartMean .OR. nVarIni.EQ.0) THEN
#endif /* EQNSYSNR != 1 */
#if PP_N==N
  CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State)
#else
  CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State)
#endif /* PP_N==N */
#if EQNSYSNR != 1
! Check if the file is a time-averaged file
ELSE
#if PP_N==N
  CALL GetDataProps(nVar_State,PP_N,nElems_State,NodeType_State,'Mean')
#else
  CALL GetDataProps(nVar_State,N_State,nElems_State,NodeType_State,'Mean')
#endif /* PP_N==N */
END IF
#endif /* EQNSYSNR != 1 */

! read options from posti parameter file
NVisu             = GETINT("NVisu",INTTOSTR(PP_N))
HighOrder         = GETLOGICAL('HighOrder')

! again read MeshFile from posti prm file (this overwrites the MeshFile read from the state file)
Meshfile          =  GETSTR("MeshFile",MeshFile_state)
IF (.NOT.FILEEXISTS(MeshFile) .OR. ((Meshfile(1:1) .NE. "/") .OR. (Meshfile(1:1) .NE. "~") .OR. (Meshfile(1:1) .NE. "."))) THEN
  !!!!!!
  ! WARNING: GETCWD is a GNU extension to the Fortran standard and will probably not work on other compilers
  CALL GETCWD(cwd)
  !!!!!!
  Meshfile          =  TRIM(cwd) // "/" // TRIM(Meshfile)
END IF
Avg2D             = GETLOGICAL("Avg2D")
#if PP_dim == 2
IF (Avg2D) THEN
  CALL PrintWarning("Avg2D not available for 2D-Posti! Switching it OFF.")
  Avg2D = .FALSE.
END IF
#endif
NodeTypeVisuPosti = GETSTR('NodeTypeVisu')
DGonly            = GETLOGICAL('DGonly')
CALL CloseDataFile()

CALL visu_getVarNamesAndFileType(statefile,'',VarnamesAll,BCNamesAll)

CALL OpenDataFile(statefile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

! HDF5 Output for avg2D
IF (Avg2D) Avg2DHDF5Output = GETLOGICAL("Avg2DHDF5Output")

! check if state, mesh, NVisu, DGonly or Avg2D changed
changedStateFile = .NOT.STRICMP(statefile,statefile_old)
changedMeshFile  = .NOT.(STRICMP(MeshFile,MeshFile_old))
changedDGonly    = (DGonly.NEQV.DGonly_old)
changedAvg2D     = (Avg2D.NEQV.Avg2D_old)

SWRITE(*,*) "state file old -> new: ", TRIM(statefile_old), " -> ",TRIM(statefile)
SWRITE(*,*) " mesh file old -> new: ", TRIM(MeshFile_old) , " -> ",TRIM(MeshFile)

! if Mesh or State changed readin some more attributes/parameters
IF (changedStateFile.OR.changedMeshFile) THEN
  IF (.NOT.STRICMP(NodeType_State, NodeType)) THEN
    CALL CollectiveStop(__STAMP__, &
        "NodeType of state does not match with NodeType the visu-posti is compiled with!")
  END IF
  CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =ProjectName)
  CALL ReadAttribute(File_ID,'Time',        1,RealScalar=OutputTime)
  ! If the polynomial degree is changing, we could need new mesh mappings.
  IF (NState_old.NE.PP_N) changedMeshFile = .TRUE.
END IF

CALL CloseDataFile()

! Polynomial degree for calculations
NCalc             = GETINT("NCalc",INTTOSTR(PP_N))
IF (NCalc.LE.0) NCalc = PP_N
changedNCalc      = NCalc.NE.NCalc_old

! Output of averaged data is only available for NVisu = PP_N and NodeTypeVisuPosti=NodeType_State
! These settings are enforced here!
IF (Avg2DHDF5Output) THEN
        NVisu = PP_N
        NodeTypeVisuPosti=NodeType_State
END IF
! Check for changed visualization basis here to take change done for average output into account
changedNVisu     = ((NVisu.NE.NVisu_old) .OR. (NodeTypeVisuPosti.NE.NodeTypeVisuPosti_old))

! set number of dependent and raw variables
SDEALLOCATE(DepTable)
SDEALLOCATE(DepSurfaceOnly)
SDEALLOCATE(DepVolumeOnly)
nVarAll=SIZE(VarnamesAll)
IF (STRICMP(FileType,'State')) THEN
  StateFileMode = .TRUE.
  nVarDep = nVarDepEOS
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  ALLOCATE(DepSurfaceOnly(nVarDep))
  ALLOCATE(DepVolumeOnly(nVarDep))
  DepTable = DepTableEOS
  DepSurfaceOnly = DepSurfaceOnlyEOS
  DepVolumeOnly  = DepVolumeOnlyEOS
ELSE
  StateFileMode = .FALSE.
  nVarDep = 0
  ALLOCATE(DepTable(nVarDep,0:nVarDep))
  ALLOCATE(DepSurfaceOnly(nVarDep))
  ALLOCATE(DepVolumeOnly(nVarDep))
  DepTable = 0
  DepSurfaceOnly = 0
  DepVolumeOnly  = 0
END IF

! build distribution of FV and DG elements, which is stored in FV_Elems_loc
IF (changedStateFile.OR.changedMeshFile.OR.changedDGonly) THEN
  CALL Build_FV_DG_distribution(&
#if FV_ENABLED
    statefile&
#endif
    )
END IF

! reset withDGOperator flag and check if it is needed due to existing FV elements
withDGOperator = .FALSE.
#if FV_RECONSTRUCT
! If what we want to visualize is a state and has FV elements, the DG operator needs to be called for reconstruction
IF (StateFileMode) THEN
  IF (hasFV_Elems) withDGOperator = .TRUE.
END IF
#endif

! build mappings of variables which must be calculated/visualized
! also set withDGOperator flag if a dependent variable requires the evaluation of the DG operator
CALL Build_mapDepToCalc_mapAllVarsToVisuVars()

IF (Avg2D) THEN
  CALL InitAverage2D()
  CALL BuildVandermonds_Avg2D(NCalc&
#if FV_ENABLED
    ,NCalc_FV&
#endif
    )
END IF

changedWithDGOperator = (withDGOperator.NEQV.withDGOperator_old)
END SUBROUTINE visu_InitFile

!===================================================================================================================================
!> Main routine of the visualization tool visu. Called either by the ParaView plugin or by the standalone program version.
!===================================================================================================================================
SUBROUTINE visu(mpi_comm_IN, prmfile, postifile, statefile)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Visu_Vars
USE MOD_HDF5_Input          ,ONLY: ISVALIDMESHFILE,ISVALIDHDF5FILE,OpenDataFile,CloseDataFile
USE MOD_Interpolation_Vars  ,ONLY: NodeType,NodeTypeVISUFVEqui
USE MOD_IO_HDF5             ,ONLY: InitMPIInfo
USE MOD_MPI                 ,ONLY: InitMPI
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_DG,CalcSurfQuantities_DG
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_DG,ConvertToSurfVisu_DG,ConvertToVisu_GenericData
USE MOD_Posti_ReadState     ,ONLY: ReadState
USE MOD_Posti_Mappings      ,ONLY: Build_mapBCSides
USE MOD_Posti_VisuMesh      ,ONLY: BuildVisuCoords,BuildSurfVisuCoords
USE MOD_Posti_VisuMesh      ,ONLY: VisualizeMesh
USE MOD_ReadInTools         ,ONLY: prms,FinalizeParameters,ExtractParameterFile,PrintDefaultParameterFile
USE MOD_Restart_Vars        ,ONLY: RestartMode
USE MOD_StringTools         ,ONLY: STRICMP,set_formatting,clear_formatting
USE MOD_Visu_Avg2D          ,ONLY: Average2D,WriteAverageToHDF5
#if FV_ENABLED
USE MOD_Posti_Calc          ,ONLY: CalcQuantities_FV,CalcSurfQuantities_FV
USE MOD_Posti_ConvertToVisu ,ONLY: ConvertToVisu_FV,ConvertToSurfVisu_FV
#endif /*FV_ENABLED*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)               :: mpi_comm_IN
CHARACTER(LEN=255),INTENT(INOUT) :: prmfile
CHARACTER(LEN=255),INTENT(INOUT) :: postifile
CHARACTER(LEN=255),INTENT(IN)    :: statefile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                          :: changedPrmFile
!===================================================================================================================================

!**********************************************************************************************
! General workflow / principles of the visu ParaView-plugin
!
! * all arrays are SDEALLOCATEd just before they are allocated. This is done to keep their
!   content during successive calls of the visu during a ParaView session. They are only
!   deallocated and reallocated, if there content should change. For example the coords of the
!   mesh file only change if the mesh, NVisu or the distribution of DG/FV elements changes.
!
! VISUALIZE MESH: Call 'VisualizeMesh' routine, which reads the mesh, interpolates it to
!   the visu grid and writes it to VTK 3D arrays.
!
! VISUALIZE STATE:
! * There are two different modes:
!   - without gradients: All quantity that should be visualized can be computed without any
!                        gradient computation (mainly the conservative and primitive quantities)
!                        If FV_RECONSTRUCT is enabled and there are FV elements present in
!                        the state, then this mode is not available.
!                        U is read from the state file directly and the EOS is initialized to
!                        perform ConsToPrim conversions based only on the conservative state U.
!
!   - with gradiens: There are quantities that require the computation of gradients or there
!                    are FV elements with FV_RECONSTRUCT enabled. In this case the DG operator
!                    'DGTimeDerivative_weakForm' is called once to fill the gradients and the
!                    reconstruction of the FV subcell method.
!                    This requires the initialization of several modules of the FLEXI.
!                    U is read via a call of 'Restart'. In the DGTimeDerivative_weakForm the
!                    primitive quantities U_Prim and gradUx/y/z as well as gradUxi/eta/zeta are
!                    filled. These are used to calculate the visu-quantities.
!
! * The calculation of derived quantities is performed on a arbitrary polynomial degree
!   NCalc and afterwards interpolated to NVisu. Default is PP_N.
!   This is not the case for FV elements with FV_RECONSTRUCT enabled.
!   These require to reconstruct the solution first to the visu grid and afterwards can
!   calculate the derived quantities on the NVisu_FV grid.
!
! * The dependencies of the visu-quantities on the state-quantities is stored in a dependency
!   integer table 'DepTable' and it corresponding row/col names 'DepNames' (see eos_vars.f90).
!   The string vector 'DepNames' contains all available quantities available for visualization
!   and the 'DepTable' contains in a row the dependencies of a quantity on other quantities.
!
! * The calculation of the visu-quantities is done in two steps:
!   1. calculate all quantities needed for the visu-quantities (stored in UCalc)
!   2. pick the visu-quantities from UCalc and interpolate them to NVisu (stored in UVisu)
!
! * Therefore two mappings from all available quantities to the calc-quantities and the
!   visu-quantities exist:
!   - 'mapAllVarsToVisuVars' is a integer array of size (1:nVarAll), where nVarAll is the total amount
!     of available quantities. This map contains a zero for all not-to-visu-quantities and for
!     all quantities the index where it is stored in 'UVisu'.
!     This mapping is filled from the 'VarName' entries in the parameter file.
!   - 'mapDepToCalc' is the same as mapAllVarsToVisuVars, but for all (intermediate) quantities stored in 'UCalc'.
!     This mapping is filled from the DepTable.
!
! CHANGED system:
! * There are different logical changedXXX variables, which indicated if XXX changed during
!   successive calls of the visu. These variables control the general workflow of the visu.
!   - changedStateFile:     new state file
!   - changedMeshFile:      new mesh file (only possible if changedStateFile==TRUE)
!   - changedVarNames:      new set of variables to visualize
!   - changedNVisu:         new NVisu, new Nodetype
!   - changedFV_Elems:      new distribution of FV/DG elements (only if changedStateFile==TRUE)
!   - changedWithDGOperator: different mode, with/without gradients
!   - changedDGonly:        the visualization of FV elements as DG elements was set or unset
!   - changedNCalc:         the polynomial degree used for calculations changed
!
! WORKFLOW:
! * The main steps are:
!   1. call InitFile for FV/DG distribution and mappings
!   2. read solution          (if changedStateFile or changedWithDGOperator or changedDGonly)
!   3. build mapping for BC sides that should be visualized (done after read solution since some
!      mesh infos are needed)
!   4. read Mesh              (if changedMeshFile)
!   5. compute UCalc          (if changedStateFile or changedVarNames or changedDGonly or changedNCalc)
!   6. convert to UVisu       (if changedStateFile or changedVarNames or changedNVisu or changedDGonly or changedNCalc)
!   7. build visu mesh        (if changedMeshFile  or changedNVisu or changedFV_Elems or changedDGonly)
!   5. - 7. are done seperately for surface variables if surface visualization is turned on
!   8. write VTK arrays       (always!)
!
!**********************************************************************************************

CALL SetStackSizeUnlimited()
postiMode = .TRUE. ! Flag used in FLEXI routines to do things only for POSTI usage
CALL InitMPI(mpi_comm_IN)
CALL InitMPIInfo()

CALL FinalizeParameters()
! Read Varnames to visualize and build calc and visu dependencies
CALL prms%SetSection("posti")
CALL prms%CreateStringOption( "MeshFile"        , "Custom mesh file ")
CALL prms%CreateStringOption( "VarName"         , "Names of variables, which should be visualized.", multiple=.TRUE.)
CALL prms%CreateLogicalOption("noVisuVars"      , "If no VarNames are given, this flags supresses visu of standard variables",&
                                                  ".FALSE.")
CALL prms%CreateIntOption(    "NVisu"           , "Polynomial degree at which solution is sampled for visualization.")
CALL prms%CreateIntOption(    "NCalc"           , "Polynomial degree at which calculations are done.")
CALL prms%CreateLogicalOption("Avg2D"           , "Average solution in z-direction",".FALSE.")
CALL prms%CreateLogicalOption("Avg2DHDF5Output" , "Write averaged solution to HDF5 file",".FALSE.")
CALL prms%CreateStringOption( "NodeTypeVisu"    , "NodeType for visualization. Visu, Gauss,Gauss-Lobatto,Visu_inner"    ,"VISU")
CALL prms%CreateLogicalOption("DGonly"          , "Visualize FV elements as DG elements."    ,".FALSE.")
CALL prms%CreateStringOption( "BoundaryName"    , "Names of boundaries for surfaces, which should be visualized.", multiple=.TRUE.)
CALL prms%CreateLogicalOption("HighOrder"       , "Write high-order element representation",".FALSE.")

IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2,statefile) !statefile string conatains --help etc!
  STOP
END IF

SWRITE(UNIT_stdOut,'(A,A)') " READING FROM: ", TRIM(statefile)

changedStateFile      = .FALSE.
changedMeshFile       = .FALSE.
changedNVisu          = .FALSE.
changedNCalc          = .FALSE.
changedVarNames       = .FALSE.
changedFV_Elems       = .FALSE.
changedWithDGOperator = .FALSE.
changedDGonly         = .FALSE.

IF (ISVALIDMESHFILE(statefile)) THEN ! visualize mesh
  SWRITE(UNIT_stdOut,'(A3,A30,A3,A33,A13)')' | ','                   Mode ',' | ','Mesh',' | HDF5    |'
  MeshFileMode = .TRUE.
  MeshFile      = statefile
  nVar_State    = 0
  withDGOperator = .FALSE.
  doSurfVisu     = .FALSE.
  CALL visu_getVarNamesAndFileType(MeshFile,'',VarNamesAll,BCNamesAll)
  CALL VisualizeMesh(postifile,MeshFile)
ELSE IF (ISVALIDHDF5FILE(statefile)) THEN ! visualize state file
  SWRITE(UNIT_stdOut,'(A3,A30,A3,A33,A13)')' | ','                   Mode ',' | ','State',' | HDF5    |'
  MeshFileMode = .FALSE.
  ! initialize state file
  CALL visu_InitFile(statefile,postifile)

  ! read solution from state file (either direct or including a evaluation of the DG operator)
  IF (LEN_TRIM(prmfile).EQ.0) THEN
    changedPrmFile = .NOT.STRICMP(prmfile_old, ".flexi.ini")
  ELSE
    changedPrmFile = (prmfile .NE. prmfile_old)
  END IF

  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') " DETECTING REQUIRED PARAMETERS..."
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | doSurfVisu              "
  CALL set_formatting(MERGE("blue ","green",doSurfVisu))             ; SWRITE(UNIT_stdOut,'(L1)') doSurfVisu             ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedStateFile        "
  CALL set_formatting(MERGE("blue ","green",changedStateFile))       ; SWRITE(UNIT_stdOut,'(L1)') changedStateFile       ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedMeshFile         "
  CALL set_formatting(MERGE("blue ","green",changedMeshFile))        ; SWRITE(UNIT_stdOut,'(L1)') changedMeshFile        ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedNVisu            "
  CALL set_formatting(MERGE("blue ","green",changedNVisu))           ; SWRITE(UNIT_stdOut,'(L1)') changedNVisu           ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedNCalc            "
  CALL set_formatting(MERGE("blue ","green",changedNCalc))           ; SWRITE(UNIT_stdOut,'(L1)') changedNCalc           ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedVarNames         "
  CALL set_formatting(MERGE("blue ","green",changedVarNames))        ; SWRITE(UNIT_stdOut,'(L1)') changedVarNames        ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedFV_Elems         "
  CALL set_formatting(MERGE("blue ","green",changedFV_Elems))        ; SWRITE(UNIT_stdOut,'(L1)') changedFV_Elems        ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedWithDGOperator   "
  CALL set_formatting(MERGE("blue ","green",changedWithDGOperator))  ; SWRITE(UNIT_stdOut,'(L1)') changedWithDGOperator  ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedDGonly           "
  CALL set_formatting(MERGE("blue ","green",changedDGonly))          ; SWRITE(UNIT_stdOut,'(L1)') changedDGonly          ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedAvg2D            "
  CALL set_formatting(MERGE("blue ","green",changedAvg2D))           ; SWRITE(UNIT_stdOut,'(L1)') changedAvg2D           ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedPrmFile          "
  CALL set_formatting(MERGE("blue ","green",changedPrmFile))         ; SWRITE(UNIT_stdOut,'(L1)') changedPrmFile         ; CALL clear_formatting()
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') " | changedBCNames          "
  CALL set_formatting(MERGE("blue ","green",changedBCNames))         ; SWRITE(UNIT_stdOut,'(L1)') changedBCNames         ; CALL clear_formatting()
  SWRITE(UNIT_StdOut,'(132("-"))')

  IF (changedStateFile.OR.changedWithDGOperator.OR.changedPrmFile.OR.changedDGonly) THEN
    CALL ReadState(prmfile,statefile)
  END IF

  ! build mappings of BC sides for surface visualization
  CALL Build_mapBCSides()

  ! ===== calc solution =====
  IF (changedStateFile.OR.changedVarNames.OR.changedDGonly.OR.changedNCalc) THEN
    CALL CalcQuantities_DG()
#if FV_ENABLED
    CALL CalcQuantities_FV()
#endif
  END IF
  IF (doSurfVisu) THEN
    ! calc surface solution
    IF (changedStateFile.OR.changedVarNames.OR.changedDGonly.OR.changedNCalc.OR.changedBCnames) THEN
      CALL CalcSurfQuantities_DG()
#if FV_ENABLED
      CALL CalcSurfQuantities_FV()
#endif
    END IF
  END IF

  ! ===== convert solution to visu grid =====
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedNCalc.OR.changedAvg2D) THEN
    ! ===== Avg2d =====
    IF (Avg2d) THEN
      SDEALLOCATE(UVisu_DG)
      SDEALLOCATE(UVisu_FV)
      ALLOCATE(UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,nElemsAvg2D_DG,nVarVisu))
      ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV,nVarVisu))
      CALL Average2D(nVarCalc,nVarCalc_FV,NCalc,NCalc_FV,nElems_DG,nElems_FV,NodeType,UCalc_DG,UCalc_FV,&
          Vdm_DGToFV,Vdm_FVToDG,Vdm_DGToVisu,Vdm_FVToVisu,1,nVarDep,mapDepToCalc,&
          UVisu_DG,UVisu_FV)
    ELSE
      CALL ConvertToVisu_DG()
#if FV_ENABLED
      CALL ConvertToVisu_FV()
#endif
    END IF
  END IF
  IF (doSurfVisu) THEN
    ! convert Surface DG solution to visu grid
    IF (changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedNCalc.OR.changedBCnames) THEN
      CALL ConvertToSurfVisu_DG()
#if FV_ENABLED
      CALL ConvertToSurfVisu_FV()
#endif
    END IF
  END IF

  ! convert generic data to visu grid
  IF (changedStateFile.OR.changedVarNames.OR.changedNVisu.OR.changedDGonly.OR.changedBCnames.OR.changedAvg2D) THEN
    CALL ConvertToVisu_GenericData(statefile)
  END IF

  IF (Avg2DHDF5Output) CALL WriteAverageToHDF5(nVarVisu,NVisu,NodeType,OutputTime,MeshFile_state,UVisu_DG&
#if FV_ENABLED
    ,NVisu_FV,UVisu_FV&
#endif /* FV_ENABLED */
    )

#if USE_MPI
   IF ((.NOT.MPIRoot).AND.(Avg2d)) THEN
     ! For parallel averaging, all data is gathered on the root. Disable output for other procs.
     nElemsAvg2D_DG = 0
     nElemsAvg2D_FV = 0
     SDEALLOCATE(UVisu_DG)
     SDEALLOCATE(UVisu_FV)
     ALLOCATE(UVisu_DG(0:NVisu   ,0:NVisu   ,0:0,nElemsAvg2D_DG,nVarVisu))
     ALLOCATE(UVisu_FV(0:NVisu_FV,0:NVisu_FV,0:0,nElemsAvg2D_FV,nVarVisu))
   END IF
#endif

  ! Convert coordinates to visu grid
  IF (changedMeshFile.OR.changedNVisu.OR.changedFV_Elems.OR.changedDGonly.OR.changedAvg2D) &
    CALL BuildVisuCoords()

  IF (doSurfVisu .AND. &
    ! Convert surface coordinates to visu grid
    (changedMeshFile.OR.changedNVisu.OR.changedFV_Elems.OR.changedDGonly.OR.changedBCnames)) &
    CALL BuildSurfVisuCoords()
END IF

MeshFile_old          = MeshFile
prmfile_old           = prmfile
statefile_old         = statefile
NVisu_old             = NVisu
NCalc_old             = NCalc
nVar_State_old        = nVar_State
withDGOperator_old    = withDGOperator
DGonly_old            = DGonly
Avg2D_old             = Avg2D
NodeTypeVisuPosti_old = NodeTypeVisuPosti
NState_old            = PP_N
RestartMode           = -1

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,A)') " Visu finished for state file: ", TRIM(statefile)
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE visu

!===================================================================================================================================
!> Deallocate arrays used by visu.
!===================================================================================================================================
SUBROUTINE FinalizeVisu()
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments,ONLY: FinalizeCommandlineArguments
USE MOD_DG                   ,ONLY: FinalizeDG
USE MOD_DG_Vars
USE MOD_Equation             ,ONLY: FinalizeEquation
USE MOD_Filter               ,ONLY: FinalizeFilter
USE MOD_Interpolation        ,ONLY: FinalizeInterpolation
USE MOD_IO_HDF5              ,ONLY: FinalizeIOHDF5
USE MOD_Mesh                 ,ONLY: FinalizeMesh
USE MOD_Mesh_Vars            ,ONLY: Elem_xGP
USE MOD_Mortar               ,ONLY: FinalizeMortar
USE MOD_Overintegration      ,ONLY: FinalizeOverintegration
USE MOD_ReadInTools          ,ONLY: FinalizeParameters
USE MOD_Restart              ,ONLY: FinalizeRestart
USE MOD_Visu_Vars
#if PARABOLIC
USE MOD_Lifting              ,ONLY: FinalizeLifting
#endif
#if FV_ENABLED
USE MOD_Indicator            ,ONLY: FinalizeIndicator
USE MOD_FV_Basis             ,ONLY: FinalizeFV_Basis
#endif /* FV_ENABLED */
#if USE_MPI
USE MOD_MPI                  ,ONLY: FinalizeMPI
#endif /* USE_MPI */
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL :: Time,SimulationTime,mins,secs,hours,days
!===================================================================================================================================

IF(MPIRoot)THEN
  IF(FILEEXISTS('.posti.ini'))THEN
    OPEN(UNIT=31, FILE='.posti.ini', STATUS='old')
    CLOSE(31, STATUS='delete')
  END IF
  IF(FILEEXISTS('.flexi.ini'))THEN
    OPEN(UNIT=31, FILE='.flexi.ini', STATUS='old')
    CLOSE(31, STATUS='delete')
  END IF
END IF

! Reset all strings and variables
prmfile_old       = ''
statefile_old     = ''
MeshFile          = ''
MeshFile_old      = ''
NodeTypeVisuPosti = 'VISU'
NodeTypeVisuPosti_old = ''
NVisu     = -1
NVisu_old = -1
nVar_State_old = -1
NState_old = -1
withDGOperator_old = .FALSE.
hasFV_Elems = .FALSE.

SDEALLOCATE(mapDepToCalc)
#if FV_ENABLED && FV_RECONSTRUCT
SDEALLOCATE(mapDepToCalc_FV)
#endif
SDEALLOCATE(mapAllVarsToVisuVars)
SDEALLOCATE(mapAllVarsToSurfVisuVars)
SDEALLOCATE(mapAllVarsToSurfVisuVars_old)
SDEALLOCATE(UCalc_DG)
SDEALLOCATE(UCalc_FV)
SDEALLOCATE(FV_Elems_loc)

SDEALLOCATE(mapDGElemsToAllElems)
SDEALLOCATE(mapFVElemsToAllElems)

SDEALLOCATE(CoordsVisu_DG)
SDEALLOCATE(UVisu_DG)
SDEALLOCATE(CoordsVisu_FV)
SDEALLOCATE(UVisu_FV)
SDEALLOCATE(U)
SDEALLOCATE(Elem_xGP)

CALL FinalizeRestart()
CALL FinalizeEquation()
CALL FinalizeDG()
CALL FinalizeOverintegration()
CALL FinalizeFilter()
CALL FinalizeCommandlineArguments()
CALL FinalizeParameters()
CALL FinalizeInterpolation()
CALL FinalizeMesh()
CALL FinalizeIOHDF5()
#if FV_ENABLED
CALL FinalizeIndicator()
CALL FinalizeFV_Basis()
#endif /* FV_ENABLED */
CALL FinalizeMortar()
#if PARABOLIC
CALL FinalizeLifting()
#endif /*PARABOLIC*/

SDEALLOCATE(VarNamesHDF5)
SDEALLOCATE(VarnamesAll)
SDEALLOCATE(BCNamesAll)
SDEALLOCATE(DepTable)
SDEALLOCATE(DepSurfaceOnly)
SDEALLOCATE(DepVolumeOnly)

SDEALLOCATE(FV_Elems_old)
SDEALLOCATE(mapDepToCalc_FV)
SDEALLOCATE(mapAllBCSidesToDGVisuBCSides)
SDEALLOCATE(mapAllBCSidesToFVVisuBCSides)
SDEALLOCATE(mapAllBCNamesToVisuBCNames_old)
SDEALLOCATE(mapAllBCNamesToVisuBCNames)
SDEALLOCATE(nSidesPerBCNameVisu_DG)
SDEALLOCATE(nSidesPerBCNameVisu_FV)

! Calculate simulation time
Time = FLEXITIME()
SimulationTime = Time-StartTime

! Get secs, mins, hours and days
secs = MOD(SimulationTime,60.)
SimulationTime = SimulationTime / 60.
mins = MOD(SimulationTime,60.)
SimulationTime = SimulationTime / 60.
hours = MOD(SimulationTime,24.)
SimulationTime = SimulationTime / 24.
!days = MOD(SimulationTime,365.) ! Use this if years are also to be displayed
days = SimulationTime

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F16.2,A)',ADVANCE='NO')  ' VISU  FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(A2,I6,A1,I0.2,A1,I0.2,A1,I0.2,A1)') ' [',INT(days),':',INT(hours),':',INT(mins),':',INT(secs),']'
SWRITE(UNIT_stdOut,'(132("="))')

#if USE_MPI
CALL FinalizeMPI()
#endif /* USE_MPI */

END SUBROUTINE FinalizeVisu

END MODULE MOD_Visu
