!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
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
!> Main program for the HIT_Filter tool. It reads in statefiles and applies a Fourier cutoff filter to the state by applying FFTs
!> and writes the filtered state to a new statefile.
!===================================================================================================================================
PROGRAM HIT_Filter
! MODULES
USE MOD_Globals
USE MOD_HIT_Filter
USE MOD_HIT_Filter_Vars
USE MOD_ReadInTools
USE MOD_Commandline_Arguments
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_HIT_FFT,                 ONLY: InitFFT,FinalizeFFT
USE MOD_HIT_FFT_Vars,            ONLY: N_FFT,NCalc
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_Mesh_Vars,               ONLY: nElems_IJK,Elem_IJK,MeshFile
USE MOD_Mesh_ReadIn,             ONLY: ReadIJKSorting
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5,FinalizeIOHDF5
USE MOD_IO_HDF5,                 ONLY: FieldOut,AddToFieldData,FinalizeFieldData
USE MOD_HDF5_Output,             ONLY: WriteState
USE MOD_HDF5_Input,              ONLY: ISVALIDMESHFILE
USE MOD_Output,                  ONLY: DefineParametersOutput
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
#if USE_OPENMP
USE OMP_Lib
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: Time
INTEGER                            :: iArg
INTEGER                            :: N_HDF5_old=0            ! Polynomial degree N of last statefile
LOGICAL                            :: changedMeshFile=.FALSE. ! True if mesh between states changed
LOGICAL                            :: changedN       =.FALSE. ! True if N between states changes
CHARACTER(LEN=255)                 :: InputStateFile          ! StateFile to be processed
CHARACTER(LEN=255)                 :: MeshFile_old = ''       ! MeshFile of last statefile
CHARACTER(LEN=255)                 :: MeshFile_prm = ''       ! Meshfile of input parameter file
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

! Define Parameters
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersOutput()
CALL DefineParametersMesh()

! Parameters for HIT_Filter
CALL prms%SetSection("HIT_Filter")
CALL prms%CreateIntOption("N_Filter" , "Maximum wavenumber for cutoff filter.")
CALL prms%CreateIntOption("NCalc"    , "Polynomial degree in each element for global DFFT basis.")

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file and min. 1 statefile is given
IF ((nArgs.LT.2).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: posti_hit_filter parameter.ini [statefile.h5, ....]')
END IF

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(11X,A)') &
    "  ____  ____ _____ _________                      _________ _____ _____    _________ _________ _______     "
SWRITE(UNIT_stdOut,'(11X,A)') &
    " |_   ||   _|_   _|  _   _  |                    |_   ___  |_   _|_   _|  |  _   _  |_   ___  |_   __ \    "
SWRITE(UNIT_stdOut,'(11X,A)') &
    "   | |__| |   | | |_/ | | \_|       ______         | |_  \_| | |   | |    |_/ | | \_| | |_  \_| | |__) |   "
SWRITE(UNIT_stdOut,'(11X,A)') &
    "   |  __  |   | |     | |          |______|        |  _|     | |   | |   _    | |     |  _|  _  |  __ /    "
SWRITE(UNIT_stdOut,'(11X,A)') &
    "  _| |  | |_ _| |_   _| |_                        _| |_     _| |_ _| |__/ |  _| |_   _| |___/ |_| |  \ \_  "
SWRITE(UNIT_stdOut,'(11X,A)') &
    " |____||____|_____| |_____|                      |_____|   |_____|________| |_____| |_________|____| |___| "
SWRITE(UNIT_stdOut,'(A)')
SWRITE(UNIT_stdOut,'(132("="))')

! Parse parameter file
CALL prms%read_options(Args(1))
ParameterFile = Args(1)

! Readin Parameters
N_Filter = GETINT('N_Filter','-1')
NCalc    = GETINT('NCalc')
MeshFile_prm = GETSTR('MeshFile','')

! Initialize IO
CALL InitIOHDF5()

! Loop over all files specified on commandline
DO iArg=2,nArgs
  Time = OMP_FLEXITIME()

  InputStateFile = Args(iArg)

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,I5,A,I5,A)') ' PROCESSING FILE ',iArg-1,' of ',nArgs-1,' FILES.'
  SWRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputStateFile),'" )'
  SWRITE(UNIT_stdOut,'(132("="))')

  ! Read attributes and solution from state file
  CALL ReadOldStateFile(InputStateFile)

  ! Check if input attributes have changed since last state file
  IF(TRIM(MeshFile).NE.TRIM(MeshFile_old)) changedMeshFile =.TRUE.
  IF(N_HDF5.NE.N_HDF5_old)                 changedN        =.TRUE.

  ! Re-initialize interpolation and re-allocate DG solution array if N has changed
  IF(changedN) THEN
    CALL FinalizeInterpolation()
    CALL InitInterpolation(N_HDF5)
  END IF

  ! Re-initialize mesh if it has changed
  IF(changedMeshFile) THEN
    CALL FinalizeMesh()
    CALL DefineParametersMesh()
    ! Only take the mesh file deposited in the state file if it is a valid mesh file.
    ! Otherwise try the mesh from the input parameter file
    IF(FILEEXISTS(MeshFile).AND.ISVALIDMESHFILE(MeshFile)) THEN
      SWRITE(UNIT_stdOut,*) "INITIALIZING MESH FROM FILE """,TRIM(MeshFile),""""
      CALL InitMesh(MeshMode=0,MeshFile_In=MeshFile)
    ELSE
      SWRITE(UNIT_stdOut,*) "WARNING: No valid mesh file is given in HDF5 attributes of current state file! &
                                    & Reading mesh from parameter file instead..."
      CALL InitMesh(MeshMode=0,MeshFile_In=MeshFile_prm)
    END IF
    CALL ReadIJKSorting() ! Read global xyz sorting of structured mesh

    ! Check if found
    IF(.NOT.ALLOCATED(Elem_IJK)) THEN
      CALL Abort(__STAMP__,'Mesh does not contain IJK sorting which is however mandatory!')
    END IF

    ! Currently only cubic meshes are allowed!
    IF(.NOT.((nElems_IJK(1).EQ.nElems_IJK(2)).AND.(nElems_IJK(1).EQ.nElems_IJK(3)))) THEN
      CALL Abort(__STAMP__,'Mesh does not have the same amount of elements in x,y and z!')
    END IF

    ! Get new number of points for fourier analysis
    N_FFT=(NCalc+1)*nElems_IJK(1)
  END IF

  IF(changedMeshFile .OR. changedN) THEN
    CALL FinalizeFFT()
    CALL InitFFT()
  END IF

  ! Transform global solution into Fourier space and apply filter there.
  CALL FourierFilter(nVar_HDF5,U)
  IF(FieldDataExists) THEN
    CALL FourierFilter(nVarField_HDF5,FieldData)
    CALL AddToFieldData(FieldOut,(/nVarField_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5+1/), &
                               'FieldData',VarNames_FieldData,RealArray=FieldData)
  END IF

  ! Write State-File
  CALL WriteNewStateFile()

  ! To determine whether meshfile or N changes
  MeshFile_old    = MeshFile
  N_HDF5_old      = N_HDF5
  changedMeshFile = .FALSE.
  changedN        = .FALSE.

  ! Deallocate DG solution array and FieldData for next file
  DEALLOCATE(U)
  IF(FieldDataExists) THEN
    DEALLOCATE(FieldData)
    DEALLOCATE(VarNames_FieldData)
    CALL FinalizeFieldData(FieldOut)
    FieldDataExists = .FALSE.
  END IF

  SWRITE(UNIT_stdOut,'(132("="))')
  SWRITE(UNIT_stdOut,'(A,A,A,F0.3,A)') ' PROCESSED FILE ',TRIM(InputStateFile),' in [',OMP_FLEXITIME()-Time,'s]'
END DO !iArg=1,nArgs

! Finalize
CALL FinalizeParameters()
CALL FinalizeInterpolation()
CALL FinalizeMesh()
CALL FinalizeFFT()
CALL FinalizeIOHDF5
#if USE_MPI
CALL FinalizeMPI()
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') ' HIT_FILTER FINISHED! '
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM HIT_Filter
