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

!==================================================================================================================================
!> Contains control routines to read high-order meshes, provide mesh data to the solver, build the metrics, partition the domain.
!==================================================================================================================================
MODULE MOD_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::FinalizeMesh
!==================================================================================================================================

PUBLIC::DefineParametersMesh
CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateStringOption(  'MeshFile',            "(relative) path to meshfile (mandatory).")
CALL prms%CreateLogicalOption( 'useCurveds',          "Controls usage of high-order information in mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption( 'interpolateFromTree', "For non-conforming meshes, built by refinement from a tree structure, "//&
                                                      "the metrics can be built from the tree geometry if it is contained "//&
                                                      "in the mesh. Can improve free-stream preservation.",&
                                                      '.TRUE.')
CALL prms%CreateRealOption(    'meshScale',           "Scale the mesh by this factor (shrink/enlarge).",&
                                                      '1.0')
CALL prms%CreateLogicalOption( 'meshdeform',          "Apply simple sine-shaped deformation on cartesion mesh (for testing).",&
                                                      '.FALSE.')
#if (PP_dim == 3)
CALL prms%CreateLogicalOption( 'crossProductMetrics', "Compute mesh metrics using cross product form. Caution: in this case "//&
                                                      "free-stream preservation is only guaranteed for N=3*NGeo.",&
                                                      '.FALSE.')
#endif
CALL prms%CreateIntOption(     'debugmesh',           "Output file with visualization and debug information for the mesh. "//&
                                                      "0: no visualization, 3: Paraview binary",'0')
CALL prms%CreateStringOption(  'BoundaryName',        "Names of boundary conditions to be set (must be present in the mesh!)."//&
                                                      "For each BoundaryName a BoundaryType needs to be specified.",&
                                                      multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType',        "Type of boundary conditions to be set. Format: (BC_TYPE,BC_STATE)",&
                                                      multiple=.TRUE.)
CALL prms%CreateLogicalOption( 'writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
CALL prms%CreateIntOption(     'NGeoOverride',        "Override switch for NGeo. Interpolate mesh to different NGeo." //&
                                                      "<1: off, >0: Interpolate",'-1')
END SUBROUTINE DefineParametersMesh


!==================================================================================================================================
!> Routine controlling the initialization of the mesh.
!> - parameter and mesh reading
!> - domain partitioning
!> - allocation of mesh arrays
!> - build mesh mappings to handle volume/surface operations
!> - compute the mesh metrics
!> - provide mesh metrics for overintegration
!==================================================================================================================================
SUBROUTINE InitMesh(meshMode,MeshFile_IN)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_HDF5_Input
USE MOD_ChangeBasisByDim   ,ONLY:ChangeBasisVolume
USE MOD_Interpolation      ,ONLY:GetVandermonde
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone,NodeType,NodeTypeVISU
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR,GETREAL,GETINT
USE MOD_Metrics,            ONLY:BuildCoords,CalcMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Mappings,           ONLY:buildMappings
#if USE_MPI
USE MOD_Prepare_Mesh,       ONLY:exchangeFlip
#endif
#if FV_ENABLED
USE MOD_FV_Metrics,         ONLY:InitFV_Metrics
#endif
USE MOD_IO_HDF5,            ONLY:AddToElemData,ElementOut
#if (PP_dim == 2)
USE MOD_2D
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode !< 0: only read and build Elem_xGP,
                               !< 1: as 0 + build connectivity, 2: as 1 + calc metrics
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: MeshFile_IN !< file name of mesh to be read
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: x(3),meshScale
REAL,POINTER      :: coords(:,:,:,:,:),coordsTmp(:,:,:,:,:),Vdm_EQNGeo_EQNGeoOverride(:,:)
INTEGER           :: iElem,i,j,k,nElemsLoc
LOGICAL           :: validMesh
INTEGER           :: firstMasterSide     ! lower side ID of array U_master/gradUx_master...
INTEGER           :: lastMasterSide      ! upper side ID of array U_master/gradUx_master...
INTEGER           :: firstSlaveSide      ! lower side ID of array U_slave/gradUx_slave...
INTEGER           :: lastSlaveSide       ! upper side ID of array U_slave/gradUx_slave...
INTEGER           :: iSide,LocSideID,SideID
INTEGER           :: NGeoOverride
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'InitMesh not ready to be called or already called.')
END IF
IF(isMortarMesh) THEN
  CALL CollectiveStop(__STAMP__,&
    'Reduced Flexi does not support Mortars.')
END IF

SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,I1,A)') ' INIT MESH IN MODE ',meshMode,'...'

! prepare pointer structure (get nElems, etc.)
IF (PRESENT(MeshFile_IN)) THEN
  MeshFile = MeshFile_IN
ELSE
  MeshFile = GETSTR('MeshFile')
END IF
validMesh = ISVALIDMESHFILE(MeshFile)
IF(.NOT.validMesh) &
    CALL CollectiveStop(__STAMP__,'ERROR - Mesh file not a valid HDF5 mesh.')

useCurveds=GETLOGICAL('useCurveds')
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo)
CALL CloseDataFile()

IF(useCurveds.AND.(PP_N.LT.NGeo))THEN
  SWRITE(UNIT_stdOut,'(A)') 'WARNING: N<NGeo, for curved hexa normals are only approximated,&
                           & can cause problems on periodic boundaries! Set N>=NGeo'
ENDIF

CALL readMesh(MeshFile) !set nElems

#if (PP_dim == 2)
! If this is a two dimensional calculation, all subsequent operations are performed on the reduced mesh.
SWRITE(UNIT_stdOut,'(A)') " RUNNING A 2D SIMULATION! "
! The mesh coordinates read in by the readMesh routine are therefore reduced by one dimension.
CALL to2D_rank5((/1,0,0,0,1/),(/3,NGeo,NGeo,NGeo,nElems/),4,NodeCoords)
NodeCoords(3,:,:,:,:) = 0.
#endif

! if trees are available: compute metrics on tree level and interpolate to elements
interpolateFromTree=.FALSE.
IF(isMortarMesh) interpolateFromTree=GETLOGICAL('interpolateFromTree')
IF(interpolateFromTree)THEN
#if (PP_dim == 2)
  CALL CollectiveStop(__STAMP__,&
      "interpolateFromTree not supported in 2D.")
#endif
  coords=>TreeCoords
  NGeo=NGeoTree
  nElemsLoc=nTrees
ELSE
  coords=>NodeCoords
  nElemsLoc=nElems
ENDIF

NGeoOverride=GETINT('NGeoOverride')
IF(NGeoOverride.GT.0)THEN
  ALLOCATE(CoordsTmp(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
  ALLOCATE(Vdm_EQNGeo_EQNGeoOverride(0:NGeoOverride,0:NGeo))
  CALL GetVandermonde(Ngeo, NodeTypeVISU, NgeoOverride, NodeTypeVISU, Vdm_EQNgeo_EQNgeoOverride, modal=.FALSE.)
  DO iElem=1,nElemsLoc
    CALL ChangeBasisVolume(3,Ngeo,NgeoOverride,Vdm_EQNGeo_EQNGeoOverride,coords(:,:,:,:,iElem),coordsTmp(:,:,:,:,iElem))
  END DO
  ! cleanup
  IF(interpolateFromTree)THEN
    DEALLOCATE(TreeCoords)
    ALLOCATE(TreeCoords(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
    coords=>TreeCoords
  ELSE
    DEALLOCATE(NodeCoords)
    ALLOCATE(NodeCoords(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
    coords=>NodeCoords
  END IF
  Coords = CoordsTmp
  NGeo = NGeoOverride
  DEALLOCATE(CoordsTmp, Vdm_EQNGeo_EQNGeoOverride)
END IF

SWRITE(UNIT_stdOut,'(a3,a30,a3,i0)')' | ','Ngeo',' | ', Ngeo

! scale and deform mesh if desired (warning: no mesh output!)
meshScale=GETREAL('meshScale')
IF(ABS(meshScale-1.).GT.1e-14)&
  Coords = Coords*meshScale

IF(GETLOGICAL('meshdeform'))THEN
  DO iElem=1,nElemsLoc
    DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
      x(:)=Coords(:,i,j,k,iElem)
#if PP_dim==3
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))*SIN(PP_Pi*x(3))
#else
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))
#endif
    END DO; END DO; END DO;
  END DO
END IF

! Build the coordinates of the solution gauss points in the volume
ALLOCATE(Elem_xGP(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
IF(interpolateFromTree)THEN
  CALL BuildCoords(NodeCoords,NodeType,PP_N,Elem_xGP,TreeCoords)
ELSE
  CALL BuildCoords(NodeCoords,NodeType,PP_N,Elem_xGP)
ENDIF
!@cuf ALLOCATE(d_Elem_xGP(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
!@cuf d_Elem_xGP=Elem_xGP

! Return if no connectivity and metrics are required (e.g. for visualization mode)
IF (meshMode.GT.0) THEN

  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
  CALL setLocalSideIDs()

#if USE_MPI
  ! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
  CALL exchangeFlip()
#endif

  !RANGES
  !-----------------|-----------------|-------------------|
  !    U_master     | U_slave         |    FLUX           |
  !-----------------|-----------------|-------------------|
  !  BCsides        |                 |    BCSides        |
  !  InnerMortars   |                 |    InnerMortars   |
  !  InnerSides     | InnerSides      |    InnerSides     |
  !  MPI_MINE sides | MPI_MINE sides  |    MPI_MINE sides |
  !                 | MPI_YOUR sides  |    MPI_YOUR sides |
  !  MPIMortars     |                 |    MPIMortars     |
  !-----------------|-----------------|-------------------|

  firstBCSide          = 1
  firstMortarInnerSide = firstBCSide         +nBCSides
  firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
  firstMPISide_MINE    = firstInnerSide      +nInnerSides
  firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
  firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR

  lastBCSide           = firstMortarInnerSide-1
  lastMortarInnerSide  = firstInnerSide    -1
  lastInnerSide        = firstMPISide_MINE -1
  lastMPISide_MINE     = firstMPISide_YOUR -1
  lastMPISide_YOUR     = firstMortarMPISide-1
  lastMortarMPISide    = nSides


  firstMasterSide = 1
  lastMasterSide  = nSides
  firstSlaveSide  = firstInnerSide
  lastSlaveSide   = lastMPISide_YOUR
  nSidesMaster    = lastMasterSide-firstMasterSide+1
  nSidesSlave     = lastSlaveSide -firstSlaveSide+1

  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8)')   'first/lastMasterSide     ', firstMasterSide,lastMasterSide
  LOGWRITE(*,'(A25,I8)')   'first/lastSlaveSide      ', firstSlaveSide, lastSlaveSide
  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8,I8)')'first/lastBCSide         ', firstBCSide         ,lastBCSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMortarInnerSide', firstMortarInnerSide,lastMortarInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastInnerSide      ', firstInnerSide      ,lastInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_MINE   ', firstMPISide_MINE   ,lastMPISide_MINE
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_YOUR   ', firstMPISide_YOUR   ,lastMPISide_YOUR
  LOGWRITE(*,'(A30,I8,I8)')'first/lastMortarMPISide  ', firstMortarMPISide  ,lastMortarMPISide
  LOGWRITE(*,*)'-------------------------------------------------------'

  ! fill ElemToSide, SideToElem,BC
  ALLOCATE(ElemToSide(3,6,nElems))
!@cuf ALLOCATE(d_ElemToSide(3,6,nElems))
  ALLOCATE(SideToElem(5,nSides))
!@cuf ALLOCATE(d_SideToElem(5,nSides))

  ALLOCATE(BC(1:nBCSides))
  ALLOCATE(AnalyzeSide(1:nSides))
  ElemToSide  = 0
  SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
  BC          = 0
  AnalyzeSide = 0

  !NOTE: nMortarSides=nMortarInnerSides+nMortarMPISides
  ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo
  ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides)) ! [1]: 1: Neighbour sides, 2: Flip, [2]: small sides
  MortarType=0
  MortarInfo=-1

  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
  CALL fillMeshInfo()

#if (PP_dim ==2)
  ! In 2D, there is only one flip for the slave sides (1)
  SideToElem(S2E_FLIP,:) = MIN(1,SideToElem(S2E_FLIP,:))
  ElemToSide(:,1,:) = -999
  ElemToSide(:,6,:) = -999
  ElemToSide(E2S_FLIP,:,:) = MIN(1,ElemToSide(E2S_FLIP,:,:))
  MortarInfo(MI_FLIP,:,:)  = MIN(1,MortarInfo(MI_FLIP,:,:))
#endif

  ! Build necessary mappings
  CALL buildMappings(PP_N,V2S=V2S,S2V=S2V,S2V2=S2V2,S2V2_inv=S2V2_inv,FS2M=FS2M,dim=PP_dim)

  ! Copy mappings to GPU
#if PP_dim==2
!@cuf ALLOCATE(d_S2V2(    2,0:PP_N,0:PP_NZ,0:1,2:5))
!@cuf ALLOCATE(d_S2V2_inv(2,0:PP_N,0:PP_NZ,0:1,2:5))
#else
!@cuf ALLOCATE(d_S2V2(    2,0:PP_N,0:PP_NZ,0:4,1:6))
!@cuf ALLOCATE(d_S2V2_inv(2,0:PP_N,0:PP_NZ,0:4,1:6))
#endif
!@cuf d_S2V2     = S2V2  
!@cuf d_S2V2_inv = S2V2_inv

END IF
! deallocate pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

IF (meshMode.GT.1) THEN

  ! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----

  ! volume data
  ALLOCATE(      dXCL_N(3,3,0:PP_N,0:PP_N,0:PP_NZ,nElems)) ! temp
  ALLOCATE(Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  ALLOCATE(Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  ALLOCATE(Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  !@cuf ALLOCATE(d_Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  !@cuf ALLOCATE(d_Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  !@cuf ALLOCATE(d_Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  ALLOCATE(            sJ(  0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  !@cuf ALLOCATE(    d_sJ(  0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
  ALLOCATE(     scaledJac(  0:PP_N,0:PP_N,0:PP_NZ,nElems))
  NGeoRef=3*NGeo ! build jacobian at higher degree
  ALLOCATE(    DetJac_Ref(1,0:NGeoRef,0:NGeoRef,0:ZDIM(NGeoRef),nElems))

  ! surface data
  ALLOCATE(      Face_xGP(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
  ALLOCATE(       NormVec(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
!@cuf ALLOCATE( d_NormVec(3,0:PP_N,0:PP_NZ,1:nSides))
  ALLOCATE(      TangVec1(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
!@cuf ALLOCATE(d_TangVec1(3,0:PP_N,0:PP_NZ,1:nSides))
  ALLOCATE(      TangVec2(3,0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
!@cuf ALLOCATE(d_TangVec2(3,0:PP_N,0:PP_NZ,1:nSides))
  ALLOCATE(      SurfElem(  0:PP_N,0:PP_NZ,0:FV_SIZE,1:nSides))
!@cuf ALLOCATE(d_SurfElem(  0:PP_N,0:PP_NZ,1:nSides))
  ALLOCATE(     Ja_Face(3,3,0:PP_N,0:PP_NZ,          1:nSides)) ! temp

#if FV_ENABLED
  ! NOTE: initialize with 1 and not with 0
  ALLOCATE(sJ_master(1,0:PP_N,0:PP_NZ,1:nSides,0:FV_SIZE))
  ALLOCATE(sJ_slave( 1,0:PP_N,0:PP_NZ,1:nSides,0:FV_SIZE))
  sJ_master = 1.
  sJ_slave  = 1.
#endif

  ! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
  ! assign 1/detJ (sJ)
  ! assign normal and tangential vectors and surfElems on faces

#if (PP_dim == 3)
  ! compute metrics using cross product instead of curl form (warning: no free stream preservation!)
  crossProductMetrics=GETLOGICAL('crossProductMetrics')
#endif
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."
  CALL CalcMetrics()     ! DG metrics
  !@cuf d_sJ             = sJ
  !@cuf d_Metrics_fTilde = Metrics_fTilde
  !@cuf d_Metrics_gTilde = Metrics_gTilde
  !@cuf d_Metrics_hTilde = Metrics_hTilde
#if FV_ENABLED
  CALL InitFV_Metrics()  ! FV metrics
#endif
  ! debugmesh: param specifies format to output, 0: no output, 1: tecplot ascii, 2: tecplot binary, 3: paraview binary
  CALL WriteDebugMesh(GETINT('debugmesh'))
END IF

IF (meshMode.GT.0) THEN
  ALLOCATE(SideToGlobalSide(nSides))
  DO iElem=1,nElems
#if PP_dim == 3
    DO LocSideID=1,6
#else
    DO LocSideID=2,5
#endif
      SideID = ElemToSide(E2S_SIDE_ID,LocSideID,iElem)
      iSide = ElemInfo(3,iElem+offsetElem) + LocSideID
      SideToGlobalSide(SideID) = ABS(SideInfo(2,iSide))
      ! Fill Master/Slave Info
      IF (SideToElem(S2E_ELEM_ID,SideID).EQ.iElem) THEN
        ElemToSide(E2S_IS_MASTER,LocSideID,iElem) = 1
      ELSE IF (SideToElem(S2E_NB_ELEM_ID,SideID).EQ.iElem) THEN
        ElemToSide(E2S_IS_MASTER,LocSideID,iElem) = 0
      ELSE
        CALL ABORT(__STAMP__,"Seems like an error in side connectivity!")
      END IF
    END DO
  END DO ! iElem
END IF

SDEALLOCATE(dXCL_N)
SDEALLOCATE(Ja_Face)
SDEALLOCATE(TreeCoords)
SDEALLOCATE(xiMinMax)
SDEALLOCATE(ElemToTree)
IF ((.NOT.postiMode).AND.(ALLOCATED(scaledJac))) DEALLOCATE(scaledJac)
CALL AddToElemData(ElementOut,'myRank',IntScalar=myRank)

! Copy more arrays to GPU
!@cuf d_NormVec    = NormVec( :,:,:,0,:)
!@cuf d_TangVec1   = TangVec1(:,:,:,0,:)
!@cuf d_TangVec2   = TangVec2(:,:,:,0,:)
!@cuf d_SurfElem   = SurfElem(  :,:,0,:)
!@cuf d_ElemToSide = ElemToSide
!@cuf d_SideToElem = SideToElem

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitMesh

!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE FinalizeMesh()
! MODULES
USE MOD_Mesh_Vars
USE MOD_Mappings  ,ONLY:FinalizeMappings
#if FV_ENABLED
USE MOD_FV_Vars   ,ONLY:FV_Elems_master
USE MOD_FV_Metrics,ONLY:FinalizeFV_Metrics
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
SDEALLOCATE(AnalyzeSide)

! mortars
SDEALLOCATE(MortarType)
SDEALLOCATE(MortarInfo)

! allocated during ReadMesh
SDEALLOCATE(NodeCoords)
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)

! Volume
SDEALLOCATE(Elem_xGP)
!@cuf SDEALLOCATE(d_Elem_xGP)
SDEALLOCATE(Metrics_fTilde)
SDEALLOCATE(Metrics_gTilde)
SDEALLOCATE(Metrics_hTilde)
SDEALLOCATE(sJ)
SDEALLOCATE(scaledJac)
SDEALLOCATE(DetJac_Ref)
#if FV_ENABLED
SDEALLOCATE(sJ_master)
SDEALLOCATE(sJ_slave)
#endif

! surface
SDEALLOCATE(Face_xGP)
SDEALLOCATE(NormVec)
SDEALLOCATE(TangVec1)
SDEALLOCATE(TangVec2)
SDEALLOCATE(SurfElem)

! ijk sorted mesh
SDEALLOCATE(Elem_IJK)
SDEALLOCATE(ElemInfo)
SDEALLOCATE(SideInfo)
SDEALLOCATE(SideToGlobalSide)

! mappings
CALL FinalizeMappings()

#if FV_ENABLED
SDEALLOCATE(FV_Elems_master) ! moved here from fv.f90
CALL FinalizeFV_Metrics()
#endif

MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
