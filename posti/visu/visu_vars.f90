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
!==================================================================================================================================
!> Contains global variables provided by the visu routines
!==================================================================================================================================
MODULE MOD_Visu_Vars
USE ISO_C_BINDING
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
CHARACTER(LEN=255)                :: fileType = ''               !< possible values:
                                                                 !< * 'State' for FLEXI state files matching the compiled EOS
                                                                 !< * 'Generic'
                                                                 !< * 'Mesh'
CHARACTER(LEN=255)                :: prmfile_old = ''            !< saves the filename of the previous FLEXI parameter file
CHARACTER(LEN=255)                :: statefile_old = ''          !< saves the filename of the previous state (*.h5)
CHARACTER(LEN=255)                :: MeshFile = ''               !< actual filename of the mesh used for visualization
CHARACTER(LEN=255)                :: MeshFile_state = ''         !< filename of the mesh given in the state file
CHARACTER(LEN=255)                :: MeshFile_old = ''           !< saves  previous MeshFile
CHARACTER(LEN=255)                :: NodeTypeVisuPosti = "VISU"  !< NodeType used for visualization output
CHARACTER(LEN=255)                :: NodeTypeVisuPosti_old = ''  !< saves previous NodeType
INTEGER                           :: NVisu                       !< polynomial degree of the visualization
INTEGER                           :: NVisu_old = -1              !< saves previous NVisu
INTEGER                           :: NVisu_FV                    !< number of output points for FV elements (always == 2*(PP_N+1))
INTEGER                           :: NCalc_FV                    !< number of calculation points for FV elements (NVisu_FV or PP_N)
INTEGER                           :: NCalc                       !< Different polynomial degree to do calculations on
INTEGER                           :: NCalc_old                   !< Different polynomial degree to do calculations on
INTEGER                           :: nVarIni                     !< number of requested variables in parameter file
INTEGER                           :: nVar_State                  !< number of variables in the state file
INTEGER                           :: nVar_State_old = -1         !< saves previous nVar_State
INTEGER                           :: nState_old = -1             !< saves previous PP_N
INTEGER                           :: nElems_DG                   !< number of DG elements in state
INTEGER                           :: nElems_FV                   !< number of FV elements in state
LOGICAL                           :: withDGOperator              !< flag indicating if call of 'DGTimeDerivative' is required
LOGICAL                           :: withDGOperator_old = .FALSE.!< saves previous withDGOperator
REAL                              :: OutputTime                  !< simulation time of actual state file
LOGICAL                           :: hasFV_Elems = .FALSE.       !< flag indicating if state file contains any FV elements
LOGICAL                           :: DGonly = .FALSE.            !< flag to force visualization of FV elements as DG elements
LOGICAL                           :: DGonly_old = .TRUE.         !< saves previous DGonly
LOGICAL                           :: HighOrder = .FALSE.         !< flag to enable high-order element representation
INTEGER,ALLOCATABLE               :: mapDGElemsToAllElems(:)     !< maps element index of DG elements to all elements
INTEGER,ALLOCATABLE               :: mapFVElemsToAllElems(:)     !< maps element index of FV elements to all elements
INTEGER,ALLOCATABLE               :: FV_Elems_loc(:)             !< current distribution of FV/DG elems
INTEGER,ALLOCATABLE               :: FV_Elems_old(:)             !< previous distribution of FV/DG elems
INTEGER                           :: meshMode_old=0              !< Used to check if InitMesh must be called again with different
                                                                 !< mesh mode
LOGICAL                           :: StateFileMode               !< Flag indicating if a state file is being visualized. Only then
                                                                 !< calculations can be performed, otherwise only pure visu.
LOGICAL                           :: MeshFileMode                !< Flag indicating a mesh file should be visualized
LOGICAL                           :: doSurfVisu                  !< Flag indicating if any surfaces need to be visualized
LOGICAL                           :: Avg2D                       !< Flag indicating if solution should be averaged in zeta dir
LOGICAL                           :: Avg2D_old = .FALSE.         !< Previus state of Avg2D flag, used to check for change
LOGICAL                           :: Avg2DHDF5Output             !< Flag indicating if the averaged solution should be written to a
                                                                 !< .h5 file
LOGICAL                           :: DependenciesOutputDone=.FALSE. !< Flag indicating if dependency table was already written

! The following flags indicate if during successive visualizations of (different) state files the respective properties
! changed. For example the mesh file of different state files in a timeseries is the same ...
LOGICAL                           :: changedStateFile            !< .h5 state file to visualize changed
LOGICAL                           :: changedMeshFile             !< Mesh file changed
LOGICAL                           :: changedNVisu                !< Polyomial degree for visualization changed
LOGICAL                           :: changedNCalc                !< Polyomial degree for calculation changed
LOGICAL                           :: changedVarNames             !< variables selected for visualization changed (ParaView plugin)
LOGICAL                           :: changedFV_Elems             !< different distribution of DG and FV elements
LOGICAL                           :: changedWithDGOperator       !< If the DG operator should be called or not changed
LOGICAL                           :: changedDGonly               !< Visualize FV cells as DG changed
LOGICAL                           :: changedBCnames              !< BCnames selected for visualization changed (ParaView plugin)
LOGICAL                           :: changedAvg2D                !< mode changed between average solution and not

CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarNamesHDF5(:)         !< varnames in state file (DG_Solution, not including generic
                                                                 !< element- or pointwise)
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: VarnamesAll(:)          !< all available varnames (state file + dependent vars + generic)
INTEGER                               :: nVarAll                 !< number of all available visu variables
INTEGER                               :: nVarDep                 !< number of dependent variables, that EOS can calculate
INTEGER                               :: nVarVisu                !< number of variables selected for visualization
INTEGER,ALLOCATABLE                   :: mapAllVarsToVisuVars(:) !< maps all available variable index to visu variable index
INTEGER,ALLOCATABLE                   :: DepTable(:,:)           !< table holding the EOS dependencies required to calculate
                                                                 !< variables, that depend on other variables (e.g. primitive ...)
                                                                 !< The i-th line of this table holds the dependency informations of
                                                                 !< the i-th quantity on the previous quantities. The j-th column
                                                                 !< is 0 if the i-th quantity does NOT depends on the j-th quantity
                                                                 !< is 1 if the i-th quantity DOES depends on the j-th quantity

REAL,ALLOCATABLE,TARGET               :: UCalc_DG(:,:,:,:,:)     !< dependet variables require the computation of intermediate
                                                                 !< variables, that may not be visualized. Therefore the whole
                                                                 !< computation process takes place on this array and is afterwards
                                                                 !< converted to the visualization array (UVisu_DG)
REAL,ALLOCATABLE,TARGET               :: UCalc_FV(:,:,:,:,:)
INTEGER                               :: nVarCalc                !< number of (intermediate) variables that must be calculated
INTEGER,ALLOCATABLE                   :: mapDepToCalc(:)         !< maps all dependend variable index to calc variable index
INTEGER                               :: nVarCalc_FV             !< since FV reconstruction is done in primitive quantities, the
INTEGER,ALLOCATABLE                   :: mapDepToCalc_FV(:)      !< dependencies are different to the DG case, where everything is
                                                                 !< based on conservative quantities

REAL(C_DOUBLE),ALLOCATABLE            :: UVisu_DG(:,:,:,:,:)     !< DG solution that is written to VTK or send to ParaView
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: UVisu_FV(:,:,:,:,:)     !< FV solution that is written to VTK or send to ParaView
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsVisu_DG(:,:,:,:,:)!< coordinates of UVisu_DG
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsVisu_FV(:,:,:,:,:)!< coordinates of UVisu_FV
INTEGER,ALLOCATABLE,TARGET            :: nodeids_DG(:)           !< nodeIDs for CoordsVisu_DG
INTEGER,ALLOCATABLE,TARGET            :: nodeids_FV(:)           !< nodeIDs for CoordsVisu_FV


! ==============================================================================================================================
! Surface visualization
! ==============================================================================================================================
INTEGER,ALLOCATABLE                   :: DepSurfaceOnly(:)       !< Mask for quantities that
                                                                 !< are exclusively available on BCs
INTEGER,ALLOCATABLE                   :: DepVolumeOnly(:)        !< Mask for quantities that
                                                                 !< are exclusively available in the volume

INTEGER                               :: nBCNamesAll                  !< number of all BC names in mesh file
CHARACTER(LEN=255),ALLOCATABLE,TARGET :: BCNamesAll(:)                !< all BC names in mesh file
INTEGER                               :: nBCNamesVisu                 !< number of BC names selected for visualization
INTEGER,ALLOCATABLE                   :: mapAllBCNamesToVisuBCNames(:)!< maps global BCName index to visu BCName index
INTEGER,ALLOCATABLE                   :: mapAllBCNamesToVisuBCNames_old(:)

INTEGER                               :: nVarSurfVisuAll                 !< number of all avail. vars that are visualized on surf.
INTEGER,ALLOCATABLE                   :: mapAllVarsToSurfVisuVars(:)     !< maps all avail. var index to surf. visu. var index
INTEGER,ALLOCATABLE                   :: mapAllVarsToSurfVisuVars_old(:) !< saves previous mapAllVarsToSurfVisuVars
INTEGER                               :: nBCSidesVisu_DG                 !< number of DG BCsides selected for visualization
INTEGER                               :: nBCSidesVisu_FV                 !< number of FV BCsides selected for visualization
INTEGER,ALLOCATABLE                   :: mapAllBCSidesToDGVisuBCSides(:) !< map global BC side index to DG visu BC sides
INTEGER,ALLOCATABLE                   :: mapAllBCSidesToFVVisuBCSides(:) !< map global BC side index to FV visu BC sides
INTEGER,ALLOCATABLE                   :: nSidesPerBCNameVisu_DG(:)       !< holds number of DG BCsides for each BCName
INTEGER,ALLOCATABLE                   :: nSidesPerBCNameVisu_FV(:)       !< holds number of FV BCsides for each BCName

REAL,ALLOCATABLE                      :: USurfCalc_DG(:,:,:,:)        !< array on which dependent DG quantities are calculated
REAL,ALLOCATABLE                      :: USurfCalc_FV(:,:,:,:)        !< array on which dependent FV quantities are calculated

REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: USurfVisu_DG(     :,:,:,:,:) !< surf. DG solution written to VTK or send to ParaView
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: USurfVisu_FV(     :,:,:,:,:) !< surf. FV solution written to VTK or send to ParaView
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsSurfVisu_DG(:,:,:,:,:) !< coordinates of DG surface solution
REAL(C_DOUBLE),ALLOCATABLE,TARGET     :: CoordsSurfVisu_FV(:,:,:,:,:) !< coordinates of FV surface solution
INTEGER,ALLOCATABLE,TARGET            :: nodeidsSurf_DG(:)            !< nodeIDs for DG surface coordinates
INTEGER,ALLOCATABLE,TARGET            :: nodeidsSurf_FV(:)            !< nodeIDs for FV surface coordinates


! ==============================================================================================================================
! Avg2D
! ==============================================================================================================================
LOGICAL                           :: IJK_exists                  !< IJK sorting of elements
INTEGER,ALLOCATABLE               :: Elem_IJK(:,:)               !< IJK sorting of elements
INTEGER,ALLOCATABLE               :: Elem_IJK_glob(:,:)          !< IJK sorting of global elements for parallel avgerage
INTEGER                           :: nElems_IJK(3)               !< Number of elements in structured direction
REAL,ALLOCATABLE                  :: FVAmountAvg2D(:,:)          !< averaged FV_elems in z-direction at i,j-th element
INTEGER                           :: nElemsAvg2D_DG              !< number of cells averaged as DG cells
INTEGER                           :: nElemsAvg2D_FV              !< number of cells averaged as FV cells
INTEGER,ALLOCATABLE               :: mapElemIJToDGElemAvg2D(:,:) !< maps i,j element index to Avg2D DG element index
INTEGER,ALLOCATABLE               :: mapElemIJToFVElemAvg2D(:,:) !< maps i,j element index to Avg2D FV element index

!> Vandermonde matrixes to interpolate between DG <--> FV  and solution <--> visu grid
REAL,ALLOCATABLE                  :: Vdm_DGToFV  (:,:)
REAL,ALLOCATABLE                  :: Vdm_FVToDG  (:,:)
REAL,ALLOCATABLE                  :: Vdm_DGToVisu(:,:)
REAL,ALLOCATABLE                  :: Vdm_FVToVisu(:,:)

END MODULE MOD_Visu_Vars
