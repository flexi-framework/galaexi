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
!==================================================================================================================================
!> Variables needed for the evaluation of the record points
!==================================================================================================================================
MODULE MOD_RecordPoints_Vars
! MODULES
#if USE_MPI
USE mpi
#endif

IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)              :: RPDefFile               !< File containing element-local parametric recordpoint coordinates and structure
LOGICAL                         :: RecordPointsInitIsDone = .FALSE. !< mark wheter recordpoints init routine is finished
LOGICAL                         :: RP_inUse  = .FALSE.     !< mark whether recordpoints should be evaluated during computation
LOGICAL                         :: RP_onProc = .FALSE.     !< marks wheter current proc has RPs
INTEGER                         :: RP_BufferSize           !< no. of time samples (size of RP_Data)
INTEGER                         :: RP_MaxBuffersize        !< max. allowed no. of time samples
INTEGER                         :: RP_SamplingOffset       !< sampling rate (each .. iterations)
INTEGER                         :: nRP                     !< no. of RP on proc
INTEGER                         :: nGlobalRP               !< total no. of RP
INTEGER                         :: offsetRP                !< offset for each proc in global RP list
INTEGER                         :: iSample=0               !< no of samples in array
INTEGER,ALLOCATABLE             :: RP_ElemID(:)            !< mapping from RP->Elem (nRP)
!@cuf INTEGER,DEVICE,ALLOCATABLE:: d_RP_ElemID(:)
#if FV_ENABLED
INTEGER,ALLOCATABLE             :: FV_RP_ijk(:,:)          !< ijk-index of FV subcell nearest to record point [1:3,nRP]
#endif
REAL,ALLOCATABLE                :: L_xi_RP(:,:)            !< Lagrange basis evaluated at RP coords (xi-dir)
REAL,ALLOCATABLE                :: L_eta_RP(:,:)           !< Lagrange basis evaluated at RP coords (eta-dir)
REAL,ALLOCATABLE                :: L_zeta_RP(:,:)          !< Lagrange basis evaluated at RP coords (zeta-dir)
!@cuf REAL,DEVICE,ALLOCATABLE   :: d_L_xi_RP(:,:)
!@cuf REAL,DEVICE,ALLOCATABLE   :: d_L_eta_RP(:,:)
!@cuf REAL,DEVICE,ALLOCATABLE   :: d_L_zeta_RP(:,:)
!@cuf REAL,DEVICE,ALLOCATABLE   :: d_U_RP(:,:)
REAL,ALLOCATABLE                :: RP_Data(:,:,:)          !< solution evaluated at RPs (nvar,nRP,RP_BufferSize)
REAL,ALLOCATABLE                :: lastSample(:,:)         !< solution evaluated at RPs (nvar,nRP,RP_BufferSize)
CHARACTER(LEN=255)              :: StrVarNames(PP_nVar)    !< RP variables names for output

!----------------------------------------------------------------------------------------------------------------------------------
! MPI Communicator for RPs
!----------------------------------------------------------------------------------------------------------------------------------
#if USE_MPI
INTEGER                         :: myRPrank                !< rank within RP communicator
INTEGER                         :: RP_COMM=MPI_COMM_NULL   !< MPI RP communicator
INTEGER                         :: nRP_Procs               !< number of procs with RPs
#endif /* USE_MPI */

END MODULE MOD_recordPoints_Vars
