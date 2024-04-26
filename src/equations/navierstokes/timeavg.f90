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
#include "eos.h"

!==================================================================================================================================
!> Routine performing time averaging of variables and the preparation to computing fluctuations
!> For the fluctuations we rely on the fact that \f$ U^{'} U^{'} = \overline{U}^2 - \overline{U^2} \f$
!> The terms computed in this routine are therefore the TimeAvg: \f$ \overline{U} \f$ and
!> the squared solution denoted by Fluc: \f$ \overline{U^2} \f$
!==================================================================================================================================
MODULE MOD_TimeAverage
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER                        :: nMaxVarAvg,nMaxVarFluc
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitCalcTimeAverage
  MODULE PROCEDURE InitCalcTimeAverage
END INTERFACE

INTERFACE FinalizeTimeAverage
  MODULE PROCEDURE FinalizeTimeAverage
END INTERFACE

INTERFACE CalcTimeAverage
  MODULE PROCEDURE CalcTimeAverage
END INTERFACE

PUBLIC::InitCalcTimeAverage, FinalizeTimeAverage, CalcTimeAverage
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Initializes the time averaging variables and builds map from fluctuation quantities to required time averaged variables
!> (e.g. if VelocityMagnitude fluctuations are to be computed the time averages of the velocity components u,v,w are computed)
!> - only variables specified in the variable list can be averaged
!==================================================================================================================================
SUBROUTINE InitCalcTimeAverage()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools, ONLY: CountOption,GETSTR,GETLOGICAL,GETINT
USE MOD_Mesh_Vars,   ONLY: nElems
USE MOD_AnalyzeEquation_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVar,iVar2
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAvgIni(:), VarNamesAvgList(:), VarNamesFlucList(:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesFlucIni(:)
LOGICAL,ALLOCATABLE            :: hasAvgVars(:)
!==================================================================================================================================
nVarAvg  = CountOption('VarNameAvg')
nVarFluc = CountOption('VarNameFluc')
IF((nVarAvg.EQ.0).AND.(nVarFluc.EQ.0))THEN
  CALL CollectiveStop(__STAMP__, &
    'No quantities for time averaging have been specified. Please specify quantities or disable time averaging!')
END IF
#if FV_ENABLED
SWRITE(UNIT_stdOut,'(A)') 'Warning: If FV is enabled, time averaging is performed on integral cell mean values.'
#endif

! Define variables to be averaged
nMaxVarAvg=15
ALLOCATE(VarNamesAvgList(nMaxVarAvg))
VarNamesAvgList(1)  ='Density'
VarNamesAvgList(2)  ='MomentumX'
VarNamesAvgList(3)  ='MomentumY'
VarNamesAvgList(4)  ='MomentumZ'
VarNamesAvgList(5)  ='EnergyStagnationDensity'
VarNamesAvgList(6)  ='VelocityX'
VarNamesAvgList(7)  ='VelocityY'
VarNamesAvgList(8)  ='VelocityZ'
VarNamesAvgList(9)  ='VelocityMagnitude'
VarNamesAvgList(10) ='Pressure'
VarNamesAvgList(11) ='VelocitySound'
VarNamesAvgList(12) ='Mach'
VarNamesAvgList(13) ='Temperature'
VarNamesAvgList(14) ='TotalTemperature'
VarNamesAvgList(15) ='TotalPressure'

nMaxVarFluc=21
ALLOCATE(VarNamesFlucList(nMaxVarFluc),hasAvgVars(nMaxVarFluc))
hasAvgVars=.TRUE.
!define fluctuation variables
VarNamesFlucList(1)  ='Density'
VarNamesFlucList(2)  ='MomentumX'
VarNamesFlucList(3)  ='MomentumY'
VarNamesFlucList(4)  ='MomentumZ'
VarNamesFlucList(5)  ='EnergyStagnationDensity'
VarNamesFlucList(6)  ='VelocityX'
VarNamesFlucList(7)  ='VelocityY'
VarNamesFlucList(8)  ='VelocityZ'
VarNamesFlucList(9)  ='VelocityMagnitude'
VarNamesFlucList(10) ='Pressure'
VarNamesFlucList(11) ='VelocitySound'
VarNamesFlucList(12) ='Mach'
VarNamesFlucList(13) ='Temperature'
!VarNamesFlucList(14) ='EnthalpyStagnation'
!VarNamesFlucList(15) ='Entropy'
!VarNamesFlucList(16) ='VorticityX'
!VarNamesFlucList(17) ='VorticityY'
!VarNamesFlucList(18) ='VorticityZ'
!VarNamesFlucList(19) ='VorticityMagnitude'
VarNamesFlucList(14) = 'uv'
VarNamesFlucList(15) = 'uw'
VarNamesFlucList(16) = 'vw'
VarNamesFlucList(17) = 'DR_u'; hasAvgVars(17)=.FALSE.
VarNamesFlucList(18) = 'DR_S'; hasAvgVars(18)=.FALSE.
VarNamesFlucList(19) = 'TKE';  hasAvgVars(19)=.FALSE.
VarNamesFlucList(20) = 'TotalTemperature'
VarNamesFlucList(21) = 'TotalPressure'

! Read VarNames from ini file
ALLOCATE(VarNamesAvgIni(nVarAvg),VarNamesFlucIni(nVarFluc))
DO iVar=1,nVarAvg
  VarNamesAvgIni(iVar)=GETSTR('VarNameAvg')
END DO
DO iVar=1,nVarFluc
  VarNamesFlucIni(iVar)=GETSTR('VarNameFluc')
END DO

! Check which variables have to be calculated and create mappings to global variable index (1:nVarout)
! CalcAvgTmp(1,:) for normal variables, CalcAvgTmp(2,:) for fluctuations
ALLOCATE(CalcAvg (nMaxVarAvg ))
ALLOCATE(CalcFluc(nMaxVarFluc))
CalcAvg  = .FALSE.
CalcFluc = .FALSE.

! check each average from ini file
DO iVar=1,nVarAvg
  ! check if avg from ini file is in avg list
  iVar2 = GETMAPBYNAME(VarNamesAvgIni(iVar),VarNamesAvgList,nMaxVarAvg)
  IF(iVar2.NE.-1)THEN
    CalcAvg(iVar2) = .TRUE.
  ELSE
    CALL CollectiveStop(__STAMP__, &
    'Specified varname does not exist: ' // VarNamesAvgIni(iVar))
  END IF
END DO

! check each fluctuation from ini file
DO iVar=1,nVarFluc
  ! check if fluc from ini file is in fluc list
  iVar2 = GETMAPBYNAME(VarNamesFlucIni(iVar),VarNamesFlucList,nMaxVarFluc)
  IF(iVar2.NE.-1)THEN
    CalcFluc(iVar2) = .TRUE.
  ELSE
    CALL CollectiveStop(__STAMP__, &
    'Specified varname does not exist: ' // VarNamesFlucIni(iVar))
  END IF

  ! if fluctuation is set also compute base variable
  iVar2 = GETMAPBYNAME(VarNamesFlucIni(iVar),VarNamesAvgList,nMaxVarAvg)
  IF(iVar2.NE.-1) CalcAvg(iVar2) = .TRUE.
END DO

! For fluctuations with mixed base vars
IF(CalcFluc(GETMAPBYNAME('uv',VarNamesFlucList,nMaxVarFluc)))THEN !uv
   CalcAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF
IF(CalcFluc(GETMAPBYNAME('uw',VarNamesFlucList,nMaxVarFluc)))THEN !uw
   CalcAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF
IF(CalcFluc(GETMAPBYNAME('vw',VarNamesFlucList,nMaxVarFluc)))THEN !vw
   CalcAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF
IF(CalcFluc(GETMAPBYNAME('DR_u',VarNamesFlucList,nMaxVarFluc)).OR.&
        CalcFluc(GETMAPBYNAME('DR_S',VarNamesFlucList,nMaxVarFluc)).OR.&
        CalcFluc(GETMAPBYNAME('TKE',VarNamesFlucList,nMaxVarFluc)))THEN !Dissipation,TKE
   CalcAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
   CalcAvg(GETMAPBYNAME  ('Density',VarNamesAvgList,nMaxVarAvg)) = .TRUE.
END IF

nVarAvg=0 ! recount nVarAvg
DO iVar=1,nMaxVarAvg
  IF(CalcAvg(iVar)) nVarAvg=nVarAvg+1
END DO

! IF(ANY(CalcAvg(6:nMaxVarAvg))) calcPrims = .TRUE.
! IF(ANY(CalcAvg(11:15)))        calcEOS   = .TRUE.

! Set indices (iAvg,iFluc) and Varnames
ALLOCATE(VarNamesFlucOut(nVarFluc),VarNamesAvgOut(nVarAvg))
ALLOCATE(iAvg(nMaxVarAvg),iFluc(nMaxVarFluc))
! iAvg     -> Mapping from VariableList to index in UAvg array
! iFluc    -> Mapping from index in UFluc array to index in UAvg array
!             (e.g. for mixed term uv: iFluc(1,1) -> u iFluc(2,1) -> v)

VarNamesFlucOut(:) = ''
VarNamesAvgOut(:)  = ''
nVarAvg  = 0
nVarFluc = 0
iAvg     = 0
iFluc    = 0
! Build map for avg
DO iVar=1,nMaxVarAvg
  IF(CalcAvg(iVar))THEN
    nVarAvg=nVarAvg+1
    iAvg(iVar)=nVarAvg
    VarNamesAvgOut(nVarAvg) = TRIM(VarNamesAvgList(iVar))
  END IF
END DO
! Build map from fluclist to calcfluc
DO iVar=1,nMaxVarFluc
  IF(CalcFluc(iVar).AND.hasAvgVars(iVar))THEN
    nVarFluc=nVarFluc+1
    iFluc(iVar)=nVarFluc
    VarNamesFlucOut(nVarFluc) = TRIM(VarNamesFlucList(iVar))
  END IF
END DO
nVarFlucHasAvg=nVarFluc
ALLOCATE(FlucAvgMap(2,nVarFlucHasAvg))
FlucAvgMap=0
DO iVar=1,nMaxVarFluc
  IF(CalcFluc(iVar).AND.(.NOT.hasAvgVars(iVar)))THEN
    nVarFluc=nVarFluc+1
    iFluc(iVar)=nVarFluc
    VarNamesFlucOut(nVarFluc) = TRIM(VarNamesFlucList(iVar))
  END IF
END DO

! set map from fluc array to avg array needed to compute fluc
DO iVar=1,nMaxVarFluc
  IF((iFluc(iVar).NE.0).AND.hasAvgVars(iVar))THEN
    iVar2 = GETMAPBYNAME(VarNamesFlucList(iVar),VarNamesAvgList,nMaxVarAvg)
    IF(iVar2.GT.0) FlucAvgMap(:,iFluc(iVar))=iAvg(iVar2)
    IF(iVar.EQ.GETMAPBYNAME('uv',VarNamesFlucList,nMaxVarFluc)) THEN !uv
      FlucAvgMap(1,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg))
      FlucAvgMap(2,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg))
    END IF
    IF(iVar.EQ.GETMAPBYNAME('vw',VarNamesFlucList,nMaxVarFluc)) THEN !vw
      FlucAvgMap(1,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityY',VarNamesAvgList,nMaxVarAvg))
      FlucAvgMap(2,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg))
    END IF
    IF(iVar.EQ.GETMAPBYNAME('uw',VarNamesFlucList,nMaxVarFluc)) THEN !uw
      FlucAvgMap(1,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityX',VarNamesAvgList,nMaxVarAvg))
      FlucAvgMap(2,iFluc(iVar))=iAvg(GETMAPBYNAME('VelocityZ',VarNamesAvgList,nMaxVarAvg))
    END IF
  END IF
END DO

#if !(PARABOLIC)
IF(CalcFluc(17).OR.CalcFluc(18))THEN
  CALL CollectiveStop(__STAMP__,&
    'Cannot compute dissipation. Not compiled with parabolic flag.')
END IF
#endif /* PARABOLIC */

! Allocate arrays
ALLOCATE(UAvg( nVarAvg ,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(UFluc(nVarFluc,0:PP_N,0:PP_N,0:PP_NZ,nElems))
UAvg  = 0.
UFluc = 0.
dtOld = 0.
dtAvg = 0.

IF (nVarAvg.GT.0) THEN
!@cuf ALLOCATE(d_UAvg( nVarAvg ,0:PP_N,0:PP_N,0:PP_N,nElems))
  d_UAvg  = UAvg
END IF ! nVarAvg.GT.0
IF (nVarFluc.GT.0) THEN
!@cuf ALLOCATE(d_UFluc(nVarFluc,0:PP_N,0:PP_N,0:PP_N,nElems))
  d_UFluc = UFluc
END IF ! nVarFluc.GT.0

IF (nMaxVarAvg.GT.0) THEN
!@cuf ALLOCATE(d_CalcAvg (nMaxVarAvg ))
!@cuf ALLOCATE(d_iAvg(nMaxVarAvg))
  d_CalcAvg  = CalcAvg
  d_iAvg     = iAvg
END IF ! nMaxVarAvg.GT.0
IF (nMaxVarFluc.GT.0) THEN
!@cuf ALLOCATE(d_CalcFluc(nMaxVarFluc))
!@cuf ALLOCATE(d_iFluc(nMaxVarFluc))
  d_CalcFluc = CalcFluc
  d_iFluc = iFluc
END IF ! nMaxVarFluc.GT.0

IF (nVarFlucHasAvg.GT.0) THEN
!@cuf ALLOCATE(d_FlucAvgMap(2,nVarFlucHasAvg))
  d_FlucAvgMap = FlucAvgMap
END IF

DEALLOCATE(VarNamesAvgList,VarNamesAvgIni,VarNamesFlucIni)
DEALLOCATE(VarNamesFlucList)
END SUBROUTINE InitCalcTimeAverage


!==================================================================================================================================
!> Return index of string VarName in array VarNameList
!==================================================================================================================================
PURE FUNCTION GETMAPBYNAME(VarName,VarNameList,nVarList)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: VarName                 !< string to be compared
INTEGER,INTENT(IN)             :: nVarList                !< length of list
CHARACTER(LEN=*),INTENT(IN)    :: VarNameList(nVarList)   !< list of strings to be searched
INTEGER                        :: GETMAPBYNAME            !< index of VarName in VarNameList
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i
!==================================================================================================================================
GETMAPBYNAME=-1
DO i=1,nVarList
  IF(TRIM(VarName).EQ.TRIM(VarNameList(i)))THEN
    GETMAPBYNAME=i
    RETURN
  END IF
END DO
END FUNCTION


!==================================================================================================================================
!> Compute time averages by trapezoidal rule
!> TODO: extend description
!==================================================================================================================================
SUBROUTINE CalcTimeAverage(Finalize,dt,t,streamID)
! MODULES
USE MOD_Globals
USE MOD_GPU
USE CUDAFOR
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: d_U,d_UPrim
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: d_GradUx,d_GradUy,d_GradUz
#endif
USE MOD_Mesh_Vars    ,ONLY: MeshFile,nElems
USE MOD_HDF5_Output  ,ONLY: WriteTimeAverage
USE MOD_EOS          ,ONLY: ConsToPrim
USE MOD_EOS_Vars     ,ONLY: Kappa
USE MOD_Analyze_Vars ,ONLY: WriteData_dt
USE MOD_AnalyzeEquation_Vars
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems,FV_Vdm
USE MOD_ChangeBasisByDim,ONLY:ChangeBasisVolume
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN)                                 :: Finalize               !< finalized trapezoidal rule and output file
REAL,INTENT(IN)                                    :: dt                     !< current timestep for averaging
REAL,INTENT(IN),OPTIONAL                           :: t                      !< current simulation time
INTEGER(KIND=CUDA_STREAM_KIND),OPTIONAL,INTENT(IN) :: streamID
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                               :: tFuture
REAL                                               :: dtStep
REAL,POINTER                                       :: Uloc(:,:,:,:)
INTEGER                                            :: FV_Elems_loc(1:nElems)
INTEGER                                            :: nThreads
INTEGER(KIND=CUDA_STREAM_KIND)                     :: mystream
INTEGER                                            :: NAvg
!----------------------------------------------------------------------------------------------------------------------------------
dtStep = (dtOld+dt)*0.5
IF(Finalize) dtStep = dt*0.5
dtAvg  = dtAvg+dtStep
dtOld  = dt

mystream=DefaultStream
IF (PRESENT(streamID)) mystream=streamID

iError=CudaDeviceSynchronize()

! PP_N might be fixed pre-processor, pass it to NAvg
NAvg = PP_N

CALL ConsToPrim(NAvg,d_UPrim,d_U,streamID=mystream)
nThreads=(NAvg+1)**3
CALL TimeAvg<<<nElems,nThreads,0,mystream>>>(NAvg,nElems,Kappa,dtStep,d_U,d_UPrim,nVarAvg,nVarFluc,nMaxVarAvg,nMaxVarFluc, &
             nVarFlucHasAvg,d_CalcAvg,d_CalcFluc,d_iAvg,d_iFluc,d_FlucAvgMap,d_UAvg,d_UFluc                                &
#if PARABOLIC
            ,d_GradUx,d_GradUy,d_GradUz                                                                                    &
#endif /*PARABOLIC*/
            )

! Calc time average and write solution to file
IF(Finalize)THEN
  ! DEVICE SYNCHRONIZE: All streams have to be ready, cannot progress on unfiltered solution
  iError=CudaDeviceSynchronize()

  UAvg  = d_UAvg
  UFluc = d_UFluc

  IF(nVarAvg .GT.0) UAvg  = UAvg /dtAvg
  IF(nVarFluc.GT.0) UFluc = UFluc/dtAvg
  tFuture      = t + WriteData_dt
  FV_Elems_loc = FV_ENABLED
  CALL WriteTimeAverage(MeshFile,t,dtAvg,FV_Elems_loc,(/PP_N+1,PP_N+1,PP_NZ+1/),&
                        nVarAvg ,VarNamesAvgOut ,UAvg ,&
                        nVarFluc,VarNamesFlucOut,UFluc,&
                        FutureTime=tFuture)
  d_UAvg  = 0.
  d_UFluc = 0.
  dtAvg = 0.
  dtOld = 0.
END IF

END SUBROUTINE CalcTimeAverage


!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE TimeAvg(NAvg,nElems,Kappa,dtStep,U,UPrim,nVarAvg,nVarFluc,nMaxVarAvg,nMaxVarFluc,nVarFlucHasAvg     &
                                     ,CalcAvg,CalcFluc,iAvg,iFluc,FlucAvgMap,UAvg,UFluc                                           &
#if PARABOLIC
                                     ,GradUx,GradUy,GradUz                                                                        &
#endif /*PARABOLIC*/
                                     )
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)        :: NAvg                                      !< Polynomial degree
INTEGER,VALUE,INTENT(IN)        :: nElems
REAL,VALUE,INTENT(IN)           :: Kappa
REAL,VALUE,INTENT(IN)           :: dtStep
REAL,DEVICE,INTENT(IN)          :: U(    CONS,0:NAvg,0:NAvg,0:NAvg,1:nElems)
REAL,DEVICE,INTENT(IN)          :: UPrim(PRIM,0:NAvg,0:NAvg,0:NAvg,1:nElems)
INTEGER,VALUE,INTENT(IN)        :: nVarAvg
INTEGER,VALUE,INTENT(IN)        :: nVarFluc
INTEGER,VALUE,INTENT(IN)        :: nMaxVarAvg
INTEGER,VALUE,INTENT(IN)        :: nMaxVarFluc
INTEGER,VALUE,INTENT(IN)        :: nVarFlucHasAvg
LOGICAL,DEVICE,INTENT(IN)       :: CalcAvg (nMaxVarAvg)
LOGICAL,DEVICE,INTENT(IN)       :: CalcFluc(nMaxVarFluc)
INTEGER,DEVICE,INTENT(IN)       :: iAvg(nMaxVarAvg)
INTEGER,DEVICE,INTENT(IN)       :: iFluc(nMaxVarFluc)
INTEGER,DEVICE,INTENT(IN)       :: FlucAvgMap(2,nVarFlucHasAvg)
REAL,DEVICE,INTENT(INOUT)       :: UAvg( nVarAvg, 0:NAvg,0:NAvg,0:NAvg,1:nElems)
REAL,DEVICE,INTENT(INOUT)       :: UFluc(nVarFluc,0:NAvg,0:NAvg,0:NAvg,1:nElems)
#if PARABOLIC
REAL,DEVICE,INTENT(IN)          :: GradUx(PP_nVarLifting,0:NAvg,0:NAvg,0:NAvg,1:nElems)
REAL,DEVICE,INTENT(IN)          :: GradUy(PP_nVarLifting,0:NAvg,0:NAvg,0:NAvg,1:nElems)
REAL,DEVICE,INTENT(IN)          :: GradUz(PP_nVarLifting,0:NAvg,0:NAvg,0:NAvg,1:nElems)
#endif /*PARABOLIC*/
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem,threadID,rest
REAL                            :: UE(PP_2Var)
REAL                            :: tmpVars(nVarAvg)
REAL                            :: a,Mach
#if PARABOLIC
INTEGER                         :: p,q
REAL                            :: GradVel(1:3,1:3),Shear(1:3,1:3)
#endif
!----------------------------------------------------------------------------------------------------------------------------------

! Get thread indices
threadID = (blockidx%x-1) * blockdim%x + threadidx%x
! Get ElemID of current thread
iElem   =        (threadID-1)/((NAvg+1)**3)+1 ! Elems are 1-indexed
rest    =  threadID-(iElem-1)*((NAvg+1)**3)
! Get jk indices of current thread
k       = (rest-1)/(NAvg+1)**2
rest    =  rest- k*(NAvg+1)**2
j       = (rest-1)/(NAvg+1)!**1
rest    =  rest- j*(NAvg+1)!**1
i       = (rest-1)!/(Nloc+1)**0
rest    =  rest- i!*(Nloc+1)**0

IF (iElem.LE.nElems) THEN
  ! Compute speed of sound
  UE(EXT_CONS) = U(    :,i,j,k,iElem)
  UE(EXT_PRIM) = UPrim(:,i,j,k,iElem)
  UE(EXT_SRHO) = 1./UPrim(DENS,i,j,k,iElem)
  a    = SPEEDOFSOUND_HE(UE)
  Mach = NORM2(UPrim(2:4,i,j,k,iElem))/a

#if PARABOLIC
  GradVel(:,1) = GradUx(LIFT_VELV,i,j,k,iElem)
  GradVel(:,2) = GradUy(LIFT_VELV,i,j,k,iElem)
  GradVel(:,3) = GradUz(LIFT_VELV,i,j,k,iElem)
#endif /*PARABOLIC*/

  ! Compute time averaged variables and fluctuations of these variables
  IF(CalcAvg(1)) &  !'Density'
    tmpVars(iAvg(1)) = U(    DENS,i,j,k,iElem)

  IF(CalcAvg(2)) &  !'MomentumX'
    tmpVars(iAvg(2)) = U(    MOM1,i,j,k,iElem)

  IF(CalcAvg(3)) &  !'MomentumY'
    tmpVars(iAvg(3)) = U(    MOM2,i,j,k,iElem)

  IF(CalcAvg(4)) &  !'MomentumZ'
    tmpVars(iAvg(4)) = U(    MOM3,i,j,k,iElem)

  IF(CalcAvg(5)) &  !'EnergyStagnationDensity'
    tmpVars(iAvg(5)) = U(    ENER,i,j,k,iElem)

  IF(CalcAvg(6)) &  !'VelocityX'
    tmpVars(iAvg(6)) = UPrim(VEL1,i,j,k,iElem)

  IF(CalcAvg(7)) &  !'VelocityY'
    tmpVars(iAvg(7)) = UPrim(VEL2,i,j,k,iElem)

  IF(CalcAvg(8)) &  !'VelocityZ'
    tmpVars(iAvg(8)) = UPrim(VEL3,i,j,k,iElem)

  IF(CalcAvg(9)) &  !'VelocityMagnitude'
    tmpVars(iAvg(9)) = NORM2(UPrim(VELV,i,j,k,iElem))

  IF(CalcAvg(10)) & !'Pressure'
    tmpVars(iAvg(10)) = UPrim(PRES,i,j,k,iElem)

  IF(CalcAvg(11)) & !'VelocitySound'
    tmpVars(iAvg(11)) = a

  IF(CalcAvg(12)) & !'Mach'
    tmpVars(iAvg(12)) = Mach

  IF(CalcAvg(13)) & !'Temperature'
    tmpVars(iAvg(13)) = UPrim(TEMP,i,j,k,iElem)

  IF(CalcAvg(14)) & !'TotalTemperature'
    tmpVars(iAvg(14)) = TOTAL_TEMPERATURE_H(UE(EXT_TEMP),Mach)

  IF(CalcAvg(15))  & !'TotalPressure
    tmpVars(iAvg(15)) = TOTAL_PRESSURE_H(UE(EXT_PRES),Mach)

  UAvg( :                 ,i,j,k,iElem) = UAvg ( :              ,i,j,k,iElem) + tmpVars(1:nVarAvg                     )*dtStep
  IF(nVarFlucHasAvg.GT.0) &
    UFluc(1:nVarFlucHasAvg,i,j,k,iElem) = UFluc(1:nVarFlucHasAvg,i,j,k,iElem) + tmpVars(FlucAvgMap(1,1:nVarFlucHasAvg)) &
                                                                              * tmpVars(FlucAvgMap(2,1:nVarFlucHasAvg))*dtStep
#if PARABOLIC
  IF(CalcFluc(17)) THEN ! DR_u
    Shear = 0.5*(Gradvel+TRANSPOSE(GradVel))
    DO p=1,3; DO q=1,3
      UFluc(iFluc(17),i,j,k,iElem) = UFluc(iFluc(17),i,j,k,iElem) + Shear(p,q)*Shear(p,q)*dtStep
    END DO; END DO
  END IF

  IF(CalcFluc(18)) THEN ! DR_S
    DO p=1,3; DO q=1,3
      UFluc(iFluc(18),i,j,k,iElem) = UFluc(iFluc(18),i,j,k,iElem) + GradVel(p,q)*GradVel(p,q)*dtStep
    END DO; END DO
  END IF
#endif /* PARABOLIC */

  IF(CalcFluc(19)) &  !'TKE'
    UFluc(iFluc(19),i,j,k,iElem)   = UFluc(iFluc(19),i,j,k,iElem) + SUM(UPrim(VELV,i,j,k,iElem)**2)*dtStep

END IF ! iElem

END SUBROUTINE TimeAvg


!==================================================================================================================================
!> Finalizes the time averaging routines
!==================================================================================================================================
SUBROUTINE FinalizeTimeAverage()
! MODULES
USE MOD_AnalyzeEquation_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(CalcAvg)
!@cuf SDEALLOCATE(d_CalcAvg)
SDEALLOCATE(CalcFluc)
!@cuf SDEALLOCATE(d_CalcFluc)
SDEALLOCATE(iAvg)
!@cuf SDEALLOCATE(d_iAvg)
SDEALLOCATE(iFluc)
!@cuf SDEALLOCATE(d_iFluc)
SDEALLOCATE(UAvg)
!@cuf SDEALLOCATE(d_UAvg)
SDEALLOCATE(UFluc)
!@cuf SDEALLOCATE(d_UFluc)
SDEALLOCATE(VarNamesAvgOut)
!@cuf SDEALLOCATE(d_FlucAvgMap)
SDEALLOCATE(VarNamesFlucOut)
SDEALLOCATE(FlucAvgMap)

END SUBROUTINE FinalizeTimeAverage

END MODULE MOD_TimeAverage
