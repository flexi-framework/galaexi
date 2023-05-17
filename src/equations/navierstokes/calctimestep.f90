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
!> This module contains the routines to calculate the equation system specific allowable timestep.
!==================================================================================================================================
MODULE MOD_CalcTimeStep
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitCalctimestep
  MODULE PROCEDURE InitCalctimestep
END INTERFACE

INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE

INTERFACE FinalizeCalctimestep
  MODULE PROCEDURE FinalizeCalctimestep
END INTERFACE


PUBLIC :: InitCalctimestep,CALCTIMESTEP,FinalizeCalctimestep
!==================================================================================================================================

REAL,DEVICE,ALLOCATABLE :: MetricsAdv(:,:,:,:,:,:)  !< support variable: NORM2(Metricsfgh)/J
#if PARABOLIC
REAL,DEVICE,ALLOCATABLE :: MetricsVisc(:,:,:,:,:,:) !< support variable: kappa/Pr*(SUM((Metricsfgh/J)**2))
#endif

CONTAINS

!==================================================================================================================================
!> Precompute some metric support variables
!==================================================================================================================================
SUBROUTINE InitCalctimestep()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:d_sJ,d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde,nElems
#if PARABOLIC
USE MOD_EOS_Vars ,ONLY:KappasPr
#endif /*PARABOLIC*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem,FVE
#if PARABOLIC
REAL                         :: KappasPr_max
#endif /*PARABOLIC*/
!==================================================================================================================================

ALLOCATE(MetricsAdv(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
!$cuf kernel do(5) <<< *, * >>>
DO FVE=0,FV_SIZE
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      MetricsAdv(1,i,j,k,iElem,FVE)=d_sJ(i,j,k,iElem,FVE)*SQRT(DOT_PRODUCT(d_Metrics_fTilde(:,i,j,k,iElem,FVE),d_Metrics_fTilde(:,i,j,k,iElem,FVE)))
      MetricsAdv(2,i,j,k,iElem,FVE)=d_sJ(i,j,k,iElem,FVE)*SQRT(DOT_PRODUCT(d_Metrics_gTilde(:,i,j,k,iElem,FVE),d_Metrics_gTilde(:,i,j,k,iElem,FVE)))
      MetricsAdv(3,i,j,k,iElem,FVE)=d_sJ(i,j,k,iElem,FVE)*SQRT(DOT_PRODUCT(d_Metrics_hTilde(:,i,j,k,iElem,FVE),d_Metrics_hTilde(:,i,j,k,iElem,FVE)))
    END DO; END DO; END DO
  END DO
END DO
#if PARABOLIC
ALLOCATE(MetricsVisc(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
KappasPr_max=KAPPASPR_MAX_TIMESTEP_H()
!$cuf kernel do(5) <<< *, * >>>
DO FVE=0,FV_SIZE
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      MetricsVisc(1,i,j,k,iElem,FVE)=KappasPR_max*(SUM((d_Metrics_fTilde(:,i,j,k,iElem,FVE)*d_sJ(i,j,k,iElem,FVE))**2))
      MetricsVisc(2,i,j,k,iElem,FVE)=KappasPR_max*(SUM((d_Metrics_gTilde(:,i,j,k,iElem,FVE)*d_sJ(i,j,k,iElem,FVE))**2))
      MetricsVisc(3,i,j,k,iElem,FVE)=KappasPR_max*(SUM((d_Metrics_hTilde(:,i,j,k,iElem,FVE)*d_sJ(i,j,k,iElem,FVE))**2))
    END DO; END DO; END DO
  END DO
END DO
#endif /*PARABOLIC*/
END SUBROUTINE


!==================================================================================================================================
!> Compute the time step for the current update of U for the Navier-Stokes-Equations
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY:d_U
USE MOD_EOS_Vars
USE MOD_Mesh_Vars    ,ONLY:d_sJ,d_Metrics_fTilde,d_Metrics_gTilde,Elem_xGP,nElems
USE MOD_TimeDisc_Vars,ONLY:d_CFLScale,ViscousTimeStep,d_dtElem
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
#if PP_dim==3
USE MOD_Mesh_Vars    ,ONLY:d_Metrics_hTilde
#endif
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY:DFLScale
USE MOD_Viscosity
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
#endif
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars, ONLY: muSGS
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep
INTEGER,INTENT(OUT)          :: errType
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k,iElem
REAL,DIMENSION(PP_2Var),DEVICE :: UE
REAL,DEVICE                    :: TimeStepConv, TimeStepVisc
REAL                           :: TimeStep(3)
REAL,DEVICE                    :: Lambda(3),c,vsJ(3)
#if EDDYVISCOSITY
REAL                           :: muSGSmax
#endif
#if PARABOLIC
REAL                           :: Max_Lambda_v(3),mu,prim(PP_nVarPrim)
#endif /*PARABOLIC*/
INTEGER,DEVICE                 :: FVE
!==================================================================================================================================
errType=0

TimeStepConv=HUGE(1.)
TimeStepVisc=HUGE(1.)
!$cuf kernel do(4) <<< *, * >>> reduce(min:TimeStepConv)
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    ! TODO: ATTENTION: Temperature of UE not filled!!!
    UE(EXT_CONS)=d_U(:,i,j,k,iElem)
    UE(EXT_SRHO)=1./UE(EXT_DENS)
    UE(EXT_VELV)=VELOCITY_HE(UE)
    UE(EXT_PRES)=PRESSURE_HE(UE)
    UE(EXT_TEMP)=TEMPERATURE_HE(UE)
    c=SPEEDOFSOUND_HE(UE)
    vsJ=UE(EXT_VELV)*d_sJ(i,j,k,iElem,FV_Elems(iElem))
    Lambda(1)=ABS(SUM(d_Metrics_fTilde(:,i,j,k,iElem,FV_Elems(iElem))*vsJ)) + &
                          c*MetricsAdv(1,i,j,k,iElem,FV_Elems(iElem))
    Lambda(2)=ABS(SUM(d_Metrics_gTilde(:,i,j,k,iElem,FV_Elems(iElem))*vsJ)) + &
                          c*MetricsAdv(2,i,j,k,iElem,FV_Elems(iElem))
    Lambda(3)=ABS(SUM(d_Metrics_hTilde(:,i,j,k,iElem,FV_Elems(iElem))*vsJ)) + &
                          c*MetricsAdv(3,i,j,k,iElem,FV_Elems(iElem))
    TimeStepConv=d_CFLScale(FV_Elems(iElem))*2./SUM(Lambda)
  END DO; END DO; END DO ! i,j,k
END DO ! iElem=1,nElems

TimeStep(1)=TimeStepConv
TimeStep(2)=TimeStepVisc
#if USE_MPI
TimeStep(3)=-errType ! reduce with timestep, minus due to MPI_MIN
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStep,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_FLEXI,iError)
errType=INT(-TimeStep(3))
#endif /*USE_MPI*/
ViscousTimeStep=(TimeStep(2) .LT. TimeStep(1))
CalcTimeStep=MINVAL(TimeStep(1:2))

END FUNCTION CALCTIMESTEP


!==================================================================================================================================
!> Deallocate CalcTimeStep arrays
!==================================================================================================================================
SUBROUTINE FinalizeCalctimestep()
! MODULES
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(MetricsAdv)
#if PARABOLIC
SDEALLOCATE(MetricsVisc)
#endif
END SUBROUTINE

END MODULE MOD_CalcTimeStep
