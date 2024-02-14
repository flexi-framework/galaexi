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

REAL,DEVICE,ALLOCATABLE :: d_MetricsAdv(:,:,:,:,:,:)  !< support variable: NORM2(Metricsfgh)/J
#if PARABOLIC
REAL,DEVICE,ALLOCATABLE :: d_MetricsVisc(:,:,:,:,:,:) !< support variable: kappa/Pr*(SUM((Metricsfgh/J)**2))
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

ALLOCATE(d_MetricsAdv(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
!$cuf kernel do(5) <<< *, * >>>
DO FVE=0,FV_SIZE
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      d_MetricsAdv(1,i,j,k,iElem,FVE)=d_sJ(i,j,k,iElem,FVE)*SQRT(DOT_PRODUCT(d_Metrics_fTilde(:,i,j,k,iElem,FVE),d_Metrics_fTilde(:,i,j,k,iElem,FVE)))
      d_MetricsAdv(2,i,j,k,iElem,FVE)=d_sJ(i,j,k,iElem,FVE)*SQRT(DOT_PRODUCT(d_Metrics_gTilde(:,i,j,k,iElem,FVE),d_Metrics_gTilde(:,i,j,k,iElem,FVE)))
      d_MetricsAdv(3,i,j,k,iElem,FVE)=d_sJ(i,j,k,iElem,FVE)*SQRT(DOT_PRODUCT(d_Metrics_hTilde(:,i,j,k,iElem,FVE),d_Metrics_hTilde(:,i,j,k,iElem,FVE)))
    END DO; END DO; END DO
  END DO
END DO
#if PARABOLIC
ALLOCATE(d_MetricsVisc(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_SIZE))
KappasPr_max=KAPPASPR_MAX_TIMESTEP_H()
!$cuf kernel do(5) <<< *, * >>>
DO FVE=0,FV_SIZE
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      d_MetricsVisc(1,i,j,k,iElem,FVE)=KappasPR_max*(SUM((d_Metrics_fTilde(:,i,j,k,iElem,FVE)*d_sJ(i,j,k,iElem,FVE))**2))
      d_MetricsVisc(2,i,j,k,iElem,FVE)=KappasPR_max*(SUM((d_Metrics_gTilde(:,i,j,k,iElem,FVE)*d_sJ(i,j,k,iElem,FVE))**2))
      d_MetricsVisc(3,i,j,k,iElem,FVE)=KappasPR_max*(SUM((d_Metrics_hTilde(:,i,j,k,iElem,FVE)*d_sJ(i,j,k,iElem,FVE))**2))
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
USE MOD_EOS_Vars     ,ONLY:d_EOS_Vars
USE MOD_Mesh_Vars    ,ONLY:d_sJ,d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde,nElems
USE MOD_TimeDisc_Vars,ONLY:d_CFLScale,ViscousTimeStep,d_dtElem
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY:d_DFLScale
#endif /*PARABOLIC*/
#if FV_ENABLED==2
USE MOD_FV_Vars,ONLY: d_FV_alpha,FV_alpha_min
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
INTEGER                        :: iElem,SHMEM_SIZE
REAL                           ::   Timestep(3),dummy
REAL,DEVICE                    :: d_Timestep(3)
REAL,DEVICE                    :: d_Lambda_max(2,nElems)
!==================================================================================================================================
errType=0
d_TimeStep = HUGE(1.)

#if PARABOLIC
SHMEM_SIZE=(PP_N+1)**2*6*SIZEOF(dummy) ! No of entries times size of real
#else
SHMEM_SIZE=(PP_N+1)**2*3*SIZEOF(dummy) ! No of entries times size of real
#endif

! Compute maximum eigenvalue within each element
CALL CalcMaxEigenvalue<<<nElems,(PP_N+1)**2,SHMEM_SIZE>>>(PP_N,nElems,d_U,d_EOS_Vars,&
                                                          d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde,&
                                                          d_MetricsAdv,&
#if PARABOLIC
                                                          d_MetricsVisc,&
#endif
                                                          d_sJ,&
                                                          d_Lambda_max)
! Compute minimum timestep
!$cuf kernel do <<< *, * >>>
DO iElem=1,nElems
#if FV_ENABLED == 2
  IF (d_FV_alpha(iElem) .LT. FV_alpha_min) THEN
    d_TimeStep(1)=MIN(d_TimeStep(1),d_CFLScale(0)*2./d_Lambda_max(1,iElem))
#if PARABOLIC
    d_TimeStep(2)=MIN(d_TimeStep(2),d_DFLScale(0)*4./d_Lambda_max(2,iElem))
#endif
  ELSE
#endif /* FV_ENABLED == 2*/
    d_TimeStep(1)=MIN(d_TimeStep(1),MINVAL(d_CFLScale(:))*2./d_Lambda_max(1,iElem))
#if PARABOLIC
    d_TimeStep(2)=MIN(d_TimeStep(2),MINVAL(d_DFLScale(:))*4./d_Lambda_max(2,iElem))
#endif
#if FV_ENABLED == 2
  END IF
#endif
END DO

TimeStep(:) = d_TimeStep

#if USE_MPI
TimeStep(3)=-errType ! reduce with timestep, minus due to MPI_MIN
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStep,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_FLEXI,iError)
errType=INT(-TimeStep(3))
#endif /*USE_MPI*/
ViscousTimeStep=(TimeStep(2) .LT. TimeStep(1))
CalcTimeStep=MINVAL(TimeStep(1:2))
END FUNCTION CALCTIMESTEP


!==================================================================================================================================
!> Compute the maximum convective/viscous Eigenvalue per element
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE CalcMaxEigenvalue(Nloc,nElems,U,EOS_Vars,&
                                                     Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,&
                                                     MetricsAdv,&
#if PARABOLIC
                                                     MetricsVisc,&
#endif
                                                     sJ,&
                                                     MaxEigenvalue)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)          :: Nloc
INTEGER,VALUE,INTENT(IN)          :: nElems
REAL,DEVICE,INTENT(IN)            :: U(PP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)
REAL,DEVICE,INTENT(IN)            :: Metrics_fTilde(3,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)
REAL,DEVICE,INTENT(IN)            :: Metrics_gTilde(3,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)
REAL,DEVICE,INTENT(IN)            :: Metrics_hTilde(3,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)
REAL,DEVICE,INTENT(IN)            :: MetricsAdv(    3,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)
#if PARABOLIC
REAL,DEVICE,INTENT(IN)            :: MetricsVisc(   3,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)
#endif
REAL,DEVICE,INTENT(IN)            :: sJ(              0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)
REAL,DEVICE,INTENT(IN)            :: EOS_Vars(PP_nVarEOS)
REAL,DEVICE,INTENT(OUT)           :: MaxEigenvalue(2,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k,ElemID
INTEGER                        :: threadID,rest
REAL                           :: UE(PP_2Var)
REAL                           :: vsJ(3)
REAL                           :: c
#if PARABOLIC
REAL                           :: mu
REAL,SHARED                    :: Max_Lambda(6,0:Nloc,0:ZDIM(Nloc))
REAL                           :: Max_LambdaTmp(6)
#else
REAL,SHARED                    :: Max_Lambda(3,0:Nloc,0:ZDIM(Nloc))
REAL                           :: Max_LambdaTmp(3)
#endif
!==================================================================================================================================
! Get thread indices
threadID = (blockidx%x-1) * blockdim%x + threadidx%x
! Get ElemID of current thread
ElemID   =        (threadID-1)/((Nloc+1)**2)+1 ! Elems are 1-indexed
rest     = threadID-(ElemID-1)*((Nloc+1)**2)
! Get jk indices of current thread
k        = (rest-1)/(Nloc+1)!**1
rest     =  rest- k*(Nloc+1)!**1
j        = (rest-1)!/(Nloc+1)**0
rest     =  rest- j!*(Nloc+1)**0

IF (ElemID.GT.nElems) RETURN

! Each Thread gets whole i-slice to compute to save shared memory space
Max_LambdaTmp=0.
DO i=0,Nloc
  ! TODO: ATTENTION: Temperature of UE not filled!!!
  UE(EXT_CONS)=U(:,i,j,k,ElemID)
  UE(EXT_SRHO)=1./UE(EXT_DENS)
  UE(EXT_VELV)=VELOCITY_HE(UE)
  UE(EXT_PRES)=PRESSURE_UE_EOS(UE,EOS_Vars)
  UE(EXT_TEMP)=TEMPERATURE_UE_EOS(UE,EOS_Vars)
  c=SPEEDOFSOUND_UE_EOS(UE,EOS_Vars)
  vsJ=UE(EXT_VELV)*sJ(i,j,k,ElemID)
  Max_LambdaTmp(  1)=MAX(Max_LambdaTmp(1),ABS(SUM(Metrics_fTilde(:,i,j,k,ElemID)*vsJ)) &
                                          + c*MetricsAdv(1,i,j,k,ElemID))
  Max_LambdaTmp(  2)=MAX(Max_LambdaTmp(2),ABS(SUM(Metrics_gTilde(:,i,j,k,ElemID)*vsJ)) &
                                          + c*MetricsAdv(2,i,j,k,ElemID))
  Max_LambdaTmp(  3)=MAX(Max_LambdaTmp(3),ABS(SUM(Metrics_hTilde(:,i,j,k,ElemID)*vsJ)) &
                                          + c*MetricsAdv(3,i,j,k,ElemID))
#if PARABOLIC
  mu=VISCOSITY_PRIM_EOS(UE(EXT_PRIM,i,j,k,ElemID),EOS_Vars)
  Max_LambdaTmp(4:6)=MAX(Max_LambdaTmp(4:6),mu*UE(EXT_SRHO)*MetricsVisc(:,i,j,k,ElemID))
#endif
END DO ! i

! Write maximum of own i-slice to shared memory
Max_Lambda(:,j,k) = Max_LambdaTmp(:)

! Ensure that every Thread of current element has written to shared memory
CALL SYNCTHREADS()

! Only Root thread computes the maximum
IF ((j.GT.0).OR.(k.GT.0)) RETURN

! Compute maximum convective eigenvalue (still has to be scaled by CFL)
MaxEigenvalue(1,ElemID)= MAXVAL(Max_Lambda(1,:,:)) &
                       + MAXVAL(Max_Lambda(2,:,:)) &
                       + MAXVAL(Max_Lambda(3,:,:))

! Compute maximum viscous eigenvalue (still has to be scaled by CFL)
#if PARABOLIC
MaxEigenvalue(2,ElemID)= MAXVAL(Max_Lambda(4,:,:)) &
                       + MAXVAL(Max_Lambda(5,:,:)) &
                       + MAXVAL(Max_Lambda(6,:,:))
#endif

END SUBROUTINE CalcMaxEigenvalue

!==================================================================================================================================
!> Deallocate CalcTimeStep arrays
!==================================================================================================================================
SUBROUTINE FinalizeCalctimestep()
! MODULES
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(d_MetricsAdv)
#if PARABOLIC
SDEALLOCATE(d_MetricsVisc)
#endif
END SUBROUTINE

END MODULE MOD_CalcTimeStep
