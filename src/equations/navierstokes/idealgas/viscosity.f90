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
!> Wrapper module related to viscosity related routines. Automatically provides correct variables depending on type of viscosity
!> to be used.
!==================================================================================================================================
MODULE MOD_Viscosity
#if PARABOLIC
! MODULES
#if   PP_VISC == 0
USE MOD_EOS_Vars,     ONLY:mu0
#elif PP_VISC == 1
USE MOD_EOS_Vars,     ONLY:R
#elif PP_VISC == 2
USE MOD_EOS_Vars,     ONLY:mu0,R,ExpoSuth
#endif

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------

#if PP_VISC == 1
INTERFACE muSuth
  MODULE PROCEDURE muSuth
  MODULE PROCEDURE muSuth_CUDA
END INTERFACE

PUBLIC::muSuth
#endif

!==================================================================================================================================

CONTAINS


#if PP_VISC == 1
!==================================================================================================================================
!> Sutherland's formula can be used to derive the dynamic viscosity of an ideal gas as a function of the temperature
!>
!> Initialization of mu0, Ts, Tref, ExpoSuth and cSuth takes place in SUBROUTINE IniEquation
!>
!> Temperatures above the Sutherlands Temperature Ts are computed according to (1)
!> 1) T >= Ts:    mu = mu0 * (T/Tref)^(expo) *  (Tref+TS)/(T+TS)
!>
!> below Ts a linear dependence is assumed, (2)
!> 2) T < Ts:    mu = mu0*T/Tref*c
!>
!> with c = (Ts/Tref)^exp*(1+(Ts/Tref))/(2(Ts/Tref)²) for steady transition from (1) to (2) at T = Ts.
!>
!> This is only valid for Temperatures in the range 0 < T < 555 K
!> For more detail see White, F. M.
!> ATTENTION!!!!! The global variable Tref=1./Tref and Ts=Ts/Tref !!!!!
!==================================================================================================================================
ELEMENTAL FUNCTION muSuth(T)
! MODULES
USE MOD_EOS_Vars, ONLY: mu0,Tref,Ts,ExpoSuth,cSuth
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                :: T !< Temperature
REAL                           :: muSuth
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: TnoDim
!==================================================================================================================================
TnoDim=T*Tref ! Tref=1./Tref !
IF(TnoDim .GE. Ts)THEN  ! Attention: only valid for T < 550K. But we don't know what to do for higher temperatures...
  muSuth=mu0*TnoDim**ExpoSuth*(1+Ts)/(TnoDim+Ts)  ! Ts=Ts/Tref !
ELSE
  muSuth=mu0*TnoDim*cSuth
END IF
END FUNCTION muSuth

!==================================================================================================================================
!> Sutherland's formula can be used to derive the dynamic viscosity of an ideal gas as a function of the temperature
!>
!> Initialization of mu0, Ts, Tref, ExpoSuth and cSuth takes place in SUBROUTINE IniEquation
!>
!> Temperatures above the Sutherlands Temperature Ts are computed according to (1)
!> 1) T >= Ts:    mu = mu0 * (T/Tref)^(expo) *  (Tref+TS)/(T+TS)
!>
!> below Ts a linear dependence is assumed, (2)
!> 2) T < Ts:    mu = mu0*T/Tref*c
!>
!> with c = (Ts/Tref)^exp*(1+(Ts/Tref))/(2(Ts/Tref)²) for steady transition from (1) to (2) at T = Ts.
!>
!> This is only valid for Temperatures in the range 0 < T < 555 K
!> For more detail see White, F. M.
!> ATTENTION!!!!! The global variable Tref=1./Tref and Ts=Ts/Tref !!!!!
!==================================================================================================================================
PPURE ATTRIBUTES(HOST,DEVICE) FUNCTION muSuth_CUDA(T,EOS_Vars)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)        :: T !< Temperature
REAL,DEVICE,INTENT(IN) :: EOS_Vars(PP_nVarEOS) !< EOS-specific variables
REAL                   :: muSuth_CUDA
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: TnoDim
!==================================================================================================================================
TnoDim=T*EOS_Vars(EOS_TREF) ! Tref=1./Tref !
IF(TnoDim .GE. EOS_Vars(EOS_TS))THEN  ! Attention: only valid for T < 550K. But we don't know what to do for higher temperatures...
  muSuth_CUDA=EOS_Vars(EOS_MU0)*TnoDim**EOS_Vars(EOS_EXPOSUTH)*(1+EOS_Vars(EOS_TS))/(TnoDim+EOS_Vars(EOS_TS))  ! Ts=Ts/Tref !
ELSE
  muSuth_CUDA=EOS_Vars(EOS_MU0)*TnoDim*EOS_Vars(EOS_CSUTH)
END IF
END FUNCTION muSuth_CUDA
#endif /*PP_VISC == 1*/

#endif /*PARABOLIC*/
END MODULE MOD_Viscosity
