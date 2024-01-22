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
#include "eos.h"

!==================================================================================================================================
!> Contains the global variables needed by the ideal gas equation of state.
!==================================================================================================================================
MODULE MOD_EOS_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
#if PARABOLIC
REAL              :: mu0               !< dynamic viscosity \f$\mu\f$
REAL              :: Pr                !< Prandtl number
REAL              :: KappasPr          !< \f$\kappa\f$/Pr
REAL              :: lambda            !< thermal conductivity
#if PP_VISC==1
REAL              :: Ts                !< Sutherland temperature
REAL              :: cSuth             !< Parameters used in muSuth
#endif
#if (PP_VISC==1) || (PP_VISC==2)
REAL              :: Tref,ExpoSuth     !< Parameters used in muSuth and power law
#endif
#endif /*PARABOLIC*/
REAL              :: cp                !< specific heat at constant pressure
REAL              :: cv                !< specific heat at constant volume
REAL              :: Kappa             !< heat capacity ratio / isentropic exponent
REAL              :: KappaM1           !< = \f$\kappa - 1\f$
REAL              :: sKappaM1          !< = \f$1/(\kappa -1)\f$
REAL              :: KappaP1           !< = \f$\kappa + 1\f$
REAL              :: sKappaP1          !< = \f$1/(\kappa +1)\f$
REAL              :: R                 !< specific gas constant
REAL              ::   EOS_Vars(PP_nVarEOS) !< Contains all EOS variables
REAL,DEVICE       :: d_EOS_Vars(PP_nVarEOS) !< Contains all EOS variables
LOGICAL           :: EosInitIsDone=.FALSE.
!==================================================================================================================================
END MODULE MOD_EOS_Vars
