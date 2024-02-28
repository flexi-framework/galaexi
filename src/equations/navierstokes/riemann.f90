!=================================================================================================================================
! Copyright (c) 2010-2017 Prof. Claus-Dieter Munz
! Copyright (c) 2016-2017 Gregor Gassner (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Florian Hindenlang (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Andrew Winters (github.com/project-fluxo/fluxo)
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
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  PPURE SUBROUTINE RiemannInt(Kappa,F_L,F_R,U_LL,U_RR,F)
    REAL,INTENT(IN)                    :: Kappa
    REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
    REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(RiemannInt),POINTER :: Riemann_pointer    !< pointer defining the standard inner Riemann solver

INTEGER,PARAMETER      :: PRM_RIEMANN_LF            = 1
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLC          = 2
INTEGER,PARAMETER      :: PRM_RIEMANN_ROE           = 3
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEL2         = 32
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEENTROPYFIX = 33
#ifndef SPLIT_DG
INTEGER,PARAMETER      :: PRM_RIEMANN_HLL           = 4
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLE          = 5
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLEM         = 6
#else
INTEGER,PARAMETER      :: PRM_RIEMANN_CH            = 7
INTEGER,PARAMETER      :: PRM_RIEMANN_Average       = 0
#endif

INTERFACE InitRiemann
  MODULE PROCEDURE InitRiemann
END INTERFACE

INTERFACE Riemann
  MODULE PROCEDURE Riemann
  MODULE PROCEDURE Riemann_Point
  MODULE PROCEDURE Riemann_Point_Device
  MODULE PROCEDURE Riemann_Side
  MODULE PROCEDURE Riemann_Sides
  !MODULE PROCEDURE Riemann_Point_CUDA
  MODULE PROCEDURE Riemann_Side_CUDA
  MODULE PROCEDURE Riemann_Sides_CUDA
END INTERFACE

INTERFACE Riemann_Solver
#if   RIEMANN==0
  MODULE PROCEDURE Riemann_LF
#elif RIEMANN==1
  MODULE PROCEDURE Riemann_Roe
#elif RIEMANN==2
  MODULE PROCEDURE Riemann_RoeL2
#elif RIEMANN==3
  MODULE PROCEDURE Riemann_RoeEntropyFix
#elif RIEMANN==4
  MODULE PROCEDURE Riemann_HLL
#elif RIEMANN==5
  MODULE PROCEDURE Riemann_HLLC
#elif RIEMANN==6
  MODULE PROCEDURE Riemann_HLLE
#elif RIEMANN==7
  MODULE PROCEDURE Riemann_HLLEM
#endif /*RIEMANN*/
END INTERFACE

#if PARABOLIC
INTERFACE ViscousFlux
  MODULE PROCEDURE ViscousFlux
  !MODULE PROCEDURE ViscousFlux_Point
  MODULE PROCEDURE ViscousFlux_Side
  MODULE PROCEDURE ViscousFlux_Sides_CUDA
END INTERFACE
#endif

INTERFACE FinalizeRiemann
  MODULE PROCEDURE FinalizeRiemann
END INTERFACE


PUBLIC::DefineParametersRiemann
PUBLIC::InitRiemann
PUBLIC::Riemann
PUBLIC::FinalizeRiemann
#if PARABOLIC
PUBLIC::ViscousFlux
#endif
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! CALL prms%SetSection("Riemann")
! CALL prms%CreateIntFromStringOption('Riemann',   "Riemann solver to be used: LF, HLLC, Roe, RoeEntropyFix, HLL, HLLE, HLLEM", &
!                                                  "RoeEntropyFix")
! CALL addStrListEntry('Riemann','lf',           PRM_RIEMANN_LF)
! CALL addStrListEntry('Riemann','hllc',         PRM_RIEMANN_HLLC)
! CALL addStrListEntry('Riemann','roe',          PRM_RIEMANN_ROE)
! CALL addStrListEntry('Riemann','roeentropyfix',PRM_RIEMANN_ROEENTROPYFIX)
! CALL addStrListEntry('Riemann','roel2',        PRM_RIEMANN_ROEL2)
! #ifndef SPLIT_DG
! CALL addStrListEntry('Riemann','hll',          PRM_RIEMANN_HLL)
! CALL addStrListEntry('Riemann','hlle',         PRM_RIEMANN_HLLE)
! CALL addStrListEntry('Riemann','hllem',        PRM_RIEMANN_HLLEM)
! #else
! CALL addStrListEntry('Riemann','ch',           PRM_RIEMANN_CH)
! CALL addStrListEntry('Riemann','avg',          PRM_RIEMANN_Average)
! #endif
END SUBROUTINE DefineParametersRiemann

!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: Riemann
!==================================================================================================================================
!Riemann = GETINTFROMSTR('Riemann')
!SELECT CASE(Riemann)
!CASE(PRM_RIEMANN_LF)
!  Riemann_pointer => Riemann_LF
!CASE(PRM_RIEMANN_HLLC)
!  Riemann_pointer => Riemann_HLLC
!CASE(PRM_RIEMANN_ROE)
!  Riemann_pointer => Riemann_Roe
!CASE(PRM_RIEMANN_ROEENTROPYFIX)
!  Riemann_pointer => Riemann_RoeEntropyFix
!CASE(PRM_RIEMANN_ROEL2)
!  Riemann_pointer => Riemann_RoeL2
!#ifndef SPLIT_DG
!CASE(PRM_RIEMANN_HLL)
!  Riemann_pointer => Riemann_HLL
!CASE(PRM_RIEMANN_HLLE)
!  Riemann_pointer => Riemann_HLLE
!CASE(PRM_RIEMANN_HLLEM)
!  Riemann_pointer => Riemann_HLLEM
!#else
!CASE(PRM_RIEMANN_CH)
!  Riemann_pointer => Riemann_CH
!CASE(PRM_RIEMANN_Average)
!  Riemann_pointer => Riemann_FluxAverage
!#endif /*SPLIT_DG*/
!CASE DEFAULT
!  CALL CollectiveStop(__STAMP__,&
!    'Riemann solver not defined!')
!END SELECT
SWRITE(*, '(A)') " | GALAEXI currently only supports selection of the Riemann solver during compilation."
SWRITE(*, '(A)') " | Any setting found in the parameter file will be ignored."
#if   RIEMANN==0
  SWRITE(*,'(A,A2)')' | Riemann is set to Lax-Friedrichs (LF)'
#elif RIEMANN==1
  SWRITE(*,'(A,A2)')' | Riemann is set to Roe'
#elif RIEMANN==2
  SWRITE(*,'(A,A2)')' | Riemann is set to RoeL2'
#elif RIEMANN==3
  SWRITE(*,'(A,A2)')' | Riemann is set to RoeEntropyFix'
#elif RIEMANN==4
  SWRITE(*,'(A,A2)')' | Riemann is set to HLL'
#elif RIEMANN==5
  SWRITE(*,'(A,A2)')' | Riemann is set to HLLC'
#elif RIEMANN==6
  SWRITE(*,'(A,A2)')' | Riemann is set to HLLE'
#elif RIEMANN==7
  SWRITE(*,'(A,A2)')' | Riemann is set to HLLEM'
#endif /*RIEMANN*/
END SUBROUTINE InitRiemann

!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann(FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,Kappa)
! MODULES
USE MOD_Flux, ONLY: EvalEulerFlux1D_fast
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)                         :: Kappa      !< ratio of specific heats
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_L        !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_R        !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_L    !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_R    !< primitive solution at right side of the interface
REAL,DIMENSION(3          ),INTENT(IN)  :: nv,t1,t2   !< normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: FOut       !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR
!==================================================================================================================================
! Momentum has to be rotated using the normal system individual for each

! left state: U_L
U_LL(EXT_DENS)=U_L(DENS)
U_LL(EXT_SRHO)=1./U_LL(EXT_DENS)
U_LL(EXT_ENER)=U_L(ENER)
U_LL(EXT_PRES)=UPrim_L(PRES)
! rotate velocity in normal and tangential direction
U_LL(EXT_VEL1)=DOT_PRODUCT(UPrim_L(VELV),nv(:))
U_LL(EXT_VEL2)=DOT_PRODUCT(UPrim_L(VELV),t1(:))
U_LL(EXT_MOM1)=U_LL(EXT_DENS)*U_LL(EXT_VEL1)
U_LL(EXT_MOM2)=U_LL(EXT_DENS)*U_LL(EXT_VEL2)
#if PP_dim==3
U_LL(EXT_VEL3)=DOT_PRODUCT(UPrim_L(VELV),t2(:))
U_LL(EXT_MOM3)=U_LL(EXT_DENS)*U_LL(EXT_VEL3)
#else
U_LL(EXT_VEL3)=0.
U_LL(EXT_MOM3)=0.
#endif

! right state: U_R
U_RR(EXT_DENS)=U_R(DENS)
U_RR(EXT_SRHO)=1./U_RR(EXT_DENS)
U_RR(EXT_ENER)=U_R(ENER)
U_RR(EXT_PRES)=UPrim_R(PRES)
! rotate momentum in normal and tangential direction
U_RR(EXT_VEL1)=DOT_PRODUCT(UPRIM_R(VELV),nv(:))
U_RR(EXT_VEL2)=DOT_PRODUCT(UPRIM_R(VELV),t1(:))
U_RR(EXT_MOM1)=U_RR(EXT_DENS)*U_RR(EXT_VEL1)
U_RR(EXT_MOM2)=U_RR(EXT_DENS)*U_RR(EXT_VEL2)
#if PP_dim==3
U_RR(EXT_VEL3)=DOT_PRODUCT(UPRIM_R(VELV),t2(:))
U_RR(EXT_MOM3)=U_RR(EXT_DENS)*U_RR(EXT_VEL3)
#else
U_RR(EXT_VEL3)=0.
U_RR(EXT_MOM3)=0.
#endif

# ifndef SPLIT_DG
CALL EvalEulerFlux1D_fast(U_LL,F_L)
CALL EvalEulerFlux1D_fast(U_RR,F_R)
#endif /*SPLIT_DG*/

CALL Riemann_Solver(F,F_L,F_R,U_LL,U_RR,Kappa)

! Back rotate the normal flux into Cartesian direction
Fout(DENS)=F(DENS)
Fout(MOMV)=nv(:)*F(MOM1)  &
          +t1(:)*F(MOM2)  &
#if PP_dim==3
          +t2(:)*F(MOM3)
#else
          +0.
#endif
Fout(ENER)=F(ENER)
END SUBROUTINE Riemann

!==================================================================================================================================
!> Computes the numerical flux for a single point calling the flux calculation.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
PPURE SUBROUTINE Riemann_Point(FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2)
! MODULES
USE MOD_EOS_Vars, ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_L       !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_R       !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_L   !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_R   !< primitive solution at right side of the interface
REAL,DIMENSION(3          ),INTENT(IN)  :: nv,t1,t2  !> normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: FOut      !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: p,q
!==================================================================================================================================
CALL Riemann(FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,Kappa)
END SUBROUTINE Riemann_Point

!==================================================================================================================================
!> Computes the numerical flux for a single point calling the flux calculation.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE) SUBROUTINE Riemann_Point_Device(FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,EOS_Vars)
! MODULES
USE MOD_EOS_Vars, ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_L       !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_R       !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_L   !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_R   !< primitive solution at right side of the interface
REAL,DIMENSION(3          ),INTENT(IN)  :: nv,t1,t2  !> normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVarEOS ),INTENT(IN)  :: EOS_Vars  !< EOS-spcific variables
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: FOut      !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: p,q
!==================================================================================================================================
CALL Riemann(FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,EOS_Vars(EOS_KAPPA))
END SUBROUTINE Riemann_Point_Device

!==================================================================================================================================
!> Computes the numerical flux for a side calling the flux calculation pointwise.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
PPURE SUBROUTINE Riemann_Side(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2)
! MODULES
USE MOD_EOS_Vars, ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                          :: Nloc      !< local polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L       !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R       !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L   !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R   !< primitive solution at right side of the interface
REAL,DIMENSION(3          ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2  !> normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: FOut      !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  CALL Riemann(FOut(:,p,q),U_L(:,p,q),U_R(:,p,q),UPrim_L(:,p,q),UPrim_R(:,p,q),nv(:,p,q),t1(:,p,q),t2(:,p,q),Kappa)
END DO; END DO
END SUBROUTINE Riemann_Side

!==================================================================================================================================
!> Computes the numerical flux for a side calling the flux calculation pointwise.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
PPURE SUBROUTINE Riemann_Sides(Nloc,nSides,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2)
! MODULES
USE MOD_EOS_Vars, ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: Nloc      !< local polynomial degree
INTEGER,INTENT(IN)        :: nSides    !< number of sides in arrays
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: U_L      !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: U_R      !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: UPrim_L  !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: UPrim_R  !< primitive solution at right side of the interface
REAL,DIMENSION(3          ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: nv,t1,t2 !> normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(OUT) :: FOut     !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: p,q,iSide
!==================================================================================================================================
DO iSide=1,nSides
  DO q=0,ZDIM(Nloc); DO p=0,Nloc
    CALL Riemann(FOut(:,p,q,iSide),U_L(:,p,q,iSide),U_R(:,p,q,iSide),UPrim_L(:,p,q,iSide),UPrim_R(:,p,q,iSide),&
                   nv(:,p,q,iSide), t1(:,p,q,iSide), t2(:,p,q,iSide),Kappa)
  END DO; END DO !p,q
END DO !iSide
END SUBROUTINE Riemann_Sides

!==================================================================================================================================
!> Kernel to start the 1D evaluation of the Rimann routine
!==================================================================================================================================
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE Riemann_CUDA_Kernel(nDOF,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,Kappa)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)  :: nDOF   !< number of degress of freedom
REAL,VALUE,INTENT(IN)     :: Kappa  !< ratio of specific heats
REAL,DIMENSION(PP_nVar    ,nDOF),INTENT(IN)  :: U_L       !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,nDOF),INTENT(IN)  :: U_R       !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,nDOF),INTENT(IN)  :: UPrim_L   !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,nDOF),INTENT(IN)  :: UPrim_R   !< primitive solution at right side of the interface
REAL,DIMENSION(3          ,nDOF),INTENT(IN)  :: nv,t1,t2  !> normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ,nDOF),INTENT(OUT) :: FOut      !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
!==================================================================================================================================
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF(i.LE.nDOF) CALL Riemann(FOut(:,i),U_L(:,i),U_R(:,i),UPrim_L(:,i),UPrim_R(:,i),nv(:,i),t1(:,i),t2(:,i),Kappa)
END SUBROUTINE Riemann_CUDA_Kernel

!==================================================================================================================================
!> Computes the numerical flux for a side calling the flux calculation pointwise.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
PPURE SUBROUTINE Riemann_Side_CUDA(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2)
! MODULES
USE MOD_EOS_Vars, ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                                 :: Nloc     !< local polynomial degree
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L      !< conservative solution at left side of the interface
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R      !< conservative solution at right side of the interface
REAL,DEVICE,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L  !< primitive solution at left side of the interface
REAL,DEVICE,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R  !< primitive solution at right side of the interface
REAL,DEVICE,DIMENSION(3          ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2 !> normal vector and tangential vectors at side
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: FOut     !< advective flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
!==================================================================================================================================
nDOF=(Nloc+1)*(ZDIM(Nloc)+1)
CALL Riemann_CUDA_Kernel<<<nDOF/nThreads+1,nThreads>>>(nDOF,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,Kappa)
END SUBROUTINE Riemann_Side_CUDA

!==================================================================================================================================
!> Computes the numerical flux for a side calling the flux calculation pointwise.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
PPURE SUBROUTINE Riemann_Sides_CUDA(Nloc,nSides,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,streamID)
! MODULES
USE CUDAFOR
USE MOD_GPU      ,ONLY:DefaultStream
USE MOD_EOS_Vars, ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: Nloc      !< local polynomial degree
INTEGER,INTENT(IN)        :: nSides    !< number of sides in arrays
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: U_L      !< conservative solution at left side of the interface
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: U_R      !< conservative solution at right side of the interface
REAL,DEVICE,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: UPrim_L  !< primitive solution at left side of the interface
REAL,DEVICE,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: UPrim_R  !< primitive solution at right side of the interface
REAL,DEVICE,DIMENSION(3          ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(IN)  :: nv,t1,t2 !> normal vector and tangential vectors at side
REAL,DEVICE,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc),nSides),INTENT(OUT) :: FOut     !< advective flux
INTEGER(KIND=CUDA_STREAM_KIND),OPTIONAL,INTENT(IN) :: streamID
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
INTEGER(KIND=CUDA_STREAM_KIND) :: mystream
!==================================================================================================================================
mystream=DefaultStream
IF (PRESENT(streamID)) mystream=streamID

nDOF=(Nloc+1)*(ZDIM(Nloc)+1)*nSides
CALL Riemann_CUDA_Kernel<<<nDOF/nThreads+1,nThreads,0,mystream>>>(nDOF,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,Kappa)
END SUBROUTINE Riemann_Sides_CUDA

#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE ViscousFlux(F,UPrim_L,UPrim_R, &
                             gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
#if EDDYVISCOSITY
                            ,muSGS_L,muSGS_R &
#endif
                            ,EOS_Vars &
                            )
! MODULES
USE MOD_Flux,ONLY: EvalDiffFlux3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                           !> solution in primitive variables at left/right side of the interface
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim_L,UPrim_R
                                           !> solution gradients in x/y/z-direction left/right of the interface
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DIMENSION(3             ),INTENT(IN)  :: nv  !< normal vector
REAL,DIMENSION(PP_nVarEOS    ),INTENT(IN)  :: EOS_Vars   !< EOS-specific variables
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: F   !< viscous flux
#if EDDYVISCOSITY
REAL,INTENT(IN)                            :: muSGS_L,muSGS_R    !> eddy viscosity left/right of the interface
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)  :: diffFluxX_L,diffFluxY_L,diffFluxZ_L
REAL,DIMENSION(PP_nVar)  :: diffFluxX_R,diffFluxY_R,diffFluxZ_R
REAL                     :: mu,lambda
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
mu=VISCOSITY_PRIM_EOS(UPrim_L(:,i),EOS_Vars)
lambda=THERMAL_CONDUCTIVITY_EOS(mu,EOS_Vars)
CALL EvalDiffFlux3D(UPrim_L,   gradUx_L,   gradUy_L,   gradUz_L  &
                           ,diffFluxX_L,diffFluxY_L,diffFluxZ_L  &
#if EDDYVISCOSITY
                   ,muSGS_L &
#endif
                   ,mu,lambda&
      )
mu=VISCOSITY_PRIM_EOS(UPrim_R(:,i),EOS_Vars)
lambda=THERMAL_CONDUCTIVITY_EOS(mu,EOS_Vars)
CALL EvalDiffFlux3D(UPrim_R,   gradUx_R,   gradUy_R,   gradUz_R  &
                           ,diffFluxX_R,diffFluxY_R,diffFluxZ_R  &
#if EDDYVISCOSITY
                   ,muSGS_R&
#endif
                   ,mu,lambda&
      )
! Arithmetic mean of the fluxes
F(:)=0.5*(nv(1)*(diffFluxX_L(:)+diffFluxX_R(:)) &
         +nv(2)*(diffFluxY_L(:)+diffFluxY_R(:)) &
         +nv(3)*(diffFluxZ_L(:)+diffFluxZ_R(:)))
END SUBROUTINE ViscousFlux

!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux_Side(Nloc,F,UPrim_L,UPrim_R, &
                            gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
#if EDDYVISCOSITY
                           ,muSGS_L,muSGS_R &
#endif
                           )
! MODULES
USE MOD_EOS_Vars,ONLY:EOS_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                             :: Nloc     !< local polynomial degree
                                                               !> solution in primitive variables at left/right side of interface
REAL,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L,UPrim_R
                                                               !> solution gradients in x/y/z-direction left/right of interface
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DIMENSION(3             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv  !< normal vector
REAL,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: F   !< viscous flux
#if EDDYVISCOSITY
REAL,DIMENSION(1             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: muSGS_L,muSGS_R   !> eddy viscosity left/right of the interface
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  CALL ViscousFlux(F(:,p,q),UPrim_L(:,p,q),UPrim_R(:,p,q), &
                   gradUx_L(:,p,q),gradUy_L(:,p,q),gradUz_L(:,p,q),gradUx_R(:,p,q),gradUy_R(:,p,q),gradUz_R(:,p,q),nv(:,p,q) &
#if EDDYVISCOSITY
                   ,muSGS_L(:,p,q),muSGS_R(:,p,q) &
#endif
                   ,EOS_Vars &
                   )
END DO; END DO
END SUBROUTINE ViscousFlux_Side

!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux and add to input flux array
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux_Sides_CUDA(Nloc,nSides,F,UPrim_L,UPrim_R, &
                                  gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
                                 ,streamID)
! MODULES
USE CUDAFOR
USE MOD_GPU     ,ONLY:DefaultStream
USE MOD_Flux    ,ONLY:EvalDiffFlux3D
USE MOD_EOS_Vars,ONLY:d_EOS_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                  :: Nloc     !< local polynomial degree
INTEGER,INTENT(IN)                  :: nSides   !< number of sides in array
                                                                               !> solution in primitive variables at left/right side of interface
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc),  nSides),INTENT(IN)    :: UPrim_L,UPrim_R
                                                                                 !> solution gradients in x/y/z-direction left/right of interface
REAL,DEVICE,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc),  nSides),INTENT(IN)    :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DEVICE,DIMENSION(3             ,0:Nloc,0:ZDIM(Nloc),  nSides),INTENT(IN)    :: nv  !< normal vector
REAL,DEVICE,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc),  nSides),INTENT(INOUT) :: F   !< viscous flux
INTEGER(KIND=CUDA_STREAM_KIND),OPTIONAL,INTENT(IN) :: streamID
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: nDOF
INTEGER,PARAMETER :: nThreads=256
INTEGER(KIND=CUDA_STREAM_KIND) :: mystream
!==================================================================================================================================
mystream=DefaultStream
IF (PRESENT(streamID)) mystream=streamID

! Don't forget the diffusion contribution, my young padawan
nDOF = (Nloc+1)*(ZDIM(Nloc)+1)*nSides

! Compute NSE Diffusion flux
CALL ViscousFlux_Kernel_CUDA<<<nDOF/nThreads+1,nThreads,0,mystream>>>(nDOF,F,UPrim_L,UPrim_R, &
                             gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
                             ,d_EOS_Vars)
END SUBROUTINE ViscousFlux_Sides_CUDA

!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux and add to input flux array
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
PPURE ATTRIBUTES(GLOBAL) SUBROUTINE ViscousFlux_Kernel_CUDA(nDOF,F,UPrim_L,UPrim_R, &
                                    gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
                                    ,EOS_Vars)
! MODULES
USE MOD_Flux, ONLY: EvalDiffFlux3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)  :: nDOF     !< local polynomial degree
                                                         !> solution in primitive variables at left/right side of interface
REAL,DEVICE,DIMENSION(PP_nVarPrim   ,nDOF),INTENT(IN)    :: UPrim_L,UPrim_R
                                                         !> solution gradients in x/y/z-direction left/right of interface
REAL,DEVICE,DIMENSION(PP_nVarLifting,nDOF),INTENT(IN)    :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DEVICE,DIMENSION(3             ,nDOF),INTENT(IN)    :: nv  !< normal vector
REAL,DEVICE,DIMENSION(PP_nVar       ,nDOF),INTENT(INOUT) :: F   !< viscous flux
REAL,DEVICE,DIMENSION(PP_nVarEOS),INTENT(IN) :: EOS_Vars !< EOS-specific variables
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
REAL    :: mu,lambda
REAL,DEVICE,DIMENSION(PP_nVar) :: normalDiffFlux
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF (i.LE.nDOF) THEN
  ! Left material parameters
  mu=VISCOSITY_PRIM_EOS(UPrim_L(:,i),EOS_Vars)
  lambda=THERMAL_CONDUCTIVITY_EOS(mu,EOS_Vars)
  ! Left flux
  CALL EvalDiffFlux3D( UPrim_L(:,i), &
                      gradUx_L(:,i), &
                      gradUy_L(:,i), &
                      gradUz_L(:,i), &
                normalDiffFlux(:  ), &
                            nv(:,i), &
                          mu,lambda)
  ! Arithmetic mean of the fluxes (Add directly to avoid additional temporary variables)
  F(:,i) = F(:,i) + 0.5*normalDiffFlux ! F = F+0.5*(Flux_L+Flux_R)

  ! Right material parameters
  mu=VISCOSITY_PRIM_EOS(UPrim_L(:,i),EOS_Vars)
  lambda=THERMAL_CONDUCTIVITY_EOS(mu,EOS_Vars)
  ! Right flux
  CALL EvalDiffFlux3D( UPrim_R(:,i), &
                      gradUx_R(:,i), &
                      gradUy_R(:,i), &
                      gradUz_R(:,i), &
                normalDiffFlux(:  ), &
                            nv(:,i), &
                          mu,lambda)
  ! Arithmetic mean of the fluxes (Add directly to avoid additional temporary variables)
  F(:,i) = F(:,i) + 0.5*normalDiffFlux ! F = F+0.5*(Flux_L+Flux_R)
END IF
END SUBROUTINE ViscousFlux_Kernel_CUDA

#endif /* PARABOLIC */

!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_LF(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: LambdaMax
!==================================================================================================================================
! Lax-Friedrichs
LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )
#ifndef SPLIT_DG
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(CONS) - U_LL(CONS)))
#else
! get split flux
CALL SplitSurfaceFlux(U_LL,U_RR,F)
! compute surface flux
F = F - 0.5*LambdaMax*(U_RR(CONS) - U_LL(CONS))
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_LF

!=================================================================================================================================
!> Harten-Lax-Van-Leer Riemann solver resolving contact discontinuity
!=================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_HLLC(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL    :: RoeVel(3),RoeH,Roec,absVel
REAL    :: Ssl,Ssr,SStar
REAL    :: U_Star(PP_nVar),EStar
REAL    :: sMu_L,sMu_R
!REAL    :: c_L,c_R
!=================================================================================================================================
! HLLC flux

! Version A: Basic Davis estimate for wave speed
!Ssl = U_LL(EXT_VEL1) - SPEEDOFSOUND_HE(U_LL)
!Ssr = U_RR(EXT_VEL1) + SPEEDOFSOUND_HE(U_RR)

! Version B: Basic Davis estimate for wave speed
!c_L = SPEEDOFSOUND_HE(U_LL)
!c_R = SPEEDOFSOUND_HE(U_RR)
!Ssl = MIN(U_LL(EXT_VEL1) - c_L,U_RR(EXT_VEL1) - c_R)
!Ssr = MAX(U_LL(EXT_VEL1) + c_L,U_RR(EXT_VEL1) + c_R)

! Version C: Better Roe estimate for wave speeds Davis, Einfeldt
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L       )     * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = SQRT((Kappa-1.)*(RoeH-0.5*absVel))
Ssl       = RoeVel(1) - Roec
Ssr       = RoeVel(1) + Roec

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  sMu_L = Ssl - U_LL(EXT_VEL1)
  sMu_R = Ssr - U_RR(EXT_VEL1)
  SStar = (U_RR(EXT_PRES) - U_LL(EXT_PRES) + U_LL(EXT_MOM1)*sMu_L - U_RR(EXT_MOM1)*sMu_R) / (U_LL(EXT_DENS)*sMu_L - U_RR(EXT_DENS)*sMu_R)
  IF ((Ssl .LE. 0.).AND.(SStar .GE. 0.)) THEN
    EStar  = TOTALENERGY_HE(U_LL) + (SStar-U_LL(EXT_VEL1))*(SStar + U_LL(EXT_PRES)*U_LL(EXT_SRHO)/sMu_L)
    U_Star = U_LL(EXT_DENS) * sMu_L/(Ssl-SStar) * (/ 1., SStar, U_LL(EXT_VEL2:EXT_VEL3), EStar /)
    F=F_L+Ssl*(U_Star-U_LL(CONS))
  ELSE
    EStar  = TOTALENERGY_HE(U_RR) + (SStar-U_RR(EXT_VEL1))*(SStar + U_RR(EXT_PRES)*U_RR(EXT_SRHO)/sMu_R)
    U_Star = U_RR(EXT_DENS) * sMu_R/(Ssr-SStar) * (/ 1., SStar, U_RR(EXT_VEL2:EXT_VEL3), EStar /)
    F=F_R+Ssr*(U_Star-U_RR(CONS))
  END IF
END IF ! subsonic case
END SUBROUTINE Riemann_HLLC

!=================================================================================================================================
!> Roe's approximate Riemann solver
!=================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_Roe(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
#endif /*SPLIT_DG*/
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL,DIMENSION(PP_nVar) :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(PP_nVar+1)
!=================================================================================================================================
! Roe flux
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(CONS)     = U_RR(CONS) - U_LL(CONS)
Delta_U(DELTA_U6) = Delta_U(DELTA_U5)-(Delta_U(DELTA_U3)-RoeVel(DELTA_U2)*Delta_U(DELTA_U1))*RoeVel(2) -&
                    (Delta_U(DELTA_U4)-RoeVel(DELTA_U3)*Delta_U(DELTA_U1))*RoeVel(DELTA_U3)
! calculate factors
Alpha3 = Delta_U(DELTA_U3) - RoeVel(DELTA_U2)*Delta_U(DELTA_U1)
Alpha4 = Delta_U(DELTA_U4) - RoeVel(DELTA_U3)*Delta_U(DELTA_U1)
Alpha2 = ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)
Alpha1 = 0.5/Roec * (Delta_U(DELTA_U1)*(RoeVel(1)+Roec) - Delta_U(DELTA_U2) - Roec*Alpha2)
Alpha5 = Delta_U(DELTA_U1) - Alpha1 - Alpha2
#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R) - &
       Alpha1*ABS(a(1))*r1 - &
       Alpha2*ABS(a(2))*r2 - &
       Alpha3*ABS(a(3))*r3 - &
       Alpha4*ABS(a(4))*r4 - &
       Alpha5*ABS(a(5))*r5)
#else
! get split flux
CALL SplitSurfaceFlux(U_LL,U_RR,F)
! assemble Roe flux
F = F - 0.5*(Alpha1*ABS(a(1))*r1 + &
             Alpha2*ABS(a(2))*r2 + &
             Alpha3*ABS(a(3))*r3 + &
             Alpha4*ABS(a(4))*r4 + &
             Alpha5*ABS(a(5))*r5)
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_Roe

!=================================================================================================================================
!> Roe's approximate Riemann solver using the Harten and Hymen II entropy fix, see
!> Pelanti, Marica & Quartapelle, Luigi & Vigevano, L & Vigevano, Luigi. (2018):
!>  A review of entropy fixes as applied to Roe's linearization.
!=================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_RoeEntropyFix(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iVar
REAL                    :: c_L,c_R
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                    :: RoeVel(3),RoeH,Roec,RoeDens
REAL,DIMENSION(PP_nVar) :: r1,r2,r3,r4,r5,a,al,ar,Delta_U,Alpha  ! Roe eigenvectors
REAL                    :: tmp,da
!=================================================================================================================================
c_L       = SPEEDOFSOUND_HE(U_LL)
c_R       = SPEEDOFSOUND_HE(U_RR)
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
RoeDens   = SQRT(U_LL(EXT_DENS)*U_RR(EXT_DENS))
! Roe+Pike version of Roe Riemann solver

! calculate jump
Delta_U(DELTA_U1)   = U_RR(EXT_DENS) - U_LL(EXT_DENS)
Delta_U(DELTA_UV)   = U_RR(EXT_VELV) - U_LL(EXT_VELV)
Delta_U(DELTA_U5)   = U_RR(EXT_PRES) - U_LL(EXT_PRES)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate wave strenghts
tmp      = 0.5/(Roec*Roec)
Alpha(1) = tmp*(Delta_U(DELTA_U5)-RoeDens*Roec*Delta_U(DELTA_U2))
Alpha(2) = Delta_U(DELTA_U1) - Delta_U(DELTA_U5)*2.*tmp
Alpha(3) = RoeDens*Delta_U(DELTA_U3)
Alpha(4) = RoeDens*Delta_U(DELTA_U4)
Alpha(5) = tmp*(Delta_U(DELTA_U5)+RoeDens*Roec*Delta_U(DELTA_U2))

! Harten+Hyman entropy fix (apply only for acoustic waves, don't fix r)

al(1) = U_LL(EXT_VEL1) - c_L
al(2) = U_LL(EXT_VEL1)
al(3) = U_LL(EXT_VEL1)
al(4) = U_LL(EXT_VEL1)
al(5) = U_LL(EXT_VEL1) + c_L
ar(1) = U_RR(EXT_VEL1) - c_R
ar(2) = U_RR(EXT_VEL1)
ar(3) = U_RR(EXT_VEL1)
ar(4) = U_RR(EXT_VEL1)
ar(5) = U_RR(EXT_VEL1) + c_R
! HH1
!IF(ABS(a(1)).LT.da1) a(1)=da1
!IF(ABS(a(5)).LT.da5) a(5)=da5
! HH2
DO iVar=1,5
  da = MAX(0.,a(iVar)-al(iVar),ar(iVar)-a(iVar))

  IF(ABS(a(iVar)).LT.da) THEN
    a(iVar)=0.5*(a(iVar)*a(iVar)/da+da)
  ELSE
    a(iVar) = ABS(a(iVar))
  END IF
END DO

#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R)        - &
       Alpha(1)*a(1)*r1 - &
       Alpha(2)*a(2)*r2 - &
       Alpha(3)*a(3)*r3 - &
       Alpha(4)*a(4)*r4 - &
       Alpha(5)*a(5)*r5)
#else
! get split flux
CALL SplitSurfaceFlux(U_LL,U_RR,F)
! for KG or PI flux eigenvalues have to be altered to ensure consistent KE dissipation
! assemble Roe flux
F= F - 0.5*(Alpha(1)*a(1)*r1 + &
            Alpha(2)*a(2)*r2 + &
            Alpha(3)*a(3)*r3 + &
            Alpha(4)*a(4)*r4 + &
            Alpha(5)*a(5)*r5)
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_RoeEntropyFix

!=================================================================================================================================
!> low mach number Roe's approximate Riemann solver according to Oßwald(2015)
!=================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_RoeL2(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL                    :: Ma_loc ! local Mach-Number
REAL,DIMENSION(PP_nVar) :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(PP_nVar+1)
!=================================================================================================================================
! Roe flux
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(CONS) = U_RR(EXT_CONS) - U_LL(EXT_CONS)
Delta_U(DELTA_U6)   = Delta_U(DELTA_U5)-(Delta_U(DELTA_U3)-RoeVel(2)*Delta_U(DELTA_U1))*RoeVel(2) - &
                      (Delta_U(DELTA_U4)-RoeVel(3)*Delta_U(DELTA_U1))*RoeVel(3)

! low Mach-Number fix
Ma_loc = SQRT(absVel)/(Roec*SQRT(kappa))
Delta_U(DELTA_UV) = Delta_U(DELTA_UV) * Ma_loc

! calculate factors
Alpha3 = Delta_U(DELTA_U3) - RoeVel(2)*Delta_U(DELTA_U1)
Alpha4 = Delta_U(DELTA_U4) - RoeVel(3)*Delta_U(DELTA_U1)
Alpha2 = ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)
Alpha1 = 0.5/Roec * (Delta_U(DELTA_U1)*(RoeVel(1)+Roec) - Delta_U(DELTA_U2) - Roec*Alpha2)
Alpha5 = Delta_U(DELTA_U1) - Alpha1 - Alpha2

#ifndef SPLIT_DG
! assemble Roe flux
F=0.5*((F_L+F_R) - &
       Alpha1*ABS(a(1))*r1 - &
       Alpha2*ABS(a(2))*r2 - &
       Alpha3*ABS(a(3))*r3 - &
       Alpha4*ABS(a(4))*r4 - &
       Alpha5*ABS(a(5))*r5)
#else
! get split flux
CALL SplitSurfaceFlux(U_LL,U_RR,F)
! assemble Roe flux
F = F - 0.5*(Alpha1*ABS(a(1))*r1 + &
             Alpha2*ABS(a(2))*r2 + &
             Alpha3*ABS(a(3))*r3 + &
             Alpha4*ABS(a(4))*r4 + &
             Alpha5*ABS(a(5))*r5)
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_RoeL2

!=================================================================================================================================
!> Standard Harten-Lax-Van-Leer Riemann solver without contact discontinuity
!=================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_HLL(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roec
REAL    :: Ssl,Ssr
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
! HLL flux
! Basic Davis estimate for wave speed
!Ssl = U_LL(EXT_VEL1) - c_L
!Ssr = U_RR(EXT_VEL1) + c_R
! Better Roe estimate for wave speeds Davis, Einfeldt
Ssl = RoeVel(1) - Roec
Ssr = RoeVel(1) + Roec
! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(EXT_CONS)-U_LL(EXT_CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLL

!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt Riemann solver
!=================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_HLLE(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roec
REAL    :: Ssl,Ssr,beta
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
! HLLE flux (positively conservative)
beta=BETA_RIEMANN_H()
SsL=MIN(RoeVel(1)-Roec,U_LL(EXT_VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(EXT_VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(EXT_CONS)-U_LL(EXT_CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLE

!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt-Munz Riemann solver
!=================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_HLLEM(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                   :: H_L,H_R
REAL                                   :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                                   :: RoeVel(3),RoeH,Roec,RoeDens
REAL                                   :: Ssl,Ssr
REAL                                   :: Alpha(2:4),delta,beta
REAL,DIMENSION(PP_nVar)                :: r2,r3,r4  ! Roe eigenvectors + jump in prims
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
RoeDens   = SQRT(U_LL(EXT_DENS)*U_RR(EXT_DENS))
! HLLEM flux (positively conservative)
beta=BETA_RIEMANN_H()
SsL=MIN(RoeVel(1)-Roec,U_LL(EXT_VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(EXT_VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  ! delta
  delta = Roec/(Roec+ABS(0.5*(Ssl+Ssr)))

  ! mean eigenvectors
  Alpha(2)   = (U_RR(EXT_DENS)-U_LL(EXT_DENS))  - (U_RR(EXT_PRES)-U_LL(EXT_PRES))/(Roec*Roec)
  Alpha(3:4) = RoeDens*(U_RR(EXT_VEL2:EXT_VEL3) - U_LL(EXT_VEL2:EXT_VEL3))
  r2 = (/ 1., RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel /)
  r3 = (/ 0., 0.,        1.,        0.,        RoeVel(2)  /)
  r4 = (/ 0., 0.,        0.,        1.,        RoeVel(3)  /)

  F=(Ssr*F_L-Ssl*F_R + Ssl*Ssr* &
     (U_RR(EXT_CONS)-U_LL(EXT_CONS) - delta*(r2*Alpha(2)+r3*Alpha(3)+r4*Alpha(4))))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLEM

#ifdef SPLIT_DG
!==================================================================================================================================
!> Riemann solver using purely the average fluxes
!==================================================================================================================================
PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_FluxAverage(F,F_L,F_R,U_LL,U_RR,Kappa)
! MODULES
USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! get split flux
CALL SplitSurfaceFlux(U_LL,U_RR,F)
END SUBROUTINE Riemann_FluxAverage

!!==================================================================================================================================
!!> kinetic energy preserving and entropy consistent flux according to Chandrashekar (2012)
!!==================================================================================================================================
!PPURE ATTRIBUTES(DEVICE,HOST) SUBROUTINE Riemann_CH(F,F_L,F_R,U_LL,U_RR,Kappa)
!! MODULES
!USE MOD_SplitFlux     ,ONLY: SplitSurfaceFlux
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT / OUTPUT VARIABLES
!REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR !> extended solution vector on the left/right side of the interface
!REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R   !> advection fluxes on the left/right side of the interface
!REAL,INTENT(IN)                    :: Kappa     !> ratio of specific heats
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL                               :: LambdaMax
!REAL                               :: beta_LL,beta_RR   ! auxiliary variables for the inverse Temperature
!REAL                               :: rhoMean           ! auxiliary variable for the mean density
!REAL                               :: uMean,vMean,wMean ! auxiliary variable for the average velocities
!REAL                               :: betaLogMean       ! auxiliary variable for the logarithmic mean inverse temperature
!!==================================================================================================================================
!! Lax-Friedrichs
!LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )
!
!! average quantities
!rhoMean = 0.5*(U_LL(EXT_DENS) + U_RR(EXT_DENS))
!uMean   = 0.5*(U_LL(EXT_VEL1) + U_RR(EXT_VEL1))
!vMean   = 0.5*(U_LL(EXT_VEL2) + U_RR(EXT_VEL2))
!wMean   = 0.5*(U_LL(EXT_VEL3) + U_RR(EXT_VEL3))
!
!! inverse temperature
!beta_LL = 0.5*U_LL(EXT_DENS)/U_LL(EXT_PRES)
!beta_RR = 0.5*U_RR(EXT_DENS)/U_RR(EXT_PRES)
!
!! logarithmic mean
!CALL GetLogMean(beta_LL,beta_RR,betaLogMean)
!
!! get split flux
!CALL SplitSurfaceFlux(U_LL,U_RR,F)
!
!!compute flux
!F(DENS:MOM3) = F(DENS:MOM3) - 0.5*LambdaMax*(U_RR(EXT_DENS:EXT_MOM3)-U_LL(EXT_DENS:EXT_MOM3))
!F(ENER)      = F(ENER)      - 0.5*LambdaMax*( &
!         (U_RR(EXT_DENS)-U_LL(EXT_DENS))*(0.5/(Kappa-1.)/betaLogMean +0.5*(U_RR(EXT_VEL1)*U_LL(EXT_VEL1)+U_RR(EXT_VEL2)*U_LL(EXT_VEL2)+U_RR(EXT_VEL3)*U_LL(EXT_VEL3))) &
!         +rhoMean*uMean*(U_RR(EXT_VEL1)-U_LL(EXT_VEL1)) + rhoMean*vMean*(U_RR(EXT_VEL2)-U_LL(EXT_VEL2)) + rhoMean*wMean*(U_RR(EXT_VEL3)-U_LL(EXT_VEL3)) &
!         +0.5*rhoMean/(Kappa-1.)*(1./beta_RR - 1./beta_LL))
!
!END SUBROUTINE Riemann_CH
#endif /*SPLIT_DG*/

!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann


END MODULE MOD_Riemann
