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

!==================================================================================================================================
!>\brief Computes the DGSEM volume integral
!> The volume integral is computed via the weak form of the DG method
!> Computes the volume integral contribution based on U and updates Ut
!> Volume integral is split into integral of advection and diffusion part
!==================================================================================================================================
#include "flexi.h"
MODULE MOD_VolInt
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE VolInt
#ifndef SPLIT_DG
  MODULE PROCEDURE VolInt_weakForm
#else
  MODULE PROCEDURE VolInt_splitForm
#endif /*SPLIT_DG*/
END INTERFACE

#if  PARABOLIC && !VOLINT_VISC
INTERFACE VolInt_Visc
  MODULE PROCEDURE VolInt_weakForm_Visc
END INTERFACE
#endif


PUBLIC::VolInt
#if  PARABOLIC && !VOLINT_VISC
PUBLIC::VolInt_Visc
#endif
!==================================================================================================================================
CONTAINS

#if PARABOLIC && !VOLINT_VISC
!==================================================================================================================================
!> Computes the viscous part volume integral of the weak DG form according to Kopriva
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is NOT overwritten, but instead added to the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolInt_weakForm_Visc(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: D_hat_T,nDOFElem,UPrim
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_Flux         ,ONLY: EvalDiffFlux3D  ! computes volume fluxes in local coordinates
USE MOD_Lifting_Vars ,ONLY: gradUx,gradUy,gradUz
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Time derivative of the volume integral (viscous part)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ) :: f,g,h     !< Volume viscous fluxes at GP
!==================================================================================================================================
! Diffusive part
DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.1) CYCLE ! FV Elem
#endif
  CALL EvalDiffFlux3D( UPrim(:,:,:,:,iElem),&
                      gradUx(:,:,:,:,iElem),&
                      gradUy(:,:,:,:,iElem),&
                      gradUz(:,:,:,:,iElem),&
                      f,g,h,iElem)

  CALL VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem,0),&
                                     Metrics_gTilde(:,:,:,:,iElem,0),&
                                     Metrics_hTilde(:,:,:,:,iElem,0))

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      ! Update the time derivative with the spatial derivatives of the transformed fluxes
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_Hat_T(l,i)*f(:,l,j,k) + &
#if PP_dim==3
                                              D_Hat_T(l,k)*h(:,i,j,l) + &
#endif
                                              D_Hat_T(l,j)*g(:,i,l,k)
    END DO ! l
  END DO; END DO; END DO !i,j,k
END DO ! iElem
END SUBROUTINE VolInt_weakForm_Visc
#endif /*PARABOLIC && !VOLINT_VISC*/


!==================================================================================================================================
!> Computes the advection and viscous part volume integral of the weak DG form according to Kopriva
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!==================================================================================================================================
#ifndef SPLIT_DG
SUBROUTINE VolInt_weakForm(d_Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars         ,ONLY: nDOFElem,d_f,d_g,d_h,nElems_Block_volInt
USE MOD_DG_Vars         ,ONLY: d_U,d_UPrim,d_D_Hat_T
USE MOD_Mesh_Vars       ,ONLY: d_Metrics_fTilde,d_Metrics_gTilde,d_Metrics_hTilde,nElems
USE MOD_Flux            ,ONLY: EvalFlux3D      ! computes volume fluxes in local coordinates
USE MOD_Flux            ,ONLY: EvalTransformedFlux3D
USE MOD_ApplyDMatrixCons,ONLY: ApplyDMatrixCons
#if VOLINT_VISC
USE MOD_Flux            ,ONLY: EvalDiffFlux3D  ! computes volume fluxes in local coordinates
USE MOD_Lifting_Vars    ,ONLY: d_gradUx,d_gradUy,d_gradUz
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DEVICE,INTENT(OUT)   :: d_Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Time derivative of the volume integral (viscous part)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem,lastElem,nElems_myBlock
INTEGER            :: nDOF
INTEGER,PARAMETER  :: nThreads=256
!==================================================================================================================================
! Diffusive part
DO iElem=1,nElems,nElems_Block_volInt
  lastElem    = MIN(nElems,iElem+nElems_Block_volInt-1)
  nElems_myBlock = lastElem-iElem+1
  nDOF = nDOFElem*nElems_myBlock

  CALL EvalTransformedFlux3D(nDOF &
                            ,d_U(    :,:,:,:,iElem:lastElem) &
                            ,d_UPrim(:,:,:,:,iElem:lastElem) &
#if PARABOLIC
                            ,d_gradUx(:,:,:,:,iElem:lastElem) &
                            ,d_gradUy(:,:,:,:,iElem:lastElem) &
                            ,d_gradUz(:,:,:,:,iElem:lastElem) &
#endif
                            ,d_f( :,:,:,:,1:nElems_myBlock)  &
                            ,d_g( :,:,:,:,1:nElems_myBlock)  &
                            ,d_h( :,:,:,:,1:nElems_myBlock)  &
                            ,d_Metrics_fTilde(:,:,:,:,iElem:lastElem,0) &
                            ,d_Metrics_gTilde(:,:,:,:,iElem:lastElem,0) &
                            ,d_Metrics_hTilde(:,:,:,:,iElem:lastElem,0))

  CALL ApplyDMatrixCons(nElems_myBlock &
                       ,d_Ut(:,:,:,:,iElem:lastElem) &
                       ,d_f( :,:,:,:,1:nElems_myBlock) &
                       ,d_g( :,:,:,:,1:nElems_myBlock) &
                       ,d_h( :,:,:,:,1:nElems_myBlock) &
                       ,d_D_Hat_T &
                       )
END DO ! iElem
END SUBROUTINE VolInt_weakForm
#endif

#ifdef SPLIT_DG
!==================================================================================================================================
!> Computes the advection and viscous part volume integral in SplitDG formulation
!>
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut is overwritten with the volume flux derivatives
!> Attention 3: the factor of 2 in front of the derivative matrix entries is incorporated into the split fluxes!
!> Attention 4: This is the strong form of the DGSEM! Substracting the inner flux is incorporated into the used D matrix, which
!>              saves performance but only works for Gauss-Lobatto points. So no changes in the surface integral or fill flux
!>              routines are necessary. For Gauss-Lobatto points, these inner fluxes cancel exactly with the non-zero (consistent)
!>              fluxes at the outer points [-1,1] of the volume integral, which is why these inner fluxes never have to be actuallys
!>              computed.
!> Attention 5: For Gauss-Lobatto points, the matrix DVolSurf will always be equal to 0 on the main diagonal for the inner points
!>              and becomes 0 for the outer points at [-1,1] after considering the inner fluxes of the strong form (Attention 4).
!>              Thus, the (consistent) fluxes will always be multiplied by zero and we don't have to take them into account at all.
!>
!> For details on the derivation see Gassner, Gregor J., Andrew R. Winters, and David A. Kopriva.
!> "Split form nodal discontinuous Galerkin schemes with summation-by-parts property for the compressible Euler equations."
!> Journal of Computational Physics 327 (2016): 39-66.
!==================================================================================================================================
SUBROUTINE VolInt_splitForm(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: DVolSurf,UPrim,U
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde,nElems
#if PP_dim==3 || VOLINT_VISC
USE MOD_Mesh_Vars    ,ONLY: Metrics_hTilde
#endif
USE MOD_SplitFlux    ,ONLY:SplitDGVolume_pointer ! computes volume fluxes in split formulation
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Time derivative of the volume integral (viscous part)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,l,iElem
REAL,DIMENSION(PP_nVar                     )  :: Flux         !< temp variable for split flux
!==================================================================================================================================
DO iElem=1,nElems

! For split DG, the matrix DVolSurf will always be equal to 0 on the main diagonal. Thus, the (consistent) fluxes will always be
! multiplied by zero and we don't have to take them into account at all.
  ! We need to nullify the Ut array
  Ut(:,:,:,:,iElem) = 0.


  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO l=i+1,PP_N
       ! compute split flux in x-direction
       CALL SplitDGVolume_pointer(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem), &
                                  U(:,l,j,k,iElem),UPrim(:,l,j,k,iElem), &
                                  Metrics_fTilde(:,i,j,k,iElem,0),Metrics_fTilde(:,l,j,k,iElem,0),Flux)
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,i)*Flux(:)
       !symmetry
       Ut(:,l,j,k,iElem) = Ut(:,l,j,k,iElem) + DVolSurf(i,l)*Flux(:)
    END DO ! l

    DO l=j+1,PP_N
       ! compute split flux in y-direction
       CALL SplitDGVolume_pointer(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem), &
                                  U(:,i,l,k,iElem),UPrim(:,i,l,k,iElem), &
                                  Metrics_gTilde(:,i,j,k,iElem,0),Metrics_gTilde(:,i,l,k,iElem,0),Flux)
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,j)*Flux(:)
       !symmetry
       Ut(:,i,l,k,iElem) = Ut(:,i,l,k,iElem) + DVolSurf(j,l)*Flux(:)
    END DO ! l

#if PP_dim==3
    DO l=k+1,PP_N
       ! compute split flux in z-direction
       CALL SplitDGVolume_pointer(U(:,i,j,k,iElem),UPrim(:,i,j,k,iElem), &
                                  U(:,i,j,l,iElem),UPrim(:,i,j,l,iElem), &
                                  Metrics_hTilde(:,i,j,k,iElem,0),Metrics_hTilde(:,i,j,l,iElem,0),Flux)
       ! add up time derivative
       Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + DVolSurf(l,k)*Flux(:)
       !symmetry
       Ut(:,i,j,l,iElem) = Ut(:,i,j,l,iElem) + DVolSurf(k,l)*Flux(:)
    END DO ! l
#endif /*PP_dim==3*/

  END DO; END DO; END DO !i,j,k
END DO ! iElem
END SUBROUTINE VolInt_splitForm
#endif /*SPLIT_DG*/

!==================================================================================================================================
!> Compute the tranformed states for all conservative variables using the metric terms
!==================================================================================================================================
PPURE SUBROUTINE VolInt_Metrics(nDOFs,f,g,h,Mf,Mg,Mh)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                          :: nDOFs    !< Number of DOFs per element
                                                        !> Metrics terms
REAL,DIMENSION(3,nDOFs),INTENT(IN)          :: Mf,Mg,Mh
                                                        !> Volume fluxes at GP to be transformed from physical to reference space
REAL,DIMENSION(PP_nVar,nDOFs),INTENT(INOUT) :: f,g,h
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                     :: i
REAL,DIMENSION(PP_nVar)                     :: fTilde,gTilde !< Auxiliary variables needed to store the fluxes at one GP
#if PP_dim==3
REAL,DIMENSION(PP_nVar)                     :: hTilde !< Auxiliary variables needed to store the fluxes at one GP
#endif
!==================================================================================================================================
DO i=1,nDOFs
  fTilde=f(:,i)
  gTilde=g(:,i)
#if PP_dim==3
  hTilde=h(:,i)

  ! Compute the transformed fluxes with the metric terms
  ! Attention 1: we store the transformed fluxes in f,g,h again
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i) + &
           hTilde*Mf(3,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i) + &
           hTilde*Mg(3,i)
  h(:,i) = fTilde*Mh(1,i) + &
           gTilde*Mh(2,i) + &
           hTilde*Mh(3,i)
#else
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i)
#endif
END DO ! i
END SUBROUTINE VolInt_Metrics

!==================================================================================================================================
!> Compute the tranformed states for all conservative variables using the metric terms
!==================================================================================================================================
ATTRIBUTES(GLOBAL) SUBROUTINE VolInt_Metrics_GPU(nDOFs,f,g,h,Mf,Mg,Mh)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,VALUE,INTENT(IN)                    :: nDOFs    !< Number of DOFs per element
                                                        !> Metrics terms
REAL,DIMENSION(3,nDOFs),INTENT(IN)          :: Mf,Mg,Mh
                                                        !> Volume fluxes at GP to be transformed from physical to reference space
REAL,DIMENSION(1:PP_nVar,nDOFs),INTENT(INOUT) :: f,g,h
!@cuf ATTRIBUTES(DEVICE) Mf,Mg,Mh,f,g,h
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                     :: i
REAL,DIMENSION(PP_nVar),DEVICE              :: fTilde,gTilde !< Auxiliary variables needed to store the fluxes at one GP
#if PP_dim==3
REAL,DIMENSION(PP_nVar),DEVICE              :: hTilde !< Auxiliary variables needed to store the fluxes at one GP
#endif
!==================================================================================================================================
i = (blockidx%x-1) * blockdim%x + threadidx%x
IF (i.LE.nDOFs) THEN
  fTilde=f(:,i)
  gTilde=g(:,i)
#if PP_dim==3
  hTilde=h(:,i)

  ! Compute the transformed fluxes with the metric terms
  ! Attention 1: we store the transformed fluxes in f,g,h again
  f(:,i) = fTilde(:)*Mf(1,i) + &
           gTilde(:)*Mf(2,i) + &
           hTilde(:)*Mf(3,i)
  g(:,i) = fTilde(:)*Mg(1,i) + &
           gTilde(:)*Mg(2,i) + &
           hTilde(:)*Mg(3,i)
  h(:,i) = fTilde(:)*Mh(1,i) + &
           gTilde(:)*Mh(2,i) + &
           hTilde(:)*Mh(3,i)
#else
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i)
#endif
ENDIF
END SUBROUTINE VolInt_Metrics_GPU


END MODULE MOD_VolInt
