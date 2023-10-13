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
!> \brief Computes the DGSEM spatial operator and updates residual Ut

!> Contains the routines to
!> - initialize and finalize the DG global variables and the DG basis
!> - compute the DG spatial operators/residuals(Ut) using U from the volume, surface and source contribution, incl.
!> lifting for the gradients and parallelization
!==================================================================================================================================
MODULE MOD_DG
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
INTERFACE FillIni
  MODULE PROCEDURE FillIni
END INTERFACE


! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitDG
  MODULE PROCEDURE InitDG
END INTERFACE


INTERFACE DGTimeDerivative_weakForm
  MODULE PROCEDURE DGTimeDerivative_weakForm
END INTERFACE


INTERFACE FinalizeDG
  MODULE PROCEDURE FinalizeDG
END INTERFACE


PUBLIC::InitDG,DGTimeDerivative_weakForm,FinalizeDG
!==================================================================================================================================



CONTAINS

!==================================================================================================================================
!> Allocate all global DG variables like U (solution in volume), U_slave/U_master (solution on faces), Flux, Ut (DG time derivative),
!> also fill the initial solution and call init DG basis. Operator building are also initialized by calling InitDGBasis.
!==================================================================================================================================
SUBROUTINE InitDG()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars
USE MOD_Interpolation_Vars,   ONLY: xGP,wGP,L_minus,L_plus
USE MOD_Interpolation_Vars,   ONLY: InterpolationInitIsDone
USE MOD_Restart_Vars,         ONLY: DoRestart,RestartInitIsDone
USE MOD_Mesh_Vars,            ONLY: nElems,nSides,Elem_xGP,MeshInitIsDone
USE MOD_ChangeBasisByDim,     ONLY: ChangeBasisVolume
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Check if all the necessary initialization is done before
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.DGInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    'InitDG not ready to be called or already called.')
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT DG...'

! Pre-compute the dg operator building blocks (differentiation matrices and prolongation operators)
CALL InitDGBasis(PP_N, xGP,wGP,L_minus,L_plus,D ,D_T ,D_Hat ,D_Hat_T ,L_HatMinus ,L_HatPlus)

! Copy LHat Vectors to GPU
!@cuf ALLOCATE(d_D_T(0:PP_N,0:PP_N))
!@cuf ALLOCATE(d_D_Hat_T(0:PP_N,0:PP_N))
!@cuf ALLOCATE(d_L_HatMinus(0:PP_N),d_L_HatPlus(0:PP_N))
!@cuf d_L_HatPlus  = L_HatPlus
!@cuf d_L_HatMinus = L_HatMinus
!@cuf d_D_Hat_T    = D_Hat_T
!@cuf d_D_T        = D_T

! Allocate the local DG solution (JU or U): element-based
ALLOCATE(U(        PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
!@cuf ALLOCATE(d_U(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
! Allocate the time derivative / solution update /residual vector dU/dt: element-based
ALLOCATE(Ut(       PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
!@cuf ALLOCATE(d_Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
U =0.
Ut=0.
d_U =U
d_Ut=Ut

! Allocate the 2D solution vectors on the sides, one array for the data belonging to the proc (the master)
! and one for the sides which belong to another proc (slaves): side-based
ALLOCATE(U_master(PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(U_slave( PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
!@cuf ALLOCATE(d_U_master(PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
!@cuf ALLOCATE(d_U_slave( PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
U_master=0.
U_slave =0.
d_U_master=U_master
d_U_slave =U_slave

! Repeat the U, U_Minus, U_Plus structure for the primitive quantities
ALLOCATE(UPrim(       PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
!@cuf ALLOCATE(d_UPrim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
!@cuf ALLOCATE(d_UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
!@cuf ALLOCATE(d_UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides))
UPrim=0.
UPrim_master=0.
UPrim_slave=0.
d_UPrim=UPrim
d_UPrim_master=UPrim_master
d_UPrim_slave =UPrim_slave

! Allocate the UPrim_boundary for the boundary fluxes
ALLOCATE(UPrim_boundary(PP_nVarPrim,0:PP_N,0:PP_NZ))

! Allocate two fluxes per side (necessary for coupling of FV and DG)
ALLOCATE(Flux_master(PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(Flux_slave (PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
!@cuf ALLOCATE(d_Flux_master(PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
!@cuf ALLOCATE(d_Flux_slave( PP_nVar,0:PP_N,0:PP_NZ,1:nSides))
Flux_master=0.
Flux_slave =0.
d_Flux_master=Flux_master
d_Flux_slave =Flux_slave

nElems_Block_volInt=nElems
ALLOCATE(d_f (PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems_Block_volInt))
ALLOCATE(d_g (PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems_Block_volInt))
ALLOCATE(d_h (PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems_Block_volInt))


! variables for performance tricks
nDOFFace=(PP_N+1)**(PP_dim-1)
nDOFElem=(PP_N+1)**PP_dim
nTotalU=PP_nVar*nDOFElem*nElems

! Fill the solution vector U with the initial solution by interpolation, if not filled through restart already
IF(.NOT.DoRestart)THEN
  CALL FillIni(PP_N,Elem_xGP,U)
END IF

DGInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT DG DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitDG



!==================================================================================================================================
!> Allocate and initialize the building blocks for the DG operator: Differentiation matrices and prolongation operators
!==================================================================================================================================
SUBROUTINE InitDGbasis(N_in,xGP,wGP,L_Minus,L_Plus,D,D_T,D_Hat,D_Hat_T,L_HatMinus,L_HatPlus)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Interpolation,    ONLY: GetNodesAndWeights
USE MOD_Basis,            ONLY: PolynomialDerivativeMatrix,LagrangeInterpolationPolys,PolynomialMassMatrix
#ifdef SPLIT_DG
USE MOD_DG_Vars,          ONLY: DVolSurf ! Transpose of differentiation matrix used for calculating the strong form
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                             :: N_in                   !< Polynomial degree
REAL,DIMENSION(0:N_in),INTENT(IN)              :: xGP                    !< Gauss/Gauss-Lobatto Nodes
REAL,DIMENSION(0:N_in),INTENT(IN)              :: wGP                    !< Gauss/Gauss-Lobatto Weights
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Minus                !< Values of lagrange polynomials at \f$ \xi = -1 \f$
REAL,DIMENSION(0:N_in),INTENT(IN)              :: L_Plus                 !< Values of lagrange polynomials at \f$ \xi = +1 \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D                      !< Differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_T                    !< Transpose of differentation matrix
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat                  !< Differentiation matrix premultiplied by mass matrix,
                                                                         !< \f$ \hat{D} = M^{-1} D^T M \f$
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)    :: D_Hat_T                !< Transpose of D_Hat matrix \f$ \hat{D}^T \f$
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatMinus             !< Values of lagrange polynomials at \f$ \xi = -1 \f$
                                                                         !< premultiplied with mass matrix
REAL,ALLOCATABLE,DIMENSION(:)  ,INTENT(OUT)    :: L_HatPlus              !< Values of lagrange polynomials at \f$ \xi = +1 \f$
                                                                         !< premultiplied with mass matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:N_in,0:N_in)              :: M,Minv
!==================================================================================================================================

ALLOCATE(L_HatMinus(0:N_in), L_HatPlus(0:N_in))
ALLOCATE(D(    0:N_in,0:N_in), D_T(    0:N_in,0:N_in))
ALLOCATE(D_Hat(0:N_in,0:N_in), D_Hat_T(0:N_in,0:N_in))
! Compute Differentiation matrix D for given Gausspoints
CALL PolynomialDerivativeMatrix(N_in,xGP,D)
D_T=TRANSPOSE(D)

!Build Mass Matrix
CALL PolynomialMassMatrix(N_in,xGP,wGP,M,Minv)

! Build D_Hat matrix. D^ = - (M^(-1) * D^T * M)
D_Hat  = -MATMUL(Minv,MATMUL(TRANSPOSE(D),M))
D_Hat_T= TRANSPOSE(D_Hat)

#ifdef SPLIT_DG
! Use a modified D matrix for the strong form volume integral, that incorporates the inner fluxes that are subtracted from the
! surfaces
ALLOCATE(DVolSurf(0:N_in,0:N_in))
DVolSurf = D_T
! Modify the D matrix here, the integral over the inner fluxes at the boundaries will then be automatically done in the volume
! integral. The factor 1/2 is needed since we incorporate a factor of 2 in the split fluxes themselves!
! For Gauss-Lobatto points, these inner flux contributions cancel exactly with entries in the DVolSurf matrix, resulting in zeros.
DVolSurf(   0,   0) = DVolSurf(   0   ,0) + 1.0/(2.0 * wGP(   0))  ! = 0. (for LGL)
DVolSurf(N_in,N_in) = DVolSurf(N_in,N_in) - 1.0/(2.0 * wGP(N_in))  ! = 0. (for LGL)
#endif /*SPLIT_DG*/

! interpolate to left and right face (1 and -1 in reference space) and pre-divide by mass matrix
L_HatPlus  = MATMUL(Minv,L_Plus)
L_HatMinus = MATMUL(Minv,L_Minus)
END SUBROUTINE InitDGbasis



!==================================================================================================================================
!> \brief Computes the residual Ut = \f$ \frac {d\vec{U}} {dt} \f$ from the current solution U employing the DG method.
!> Computes the weak DGSEM space operator from surface, volume and source contributions. To do this we need to:
!> - Prolong the solution from the volume to the interface
!> - Invoke the lifting operator to calculate the gradients
!> - Perform the volume integral
!> - Perform the surface integral
!> - If needed, add source terms to the residual
!==================================================================================================================================
SUBROUTINE DGTimeDerivative_weakForm(t)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_DG_Vars             ,ONLY: Ut,U,U_slave,U_master,Flux_master,Flux_slave,L_HatMinus,L_HatPlus,nTotalU
!@cuf USE MOD_DG_Vars          ,ONLY: d_L_HatPlus,d_L_HatMinus
USE MOD_DG_Vars             ,ONLY: UPrim,UPrim_master,UPrim_slave,nDOFElem,nDOFFace
!@cuf USE MOD_DG_Vars          ,ONLY: d_U,d_UPrim,d_Ut
!@cuf USE MOD_DG_Vars          ,ONLY: d_U_master,d_U_slave,d_UPrim_master,d_UPrim_Slave
!@cuf USE MOD_DG_Vars          ,ONLY: d_Flux_master,d_Flux_slave,d_UPrim_master,d_UPrim_Slave
USE MOD_VolInt
USE MOD_SurfIntCons         ,ONLY: SurfIntCons
USE MOD_ProlongToFaceCons   ,ONLY: ProlongToFaceCons
USE MOD_FillFlux            ,ONLY: FillFlux
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
!@cuf USE MOD_Interpolation_Vars  ,ONLY: d_L_Minus,d_L_Plus
USE MOD_Overintegration_Vars,ONLY: OverintegrationType
USE MOD_Overintegration,     ONLY: Overintegration
USE MOD_ChangeBasisByDim    ,ONLY: ChangeBasisVolume
USE MOD_TestCase            ,ONLY: TestcaseSource
USE MOD_TestCase_Vars       ,ONLY: doTCSource
USE MOD_Equation            ,ONLY: GetPrimitiveStateSurface,GetConservativeStateSurface
USE MOD_EOS                 ,ONLY: ConsToPrim
USE MOD_Exactfunc           ,ONLY: CalcSource
USE MOD_Equation_Vars       ,ONLY: doCalcSource
USE MOD_Sponge              ,ONLY: Sponge
USE MOD_Sponge_Vars         ,ONLY: doSponge
USE MOD_Filter              ,ONLY: Filter_Pointer
USE MOD_Filter_Vars         ,ONLY: FilterType,FilterMat
USE MOD_Mesh_Vars           ,ONLY: nElems,nSides,lastInnerSide,firstMPISide_MINE,lastMPISide_MINE,firstMPISide_YOUR,lastMPISide_YOUR
USE NVTX
#if PARABOLIC
USE MOD_Lifting             ,ONLY: Lifting
#endif
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI                 ,ONLY: StartReceiveMPIData_GPU,StartSendMPIData_GPU,FinishExchangeMPIData
#endif /*USE_MPI*/
!@cuf USE MOD_Mesh_Vars     ,ONLY: d_sJ
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< Current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iSide,iElem
!==================================================================================================================================

! -----------------------------------------------------------------------------
! MAIN STEPS        []=FV only
! -----------------------------------------------------------------------------
! 1.   Filter solution vector
! 2.   Prolong to face (fill U_master/slave)
! 3.   Convert volume solution to primitive
! 4.   Finish communication of prolongated face data (2.)
! 5.   ConsToPrim of face data (U_master/slave)
! 6.   Lifting
! 7.   IF EDDYVISCOSITY: Prolong muSGS to face and send from slave to master
! 8.   Volume integral (DG only)
![9. ] FV volume integral
![10.] Volume integral (viscous contribution) if FV-blending
! 11.  Fill flux (Riemann solver) + surface integral
! 12.  Ut = -Ut
! 13.  Sponge and source terms
! 14.  Perform overintegration and apply Jacobian
! -----------------------------------------------------------------------------

! 2. Prolong the solution to the face integration points for flux computation (and do overlapping communication)
#if USE_MPI
CALL StartReceiveMPIData_GPU(d_U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE / U_slave: slave -> master
#endif
CALL ProlongToFaceCons(PP_N,d_U,d_U_master,d_U_slave,d_L_Minus,d_L_Plus)
#if USE_MPI
CALL StartSendMPIData_GPU(   d_U_slave,DataSizeSide,1,nSides,MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR / U_slave: slave -> master
#endif

! 3. Convert Volume solution to primitive for latency hiding
CALL ConsToPrim(PP_N,d_UPrim,d_U)

#if USE_MPI
! 4. Finish communication of face solution
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! U_slave: slave -> master
#endif

! 5. Convert face data from conservative to primitive variables
!    Attention: For FV with 2nd order reconstruction U_master/slave and therewith UPrim_master/slave are still only 1st order
! TODO: Linadv?
!CALL GetPrimitiveStateSurface(U_master,U_slave,UPrim_master,UPrim_slave)
CALL ConsToPrim(PP_N,nSides,d_UPrim_master,d_U_master) ! TODO: Skip non-filled MPI sides
CALL ConsToPrim(PP_N,nSides,d_UPrim_slave ,d_U_slave )

#if PARABOLIC
! 6. Lifting
! Compute the gradients using Lifting (BR1 scheme,BR2 scheme ...)
! The communication of the gradients is initialized within the lifting routines
CALL Lifting(d_UPrim,d_UPrim_master,d_UPrim_slave,t)
#endif /* PARABOLIC */

! 8. Compute volume integral contribution and add to Ut
CALL VolInt(d_Ut)

#if PARABOLIC && USE_MPI
! Complete send / receive for gradUx, gradUy, gradUz, started in the lifting routines
CALL FinishExchangeMPIData(6*nNbProcs,MPIRequest_gradU) ! gradUx,y,z: slave -> master
#endif /*PARABOLIC && USE_MPI*/

! 11. Fill flux and Surface integral
#if USE_MPI
CALL StartReceiveMPIData_GPU(d_Flux_slave, DataSizeSide, 1,nSides,MPIRequest_Flux( :,SEND),SendID=1)
CALL FillFlux(t,d_Flux_master,d_Flux_slave,d_U_master,d_U_slave,d_UPrim_master,d_UPrim_slave,doMPISides=.TRUE.)
CALL StartSendMPIData_GPU(   d_Flux_slave, DataSizeSide, 1,nSides,MPIRequest_Flux( :,RECV),SendID=1)
#endif
CALL FillFlux(t,d_Flux_master,d_Flux_slave,d_U_master,d_U_slave,d_UPrim_master,d_UPrim_slave,doMPISides=.FALSE.)
#if USE_MPI
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux )                       ! Flux_slave: master -> slave
#endif

! 11.5)
CALL SurfIntCons(PP_N,d_Flux_master,d_Flux_slave,d_Ut,d_L_HatMinus,d_L_hatPlus)

! 12. Swap to right sign :)
CALL VAX_GPU(nTotalU,d_Ut,-1.) ! Multiply array by -1

! 13. Compute source terms and sponge (in physical space, conversion to reference space inside routines)
!IF(doCalcSource) CALL CalcSource(Ut,t)
!IF(doSponge)     CALL Sponge(Ut)
!IF(doTCSource)   CALL TestcaseSource(Ut)

! 14. apply Jacobian
CALL ApplyJacobianCons(d_Ut,toPhysical=.TRUE.)

END SUBROUTINE DGTimeDerivative_weakForm



!==================================================================================================================================
!> Fills the solution array U with a initial solution provided by the ExactFunc subroutine though interpolation
!==================================================================================================================================
SUBROUTINE FillIni(Nloc,xGP,U)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Exactfunc     ,ONLY: ExactFunc
USE MOD_Mesh_Vars     ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: Nloc                                    !< Polynomial degree of solution
REAL,INTENT(IN)                 :: xGP(3,    0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)  !< Coordinates of Gauss-points
REAL,INTENT(OUT)                :: U(PP_nVar,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)  !< Solution array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
!==================================================================================================================================

! Evaluate the initial solution at the nodes and fill the solutin vector U.
DO iElem=1,nElems
  DO k=0,ZDIM(Nloc); DO j=0,Nloc; DO i=0,Nloc
    CALL ExactFunc(IniExactFunc,0.,xGP(1:3,i,j,k,iElem),U(:,i,j,k,iElem))
  END DO; END DO; END DO
END DO
END SUBROUTINE FillIni



!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeDG()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_DG_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(D)
SDEALLOCATE(D_T)
SDEALLOCATE(D_Hat)
SDEALLOCATE(D_Hat_T)
#if SPLIT_DG
SDEALLOCATE(DVolSurf)
#endif
SDEALLOCATE(L_HatMinus)
SDEALLOCATE(L_HatPlus)
SDEALLOCATE(U)
SDEALLOCATE(Ut)
SDEALLOCATE(U_master)
SDEALLOCATE(U_slave)
SDEALLOCATE(Flux_master)
SDEALLOCATE(Flux_slave)
SDEALLOCATE(UPrim)
SDEALLOCATE(UPrim_master)
SDEALLOCATE(UPrim_slave)
SDEALLOCATE(UPrim_boundary)
SDEALLOCATE(d_f)
SDEALLOCATE(d_g)
SDEALLOCATE(d_h)

DGInitIsDone = .FALSE.
END SUBROUTINE FinalizeDG


END MODULE MOD_DG
