!=================================================================================================================================
! Copyright (c) 2010-2022  Prof. Claus-Dieter Munz
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
!> \brief This module contains routines for computing the geometries volume and surface metric terms.
!>
!> Compute the volume and surface metric terms:
!>     Metrics_fTilde(n=1:3,i,j,k,iElem)=Ja_n^1
!>     Metrics_gTilde(n=1:3,i,j,k,iElem)=Ja_n^2
!>     Metrics_hTilde(n=1:3,i,j,k,iElem)=Ja_n^3
!>
!>   Per Element we do:
!>   1.) a.) Preparation: the geometry (equidistant nodal basis, NGeo+1 points/dir) is interpolated to a high precision
!>           mapping X_n(xi_i) using a Chebyshev-Lobatto basis and stored in XCL_NGeo(1:3,i,j,k,iElem) i,j,k=[0:NGeo]
!>       b.) Computing the gradients: compute the derivative of the mapping XCL_NGeo in \f$ (xi_1,xi_2,xi_3) \f$ direction,
!>           using a polynomial derivative Matrix at degree NGeo.
!>       c.) Computing the Jacobian: compute Jacobian JRef at a degree of NGeoRef=3*NGeo (exact).
!>                                   For this gradients have to be interpolated to NGeoRef first.
!>                                   Then project JRef down to degree N. Finally check for negative Jacobians.
!>       d.) For computing Ja the gradients at degree N are required: if N>=NGeo directly interpolate dXCL_NGeo to dXCL_N,
!>                                                                    else compute dXCL_N from XCL_N directly.
!>
!>   2.) for each direction n
!>       a.) compute the nth vector and for each Chebyshev point (:,i,j,k)
!>          \f$(dXCL_n^1,dXCL_n^2,dXCL_n^3)^T=(X_l grad_xi (X_m) )\f$ for n=1,2,3 and (n,m,l) cyclic
!>       b.) interpolate the dXCL_n vector defined primarily on (NGeo+1)x(NGeo+1)x(NGeo+1) Chebyshev-Lobatto points to
!>             (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points and write to Ja_n(1:3,i,j,k) i,j,k=[0:N]
!>       c.) compute the curl of vector Ja_n(1:3,i,j,k) using the derivative Matrix DCL_N [NxN]
!>       d.) interpolate from (N+1)x(N+1)x(N+1) Chebyshev-Lobatto points to  Gauss-Points (N+1)x(N+1)x(N+1) (exact!)
!>       e.) store Ja_n in the Metrics arrays
!>
!>   3.) Compute the surface metrics (normal/tangential vectors, surface area) from volume metrics for each side.
!>
!>  Special case if non-conforming meshes with octree mappings are used. Then compute ALL volume quantities on tree (macro element)
!>  level and interpolate down to small actual elements. This will ensure watertight meshes and free-stream preservation.
!==================================================================================================================================
MODULE MOD_Metrics
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE BuildCoords
  MODULE PROCEDURE BuildCoords
END INTERFACE

INTERFACE CalcMetrics
  MODULE PROCEDURE CalcMetrics
END INTERFACE

INTERFACE CalcSurfMetrics
  MODULE PROCEDURE CalcSurfMetrics
END INTERFACE

INTERFACE SurfMetricsFromJa
  MODULE PROCEDURE SurfMetricsFromJa
END INTERFACE

PUBLIC::BuildCoords
PUBLIC::CalcMetrics
PUBLIC::CalcSurfMetrics
PUBLIC::SurfMetricsFromJa
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> This routine takes the equidistant node coordinats of the mesh (on NGeo+1 points) and uses them to build the coordinates
!> of solution/interpolation points of type NodeType on polynomial degree Nloc (Nloc+1 points per direction).
!> The coordinates (for a non-conforming mesh) can also be built from an octree if the mesh is based on a conforming baseline mesh.
!==================================================================================================================================
SUBROUTINE BuildCoords(NodeCoords,NodeType,Nloc,VolumeCoords,TreeCoords)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars          ,ONLY: NGeo,nElems
USE MOD_Mesh_Vars          ,ONLY: ElemToTree,xiMinMax,nTrees,NGeoTree
USE MOD_Interpolation_Vars ,ONLY: NodeTypeCL,NodeTypeVISU
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights
#if (PP_dim == 3)
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D_XYZ
#else
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D_XYZ
#endif
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)            :: NodeCoords(3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems)       !< Equidistant mesh coordinates
CHARACTER(LEN=255),INTENT(IN) :: NodeType                                              !< Type of node that should be converted to
INTEGER,INTENT(IN)            :: Nloc                                                  !< Convert to Nloc+1 points per direction
REAL,INTENT(OUT)              :: VolumeCoords(3,0:Nloc,0:Nloc,0:ZDIM(Nloc),nElems)     !< OUT: Coordinates of solution/interpolation
                                                                                       !< points
REAL,INTENT(INOUT),OPTIONAL   :: TreeCoords(3,0:NGeoTree,0:NGeoTree,0:NGeoTree,nTrees) !< coordinates of nodes of tree-elements
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                       :: i,iElem
REAL                          :: XCL_Nloc(3,0:Nloc,0:Nloc,0:Nloc)
REAL,DIMENSION(0:Nloc,0:Nloc) :: Vdm_xi_N,Vdm_eta_N
#if (PP_dim == 3)
REAL,DIMENSION(0:Nloc,0:Nloc) :: Vdm_zeta_N
#endif
REAL                          :: Vdm_EQNGeo_CLNloc(0:Nloc ,0:NGeo)
REAL                          :: Vdm_CLNloc_Nloc  (0:Nloc ,0:Nloc)
REAL                          :: xi0(3),dxi(3),length(3)
REAL                          :: xiCL_Nloc(0:Nloc),wBaryCL_Nloc(0:Nloc)
REAL                          :: xiNloc(0:Nloc)
!==================================================================================================================================

CALL GetVandermonde(    NGeo, NodeTypeVISU, NLoc, NodeTypeCL, Vdm_EQNGeo_CLNloc,  modal=.FALSE.)
CALL GetVandermonde(    Nloc, NodeTypeCL  , Nloc, NodeType  , Vdm_CLNloc_Nloc,     modal=.FALSE.)

! NOTE: Transform intermediately to CL points, to be consistent with metrics being built with CL
!       Important for curved meshes if NGeo<N, no effect for N>=NGeo

!1.a) Transform from EQUI_NGeo to solution points on Nloc
IF(PRESENT(TreeCoords))THEN
  CALL GetNodesAndWeights(Nloc, NodeTypeCL  , xiCL_Nloc  , wIPBary=wBaryCL_Nloc)
  CALL GetNodesAndWeights(Nloc, NodeType  ,   xiNloc)
  DO iElem=1,nElems
    xi0   =xiMinMax(:,1,iElem)
    length=xiMinMax(:,2,iElem)-xi0
    CALL ChangeBasisVolume(3,NGeo,Nloc,Vdm_EQNGeo_CLNloc,TreeCoords(:,:,:,:,ElemToTree(iElem)),XCL_Nloc)
    DO i=0,Nloc
      dxi=0.5*(xiNloc(i)+1.)*length
      CALL LagrangeInterpolationPolys(xi0(1) + dxi(1),Nloc,xiCL_Nloc,wBaryCL_Nloc,Vdm_xi_N(  i,:))
      CALL LagrangeInterpolationPolys(xi0(2) + dxi(2),Nloc,xiCL_Nloc,wBaryCL_Nloc,Vdm_eta_N( i,:))
#if (PP_dim == 3)
      CALL LagrangeInterpolationPolys(xi0(3) + dxi(3),Nloc,xiCL_Nloc,wBaryCL_Nloc,Vdm_zeta_N(i,:))
#endif
    END DO
#if (PP_dim == 3)
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,XCL_Nloc,VolumeCoords(:,:,:,:,iElem))
#else
    CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,XCL_Nloc(:,:,:,0),VolumeCoords(:,:,:,0,iElem))
#endif
  END DO
ELSE
  Vdm_EQNGeo_CLNloc=MATMUL(Vdm_CLNloc_Nloc,Vdm_EQNGeo_CLNloc)
  DO iElem=1,nElems
    CALL ChangeBasisVolume(3,NGeo,Nloc,Vdm_EQNGeo_CLNloc,NodeCoords(:,:,:,:,iElem),VolumeCoords(:,:,:,:,iElem))
  END DO
END IF

#if (PP_dim == 2)
! Nullify coordinates in the third dimension
VolumeCoords(3,:,:,:,:) = 0.
#endif


END SUBROUTINE BuildCoords

!==================================================================================================================================
!> This routine computes the geometries volume metric terms.
!==================================================================================================================================
SUBROUTINE CalcMetrics()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Basis              ,ONLY: LagrangeInterpolationPolys
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_Interpolation_Vars
USE MOD_Interpolation      ,ONLY: GetVandermonde,GetNodesAndWeights,GetDerivativeMatrix
USE MOD_Mesh_Vars          ,ONLY: NGeo,NGeoRef,nElems,offsetElem
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,dXCL_N
USE MOD_Mesh_Vars          ,ONLY: sJ,detJac_Ref,Ja_Face
USE MOD_Mesh_Vars          ,ONLY: NodeCoords,TreeCoords,Elem_xGP
USE MOD_Mesh_Vars          ,ONLY: ElemToTree,xiMinMax,interpolateFromTree
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2,SurfElem,Face_xGP
USE MOD_Mesh_Vars          ,ONLY: scaledJac
USE MOD_Output_Vars        ,ONLY: doPrintStatusLine
USE MOD_ReadInTools        ,ONLY: prms
#if PP_dim == 3
USE MOD_Mesh_Vars          ,ONLY: crossProductMetrics
#endif
#if FV_ENABLED
USE MOD_Mesh_Vars          ,ONLY: sJ_master,sJ_slave
USE MOD_ProlongToFace1     ,ONLY: ProlongToFace1
USE MOD_FillMortar1        ,ONLY: U_Mortar1
#endif /*FV_ENABLED*/
#if (PP_dim == 3)
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D_XYZ
#else
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D_XYZ
#endif
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY: firstMPISide_MINE,firstMPISide_YOUR,lastMPISide_YOUR,nSides
USE MOD_MPI_Vars           ,ONLY: nNbProcs
USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
#if PP_dim == 3
INTEGER :: q
#endif
INTEGER :: ll
! Jacobian on CL N and NGeoRef
REAL    :: DetJac_N( 1,0:PP_N,   0:PP_N,   0:PP_NZ,nElems)
REAL    :: tmp(      1,0:NGeoRef,0:NGeoRef,0:ZDIM(NGeoRef))
!REAL    :: tmp2(     1,0:NGeo,0:NGeo,0:NGeo)
! interpolation points and derivatives on CL N
REAL    :: XCL_N(      3,  0:PP_N,0:PP_N,0:PP_NZ)               ! mapping X(xi) P\in N
REAL    :: XCL_Ngeo(   3,  0:Ngeo,0:Ngeo,0:ZDIM(NGeo))          ! mapping X(xi) P\in NGeo
REAL    :: XCL_N_quad( 3,  0:PP_N,0:PP_N,0:PP_NZ)               ! mapping X(xi) P\in N
REAL    :: dXCL_NGeo(  3,3,0:NGeo,0:NGeo,0:ZDIM(NGeo))          ! jacobi matrix on CL NGeo
REAL    :: dX_NGeoRef( 3,3,0:NGeoRef,0:NGeoRef,0:ZDIM(NGeoRef)) ! jacobi matrix on SOL NGeoRef

#if PP_dim == 3
REAL    :: R_CL_N(     3,3,0:PP_N,0:PP_N,0:PP_NZ)    ! buffer for metric terms, uses XCL_N,dXCL_N
#endif
REAL    :: JaCL_N(     3,3,0:PP_N,0:PP_N,0:PP_NZ)    ! metric terms P\in N
REAL    :: JaCL_N_quad(3,3,0:PP_N,0:PP_N,0:PP_NZ)    ! metric terms P\in N

! Polynomial derivativion matrices
REAL    :: DCL_NGeo(0:NGeo,0:NGeo)
REAL    :: DCL_N(   0:PP_N,0:PP_N)

! Vandermonde matrices (N_OUT,N_IN)
REAL    :: Vdm_EQNGeo_CLNGeo( 0:NGeo   ,0:NGeo)
REAL    :: Vdm_CLNGeo_NGeoRef(0:NGeoRef,0:NGeo)
REAL    :: Vdm_NGeoRef_N(     0:PP_N   ,0:NGeoRef)
REAL    :: Vdm_CLNGeo_CLN(    0:PP_N   ,0:NGeo)
REAL    :: Vdm_CLN_N(         0:PP_N   ,0:PP_N)

! 3D Vandermonde matrices and lengths,nodes,weights
REAL,DIMENSION(0:NGeoRef,0:NGeoRef) :: Vdm_xi_Ref,Vdm_eta_Ref
REAL,DIMENSION(0:PP_N   ,0:PP_N)    :: Vdm_xi_N  ,Vdm_eta_N
#if PP_dim == 3
REAL,DIMENSION(0:NGeoRef,0:NGeoRef) :: Vdm_zeta_Ref
REAL,DIMENSION(0:PP_N   ,0:PP_N)    :: Vdm_zeta_N
#endif
REAL    :: xiRef( 0:NGeoRef),wBaryRef( 0:NGeoRef)
REAL    :: xiCL_N(0:PP_N)   ,wBaryCL_N(0:PP_N)
REAL    :: xi0(3),dxi(3),length(3)

#if USE_MPI
INTEGER           :: MPIRequest_Geo(nNbProcs,2)
REAL,ALLOCATABLE  :: Geo(:,:,:,:,:)
#endif

! Output
REAL              :: percent
CHARACTER(LEN=20) :: fmtName
!==================================================================================================================================
! Prerequisites
Metrics_fTilde=0.
Metrics_gTilde=0.
Metrics_hTilde=0.

! 1.) Compute the Jacobian, at this point first the Vandermonde and D matrices are prepared

! 1.a) Mapping NodeCoords defined on equidistant points (NodeTypeVISU) with NGeo (EQNGeo) to CL points with NGeo (CLNGeo) -> l_k^{EQUI} (s^{CLN})
CALL GetVandermonde(    NGeo   , NodeTypeVISU, NGeo    , NodeTypeCL, Vdm_EQNGeo_CLNGeo , modal=.FALSE.)

! 1.b) Switching from NGeo to PP_N still on CL points
CALL GetVandermonde(    NGeo   , NodeTypeCL  , PP_N    , NodeTypeCL, Vdm_CLNGeo_CLN    , modal=.FALSE.)

! 1.c) dXCL_NGeo: l_k'^{EQUI} (s^{CLN})
! Get the derivative matrix at the NGeo points
CALL GetDerivativeMatrix(NGeo  , NodeTypeCL  , DCL_NGeo)

! 1.d) Interpolation from CL points to LG/LGL points (dX/dXi): Jacobian of the mapping
CALL GetVandermonde(    NGeo   , NodeTypeCL  , NGeoRef , NodeType  , Vdm_CLNGeo_NGeoRef, modal=.FALSE.)

! 1.f) Projection of the Jacobian from NGeo to PP_N
! Switch from NGeoRef to PP_N still on LG/LGL points but we use a modal Vandermonde
! This has to be done to conserve the Jacobian if N_out>N_in as PP_N > NGeo
CALL GetVandermonde(    NGeoRef, NodeType    , PP_N    , NodeType  , Vdm_NGeoRef_N     , modal=.TRUE.)
! Still 1.f) but only required for interpolate from tree
! Compute the integration nodes and weights for NGeoRef on LG/LGL points
CALL GetNodesAndWeights(NGeoRef, NodeType    , xiRef   , wIPBary=wBaryRef)


! 2.) Compute metrics, at this point first the Vandermonde and D matrices are prepared
! Get the derivative matrix at the CL points with PP_N points
CALL GetDerivativeMatrix(PP_N  , NodeTypeCL  , DCL_N)

! 2.d) derivatives (dXCL) by projection or by direct derivation (D_CL):
! Interpolation from CL to LG/LGL points (PP_N)
CALL GetVandermonde(    PP_N   , NodeTypeCL  , PP_N    , NodeType,   Vdm_CLN_N         , modal=.FALSE.)
! Compute the integration nodes and weights for PP_N on CL points
CALL GetNodesAndWeights(PP_N   , NodeTypeCL  , xiCL_N  , wIPBary=wBaryCL_N)

! Outer loop over all elements
detJac_Ref=0.
dXCL_N=0.
DO iElem=1,nElems
  !1.a) Transform NodeCoords [\Gamma (s_k^{EQNGeo})] from EQNGeo to CLNGeo using Vdm_EQNGeo_CLNGeo to obtain XCL_NGeo
  ! Equation: I_N \Gamma (s^{CLNGeo}) = \sum_{k=0}^{NGeo} \Gamma (s_k^{EQNGeo}) l_k^{EQNGeo} (s^{CLNGeo})
  IF(interpolateFromTree)THEN
    xi0   =xiMinMax(:,1,iElem)
    length=xiMinMax(:,2,iElem)-xi0
#if (PP_dim == 2)
    length(3) = 1.
#endif
    CALL ChangeBasisVolume(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,TreeCoords(:,:,:,:,ElemToTree(iElem)),XCL_NGeo)
  ELSE
    CALL ChangeBasisVolume(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem)            ,XCL_NGeo)
  END IF
  !1.b) Switching XCL_NGeo from NGeo to PP_N to obtain XCL_N using Vdm_CLNGeo_CLN; still on CL points
  CALL   ChangeBasisVolume(3,NGeo,PP_N,Vdm_CLNGeo_CLN,   XCL_NGeo                             ,XCL_N)

  !1.c) Jacobi Matrix of dX/dXi_dd(XCL_NGeo_nn): dXCL_NGeo(dd,nn,i,j,k))
  ! d/dXi_dd(XCL_NGeo_nn) = DCL_NGeo
  dXCL_NGeo=0.
  DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
    ! Matrix-vector multiplication
    DO ll=0,NGeo
      dXCL_NGeo(1,1:PP_dim,i,j,k)=dXCL_NGeo(1,1:PP_dim,i,j,k) + DCL_NGeo(i,ll)*XCL_NGeo(1:PP_dim,ll,j,k)
      dXCL_NGeo(2,1:PP_dim,i,j,k)=dXCL_NGeo(2,1:PP_dim,i,j,k) + DCL_NGeo(j,ll)*XCL_NGeo(1:PP_dim,i,ll,k)
#if (PP_dim == 3)
      dXCL_NGeo(3,:,i,j,k)=dXCL_NGeo(3,:,i,j,k) + DCL_NGeo(k,ll)*XCL_NGeo(:,i,j,ll)
#endif
    END DO !l=0,N
  END DO; END DO; END DO !i,j,k=0,NGeo

  !1.d) Compute the Jacobian of the mapping (dX/dXi) from CL points to LG/LGL points: dX_NGeoRef
  CALL ChangeBasisVolume(PP_dim,NGeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo(1:PP_dim,1,:,:,:),dX_NGeoRef(1:PP_dim,1,:,:,:))
  CALL ChangeBasisVolume(PP_dim,NGeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo(1:PP_dim,2,:,:,:),dX_NGeoRef(1:PP_dim,2,:,:,:))
#if (PP_dim == 3)
  CALL ChangeBasisVolume(3,NGeo,NGeoRef,Vdm_CLNGeo_NGeoRef,dXCL_NGeo(:,3,:,:,:),dX_NGeoRef(:,3,:,:,:))
#endif
  !1.e) Jacobian on NGeo = a_i (a_j x a_k) = d(X_1)/dXi_i (d(X_2)/dXi_j x d(X_3)/dXi_k): detJac_Ref
  DO k=0,ZDIM(NGeoRef)  ; DO j=0,NGeoRef; DO i=0,NGeoRef
#if (PP_dim == 3)
    detJac_Ref(1,i,j,k,iElem)=detJac_Ref(1,i,j,k,iElem) &
      + dX_NGeoRef(1,1,i,j,k)*(dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(3,3,i,j,k) - dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(2,3,i,j,k))  &
      + dX_NGeoRef(2,1,i,j,k)*(dX_NGeoRef(3,2,i,j,k)*dX_NGeoRef(1,3,i,j,k) - dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(3,3,i,j,k))  &
      + dX_NGeoRef(3,1,i,j,k)*(dX_NGeoRef(1,2,i,j,k)*dX_NGeoRef(2,3,i,j,k) - dX_NGeoRef(2,2,i,j,k)*dX_NGeoRef(1,3,i,j,k))
#else
      detJac_Ref(1,i,j,k,iElem)=detJac_Ref(1,i,j,k,iElem) &
        + dX_NGeoRef(1,1,i,j,k)*dX_NGeoRef(2,2,i,j,k) - dX_NGeoRef(2,1,i,j,k)*dX_NGeoRef(1,2,i,j,k)
#endif
  END DO; END DO; END DO !i,j,k=0,NGeoRef

  !1.f) Project the Jacobian from NGeo to PP_N using the modal Vandermonde to guarantee conservativity if N > NGeo
  IF(interpolateFromTree)THEN
    !interpolate detJac to the GaussPoints
    DO i=0,NGeoRef
      dxi=0.5*(xiRef(i)+1.)*Length
      CALL LagrangeInterpolationPolys(xi0(1) + dxi(1),NGeoRef,xiRef,wBaryRef,Vdm_xi_Ref(  i,:))
      CALL LagrangeInterpolationPolys(xi0(2) + dxi(2),NGeoRef,xiRef,wBaryRef,Vdm_eta_Ref( i,:))
#if (PP_dim == 3)
      CALL LagrangeInterpolationPolys(xi0(3) + dxi(3),NGeoRef,xiRef,wBaryRef,Vdm_zeta_Ref(i,:))
#endif
    END DO
    tmp=DetJac_Ref(:,:,:,:,iElem)
#if (PP_dim == 3)
    CALL ChangeBasis3D_XYZ(1,NGeoRef,NGeoRef,Vdm_xi_Ref,Vdm_eta_Ref,Vdm_zeta_Ref,&
                           tmp,DetJac_Ref(:,:,:,:,iElem))
#else
    CALL ChangeBasis2D_XYZ(1,NGeoRef,NGeoRef,Vdm_xi_Ref,Vdm_eta_Ref,&
                           tmp(:,:,:,0),DetJac_Ref(:,:,:,0,iElem))
#endif
  END IF
  CALL ChangeBasisVolume(1,NGeoRef,PP_N,Vdm_NGeoRef_N,DetJac_Ref(:,:,:,:,iElem),DetJac_N(:,:,:,:,iElem))

  ! assign to global Variable sJ
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    sJ(i,j,k,iElem,0)=1./DetJac_N(1,i,j,k,iElem)
  END DO; END DO; END DO !i,j,k=0,PP_N

  ! check for negative Jacobians
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    IF(detJac_N(1,i,j,k,iElem).LE.0.)&
      WRITE(UNIT_stdOut,*) 'Negative Jacobian found on Gauss point. Coords:', Elem_xGP(:,i,j,k,iElem)
    ! check scaled Jacobians
    scaledJac(i,j,k,iElem)=detJac_N(1,i,j,k,iElem)/MAXVAL(detJac_N(1,:,:,:,iElem))
    IF(scaledJac(i,j,k,iElem).LT.0.01) THEN
      WRITE(UNIT_stdOut,*) 'Too small scaled Jacobians found (CL/Gauss):', scaledJac(i,j,k,iElem)
      CALL Abort(__STAMP__,&
        'Scaled Jacobian lower then tolerance in global element:',iElem+offsetElem)
    END IF
  END DO; END DO; END DO !i,j,k=0,N

  !2.a) Jacobi Matrix of dX/dxi_dd(X_nn): covariant basis vector a_dd = dXCL_N(dd,nn,i,j,k))
  ! N>=NGeo: interpolate from dXCL_NGeo (default)
  ! N< NGeo: directly derive XCL_N
  IF(PP_N.GE.NGeo)THEN !compute first derivative on NGeo and then interpolate
    CALL ChangeBasisVolume(PP_dim,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(1:PP_dim,1,:,:,:),dXCL_N(1:PP_dim,1,:,:,:,iElem))
    CALL ChangeBasisVolume(PP_dim,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(1:PP_dim,2,:,:,:),dXCL_N(1:PP_dim,2,:,:,:,iElem))
#if (PP_dim == 3)
    CALL ChangeBasisVolume(3,NGeo,PP_N,Vdm_CLNGeo_CLN,dXCL_NGeo(:,3,:,:,:),dXCL_N(:,3,:,:,:,iElem))
#endif
  ELSE  !N<NGeo: first interpolate and then compute derivative (important if curved&periodic)
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      ! Matrix-vector multiplication
      ASSOCIATE(dXCL => dXCL_N(:,:,i,j,k,iElem))
      DO ll=0,PP_N
        dXCL(1,1:PP_dim)=dXCL(1,1:PP_dim) + DCL_N(i,ll)*XCL_N(1:PP_dim,ll,j,k)
        dXCL(2,1:PP_dim)=dXCL(2,1:PP_dim) + DCL_N(j,ll)*XCL_N(1:PP_dim,i,ll,k)
#if (PP_dim == 3)
        dXCL(3,:)=dXCL(3,:) + DCL_N(k,ll)*XCL_N(:,i,j,ll)
#endif
      END DO !l=0,N
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  END IF !N>=NGeo

  !2.b) Compute metrics: Jacobian * contravariant basis vector = J*a^i
  JaCL_N=0.
#if (PP_dim == 2)
  ! No need to differentiate between curl and cross product metrics in 2D, we can directly use the calculated derivatives
  ! a^1 = JaCL_N(1,1:PP_dim,:,:,:); a^2 = JaCL_N(2,1:PP_dim,:,:,:)
  DO k=0,0; DO j=0,PP_N; DO i=0,PP_N
    JaCL_N(1,1,i,j,k)=  dXCL_N(2,2,i,j,k,iElem)
    JaCL_N(1,2,i,j,k)=- dXCL_N(2,1,i,j,k,iElem)
    JaCL_N(2,1,i,j,k)=- dXCL_N(1,2,i,j,k,iElem)
    JaCL_N(2,2,i,j,k)=  dXCL_N(1,1,i,j,k,iElem)
  END DO; END DO; END DO !i,j,k=0,N
#else
  IF(crossProductMetrics)THEN
    ! exact (cross-product) form
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(dXCL => dXCL_N(:,:,i,j,k,iElem))
      ! exact (cross-product) form
      ! Ja(:)^nn = ( d/dxi_(nn+1) XCL_N(:) ) x (d/xi_(nn+2) XCL_N(:))
      !
      ! JaCL_N(dd,nn) = dXCL_N(dd+1,nn+1)*dXCL_N(dd+2,nn+2) -dXCL_N(dd+1,nn+2)*dXCL_N(dd+2,nn+1)
      JaCL_N(1,1,i,j,k)=dXCL(2,2)*dXCL(3,3) - dXCL(2,3)*dXCL(3,2)
      JaCL_N(2,1,i,j,k)=dXCL(3,2)*dXCL(1,3) - dXCL(3,3)*dXCL(1,2)
      JaCL_N(3,1,i,j,k)=dXCL(1,2)*dXCL(2,3) - dXCL(1,3)*dXCL(2,2)
      JaCL_N(1,2,i,j,k)=dXCL(2,3)*dXCL(3,1) - dXCL(2,1)*dXCL(3,3)
      JaCL_N(2,2,i,j,k)=dXCL(3,3)*dXCL(1,1) - dXCL(3,1)*dXCL(1,3)
      JaCL_N(3,2,i,j,k)=dXCL(1,3)*dXCL(2,1) - dXCL(1,1)*dXCL(2,3)
      JaCL_N(1,3,i,j,k)=dXCL(2,1)*dXCL(3,2) - dXCL(2,2)*dXCL(3,1)
      JaCL_N(2,3,i,j,k)=dXCL(3,1)*dXCL(1,2) - dXCL(3,2)*dXCL(1,1)
      JaCL_N(3,3,i,j,k)=dXCL(1,1)*dXCL(2,2) - dXCL(1,2)*dXCL(2,1)
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
  ELSE ! curl metrics
    ! invariant curl form, as cross product: R^dd = 1/2( XCL_N(:) x (d/dxi_dd XCL_N(:)))
    !
    !R_CL_N(dd,nn)=1/2*( XCL_N(nn+2)* d/dxi_dd XCL_N(nn+1) - XCL_N(nn+1)* d/dxi_dd XCL_N(nn+2))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(dXCL => dXCL_N(:,:,i,j,k,iElem))
      R_CL_N(:,1,i,j,k)=0.5*(XCL_N(3,i,j,k)*dXCL(:,2) - XCL_N(2,i,j,k)*dXCL(:,3) )
      R_CL_N(:,2,i,j,k)=0.5*(XCL_N(1,i,j,k)*dXCL(:,3) - XCL_N(3,i,j,k)*dXCL(:,1) )
      R_CL_N(:,3,i,j,k)=0.5*(XCL_N(2,i,j,k)*dXCL(:,1) - XCL_N(1,i,j,k)*dXCL(:,2) )
      END ASSOCIATE
    END DO; END DO; END DO !i,j,k=0,N
    ! Metrics are the curl of R:  Ja(:)^nn = -(curl R_CL(:,nn))
    ! JaCL_N(dd,nn)= -[d/dxi_(dd+1) RCL(dd+2,nn) - d/dxi_(dd+2) RCL(dd+1,nn) ]
    !              =   d/dxi_(dd+2) RCL(dd+1,nn) - d/dxi_(dd+1) RCL(dd+2,nn)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      ASSOCIATE(JaCL => JaCL_N(:,:,i,j,k))
      DO q=0,PP_N
        JaCL(1,:)=JaCL(1,:) - DCL_N(j,q)*R_CL_N(3,:,i,q,k)
        JaCL(2,:)=JaCL(2,:) - DCL_N(k,q)*R_CL_N(1,:,i,j,q)
        JaCL(3,:)=JaCL(3,:) - DCL_N(i,q)*R_CL_N(2,:,q,j,k)
      END DO!q=0,PP_N
      DO q=0,PP_N
        JaCL(1,:)=JaCL(1,:) + DCL_N(k,q)*R_CL_N(2,:,i,j,q)
        JaCL(2,:)=JaCL(2,:) + DCL_N(i,q)*R_CL_N(3,:,q,j,k)
        JaCL(3,:)=JaCL(3,:) + DCL_N(j,q)*R_CL_N(1,:,i,q,k)
      END DO!q=0,PP_N
      END ASSOCIATE
! same with only one loop, gives different roundoff ...
!      DO q=0,PP_N
!        JaCL_N(1,:,i,j,k)=JaCL_N(1,:,i,j,k) - DCL_N(j,q)*R_CL_N(3,:,i,q,k) + DCL_N(k,q)*R_CL_N(2,:,i,j,q)
!        JaCL_N(2,:,i,j,k)=JaCL_N(2,:,i,j,k) - DCL_N(k,q)*R_CL_N(1,:,i,j,q) + DCL_N(i,q)*R_CL_N(3,:,q,j,k)
!        JaCL_N(3,:,i,j,k)=JaCL_N(3,:,i,j,k) - DCL_N(i,q)*R_CL_N(2,:,q,j,k) + DCL_N(j,q)*R_CL_N(1,:,i,q,k)
!      END DO!q=0,PP_N
    END DO; END DO; END DO !i,j,k=0,N
  END IF !crossProductMetrics
#endif

  !2.c) Interpolate J*a^i (JaCL) from CL to LG/LGL points: Metrics_fTilde for i, ...
  ! Calc metrics at the surface: NormVec, TangVec, SurfElem, ...
  IF(interpolateFromTree)THEN
    ! interpolate Metrics from Cheb-Lobatto N on tree level onto GaussPoints N on quad level
    DO i=0,PP_N
      dxi=0.5*(xGP(i)+1.)*length
      CALL LagrangeInterpolationPolys(xi0(1) + dxi(1),PP_N,xiCL_N,wBaryCL_N,Vdm_xi_N(  i,:))
      CALL LagrangeInterpolationPolys(xi0(2) + dxi(2),PP_N,xiCL_N,wBaryCL_N,Vdm_eta_N( i,:))
#if (PP_dim == 3)
      CALL LagrangeInterpolationPolys(xi0(3) + dxi(3),PP_N,xiCL_N,wBaryCL_N,Vdm_zeta_N(i,:))
#endif
    END DO
#if (PP_dim == 3)
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(1,:,:,:,:),Metrics_fTilde(:,:,:,:,iElem,0))
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(2,:,:,:,:),Metrics_gTilde(:,:,:,:,iElem,0))
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(3,:,:,:,:),Metrics_hTilde(:,:,:,:,iElem,0))
#else
    CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,JaCL_N(1,:,:,:,0),Metrics_fTilde(:,:,:,0,iElem,0))
    CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,JaCL_N(2,:,:,:,0),Metrics_gTilde(:,:,:,0,iElem,0))
#endif
    ! for the metrics and the jacobian, we have to take into account the level !!!!!
    Metrics_fTilde(:,:,:,:,iElem,0)=(length(1)/2.)**2*Metrics_fTilde(:,:,:,:,iElem,0)
    Metrics_gTilde(:,:,:,:,iElem,0)=(length(2)/2.)**2*Metrics_gTilde(:,:,:,:,iElem,0)
#if (PP_dim == 3)
    Metrics_hTilde(:,:,:,:,iElem,0)=(length(3)/2.)**2*Metrics_hTilde(:,:,:,:,iElem,0)
#endif
    sJ(:,:,:,iElem,0)=(2.**PP_dim/PRODUCT(length))*sJ(:,:,:,iElem,0) ! scale down sJ

    ! interpolate Metrics and grid to Cheb-Lobatto on quadrant level for Surface metrics
    DO i=0,PP_N
      dxi=0.5*(xiCL_N(i)+1.)*length
      CALL LagrangeInterpolationPolys(xi0(1) + dxi(1),PP_N,xiCL_N,wBaryCL_N,Vdm_xi_N(  i,:))
      CALL LagrangeInterpolationPolys(xi0(2) + dxi(2),PP_N,xiCL_N,wBaryCL_N,Vdm_eta_N( i,:))
#if (PP_dim == 3)
      CALL LagrangeInterpolationPolys(xi0(3) + dxi(3),PP_N,xiCL_N,wBaryCL_N,Vdm_zeta_N(i,:))
#endif
    END DO
#if (PP_dim == 3)
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,XCL_N            ,XCL_N_quad            )
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(1,:,:,:,:),JaCL_N_quad(1,:,:,:,:))
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(2,:,:,:,:),JaCL_N_quad(2,:,:,:,:))
    CALL ChangeBasis3D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,Vdm_zeta_N,JaCL_N(3,:,:,:,:),JaCL_N_quad(3,:,:,:,:))
#else
    CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,XCL_N(   :,:,:,0),XCL_N_quad(   :,:,:,0))
    CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,JaCL_N(1,:,:,:,0),JaCL_N_quad(1,:,:,:,0))
    CALL ChangeBasis2D_XYZ(3,PP_N,PP_N,Vdm_xi_N,Vdm_eta_N,JaCL_N(2,:,:,:,0),JaCL_N_quad(2,:,:,:,0))
#endif
    !TODO: scale Ja for anisotropic
    JaCL_N_quad(:,1,:,:,:)=(length(2)*length(3)/(2.**(PP_dim-1)))*JaCL_N_quad(:,1,:,:,:)
    JaCL_N_quad(:,2,:,:,:)=(length(1)*length(3)/(2.**(PP_dim-1)))*JaCL_N_quad(:,2,:,:,:)
#if (PP_dim == 3)
    JaCL_N_quad(:,3,:,:,:)=(length(1)*length(2)/4.)*JaCL_N_quad(:,3,:,:,:)
#endif
    CALL CalcSurfMetrics(PP_N,FV_SIZE,JaCL_N_quad,XCL_N_quad,Vdm_CLN_N,iElem,&
                         NormVec,TangVec1,TangVec2,SurfElem,Face_xGP,Ja_Face)
  ELSE
    ! interpolate Metrics from Cheb-Lobatto N onto GaussPoints N
    CALL ChangeBasisVolume(PP_dim,PP_N,PP_N,Vdm_CLN_N,JaCL_N(1,1:PP_dim,:,:,:),Metrics_fTilde(1:PP_dim,:,:,:,iElem,0))
    CALL ChangeBasisVolume(PP_dim,PP_N,PP_N,Vdm_CLN_N,JaCL_N(2,1:PP_dim,:,:,:),Metrics_gTilde(1:PP_dim,:,:,:,iElem,0))
#if (PP_dim == 3)
    CALL ChangeBasisVolume(3,PP_N,PP_N,Vdm_CLN_N,JaCL_N(3,:,:,:,:),Metrics_hTilde(:,:,:,:,iElem,0))
#endif
    CALL CalcSurfMetrics(PP_N,FV_SIZE,JaCL_N,XCL_N,Vdm_CLN_N,iElem,&
                         NormVec,TangVec1,TangVec2,SurfElem,Face_xGP,Ja_Face)
  END IF

  ! Print CalcMetrics progress
  IF (doPrintStatusLine .AND. MPIRoot) THEN
    percent = REAL(iElem)/REAL(nElems)*100.
    WRITE(fmtName,*) prms%maxNameLen
    WRITE(UNIT_stdOut,'(A3,A'//ADJUSTL(fmtName)//',A2,A,A1,A,A3,F6.2,A3,A1)',ADVANCE='NO')    &
    ' | ','Calculating Metrics',' |',     &
      REPEAT('=',MAX(CEILING(percent*(prms%maxValueLen+2)/100.)-1,0)),'>',&
      REPEAT(' ',(prms%maxValueLen+2)-MAX(CEILING(percent*(prms%maxValueLen+2)/100.),0)),'| [',percent,'%] ',&
    ACHAR(13) ! ACHAR(13) is carriage return
  END IF
END DO !iElem=1,nElems

#if USE_MPI
! Send surface geomtry informations from mpi master to mpi slave
ALLOCATE(Geo(10,0:PP_N,0:PP_NZ,0:FV_SIZE,firstMPISide_MINE:nSides))
Geo=0.
Geo(1,:,:,:,:)   =SurfElem(  :,0:PP_NZ,:,firstMPISide_MINE:nSides)
Geo(2:4,:,:,:,:) =NormVec (:,:,0:PP_NZ,:,firstMPISide_MINE:nSides)
Geo(5:7,:,:,:,:) =TangVec1(:,:,0:PP_NZ,:,firstMPISide_MINE:nSides)
Geo(8:10,:,:,:,:)=TangVec2(:,:,0:PP_NZ,:,firstMPISide_MINE:nSides)
MPIRequest_Geo=MPI_REQUEST_NULL
CALL StartReceiveMPIData(Geo,10*(PP_N+1)**(PP_dim-1)*(FV_SIZE+1),firstMPISide_MINE,nSides,MPIRequest_Geo(:,RECV),SendID=1) ! Receive YOUR / Geo: master -> slave
CALL StartSendMPIData(   Geo,10*(PP_N+1)**(PP_dim-1)*(FV_SIZE+1),firstMPISide_MINE,nSides,MPIRequest_Geo(:,SEND),SendID=1) ! SEND MINE / Geo: master -> slave
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Geo)
SurfElem  (:,0:PP_NZ,:,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(1   ,:,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
NormVec (:,:,0:PP_NZ,:,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(2:4 ,:,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
TangVec1(:,:,0:PP_NZ,:,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(5:7 ,:,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
TangVec2(:,:,0:PP_NZ,:,firstMPISide_YOUR:lastMPISide_YOUR)= Geo(8:10,:,:,:,firstMPISide_YOUR:lastMPISide_YOUR)
DEALLOCATE(Geo)
#endif /*USE_MPI*/

!#if FV_ENABLED
!#if USE_MPI
!MPIRequest_Geo=MPI_REQUEST_NULL
!CALL StartReceiveMPIData(sJ_slave(:,:,:,:,0),(PP_N+1)*(PP_NZ+1),1,nSides,MPIRequest_Geo(:,SEND),SendID=2)
!CALL ProlongToFace1(PP_N,detJac_N,sJ_master(:,:,:,:,0),sJ_slave(:,:,:,:,0),L_Minus,L_Plus,doMPISides=.TRUE.,pureDG=.TRUE.)
!CALL U_Mortar1(sJ_master(:,:,:,:,0),sJ_slave(:,:,:,:,0),doMPISides=.TRUE.,pureDG=.TRUE.)
!CALL StartSendMPIData(   sJ_slave(:,:,:,:,0),(PP_N+1)*(PP_NZ+1),1,nSides,MPIRequest_Geo(:,RECV),SendID=2)
!#endif
!CALL ProlongToFace1(PP_N,detJac_N,sJ_master(:,:,:,:,0),sJ_slave(:,:,:,:,0),L_Minus,L_Plus,doMPISides=.FALSE.,pureDG=.TRUE.)
!CALL U_Mortar1(sJ_master(:,:,:,:,0),sJ_slave(:,:,:,:,0),doMPISides=.FALSE.,pureDG=.TRUE.)
!#if USE_MPI
!CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Geo)
!#endif
!sJ_slave( :,:,:,:,0) = 1./sJ_slave( :,:,:,:,0)
!sJ_master(:,:,:,:,0) = 1./sJ_master(:,:,:,:,0)
!#endif /* FV_ENABLED */

END SUBROUTINE CalcMetrics



!==================================================================================================================================
!> Prepares computation of the faces' normal, tangential vectors, surface area and Gauss points from volume metrics.
!> Input is JaCL_N, the 3D element metrics on Cebychev-Lobatto points.
!> For each side the volume metrics are interpolated to the surface and rotated into the side reference frame.
!==================================================================================================================================
SUBROUTINE CalcSurfMetrics(Nloc,FVE,JaCL_N,XCL_N,Vdm_CLN_N,iElem,NormVec,TangVec1,TangVec2,SurfElem,Face_xGP,Ja_Face)
! MODULES
USE MOD_Mathtools        ,ONLY: CROSS
USE MOD_Mesh_Vars        ,ONLY: ElemToSide,nSides,meshHasMortars,MortarType
#if PP_dim == 2
USE MOD_Mesh_Vars        ,ONLY: MortarInfo
#endif
USE MOD_Mesh_Vars        ,ONLY: NormalDirs,TangDirs,NormalSigns
USE MOD_Mappings         ,ONLY: SideToVol2
USE MOD_ChangeBasis      ,ONLY: ChangeBasis2D
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
USE MOD_Mortar_Metrics   ,ONLY: Mortar_CalcSurfMetrics
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                                             !< (IN) polynomial degree
INTEGER,INTENT(IN) :: FVE                                              !< (IN) Finite Volume enabled
INTEGER,INTENT(IN) :: iElem                                            !< (IN) element index
REAL,INTENT(IN)    :: JaCL_N(  3,3,0:Nloc,0:Nloc,0:ZDIM(Nloc))         !< (IN) volume metrics of element
REAL,INTENT(IN)    :: XCL_N(     3,0:Nloc,0:Nloc,0:ZDIM(Nloc))         !< (IN) element geo. interpolation points (CL)
REAL,INTENT(IN)    :: Vdm_CLN_N(   0:Nloc,0:Nloc)                      !< (IN) Vandermonde matrix from Cheby-Lob on N to final nodeset on N
REAL,INTENT(OUT)   ::    NormVec(3,0:Nloc,0:ZDIM(Nloc),0:FVE,1:nSides) !< (OUT) element face normal vectors
REAL,INTENT(OUT)   ::   TangVec1(3,0:Nloc,0:ZDIM(Nloc),0:FVE,1:nSides) !< (OUT) element face tangential vectors
REAL,INTENT(OUT)   ::   TangVec2(3,0:Nloc,0:ZDIM(Nloc),0:FVE,1:nSides) !< (OUT) element face tangential vectors
REAL,INTENT(OUT)   ::   SurfElem(  0:Nloc,0:ZDIM(Nloc),0:FVE,1:nSides) !< (OUT) element face surface area
REAL,INTENT(OUT)   ::   Face_xGP(3,0:Nloc,0:ZDIM(Nloc),0:FVE,1:nSides) !< (OUT) element face interpolation points
REAL,INTENT(OUT),OPTIONAL :: Ja_Face(3,3,0:Nloc,0:ZDIM(Nloc),1:nSides) !< (OUT) surface metrics
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2),dd,iLocSide,SideID,SideID2,iMortar,nbSideIDs(4),flip
#if PP_dim == 2
INTEGER            :: nMortars,tmp_MI(1:2),SideID_Mortar
#endif
INTEGER            :: NormalDir,TangDir
REAL               :: NormalSign
REAL               :: Ja_Face_l(3,3,0:Nloc,0:ZDIM(Nloc))
REAL               :: Mortar_xGP( 3,0:Nloc,0:ZDIM(Nloc),4)
REAL               :: tmp(        3,0:Nloc,0:ZDIM(Nloc))
REAL               :: tmp2(       3,0:Nloc,0:ZDIM(Nloc))
! Mortars
REAL,ALLOCATABLE   :: Mortar_Ja(:,:,:,:,:)
!==================================================================================================================================

IF (meshHasMortars) ALLOCATE(Mortar_Ja(3,3,0:Nloc,0:ZDIM(Nloc),4))

#if PP_dim == 3
DO iLocSide=1,6
#else
DO iLocSide=2,5
#endif
  flip = ElemToSide(E2S_FLIP,iLocSide,iElem)
  IF(flip.NE.0) CYCLE ! only master sides with flip=0
  SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)

  SELECT CASE(iLocSide)
  CASE(XI_MINUS)
    tmp=XCL_N(1:3,0   ,:   ,:   )
  CASE(XI_PLUS)
    tmp=XCL_N(1:3,Nloc,:   ,:   )
  CASE(ETA_MINUS)
    tmp=XCL_N(1:3,:   ,0   ,:   )
  CASE(ETA_PLUS)
    tmp=XCL_N(1:3,:   ,Nloc,:   )
  CASE(ZETA_MINUS)
    tmp=XCL_N(1:3,:   ,:   ,0   )
  CASE(ZETA_PLUS)
    tmp=XCL_N(1:3,:   ,:   ,Nloc)
  END SELECT
  CALL ChangeBasisSurf(3,Nloc,Nloc,Vdm_CLN_N,tmp,tmp2)
  ! turn into right hand system of side
  DO q=0,ZDIM(Nloc); DO p=0,Nloc
    pq=SideToVol2(Nloc,p,q,0,iLocSide,PP_dim)
    ! Compute Face_xGP for sides
    Face_xGP(1:3,p,q,0,sideID)=tmp2(:,pq(1),pq(2))
  END DO; END DO ! p,q

  Ja_Face_l=0.
  DO dd=1,PP_dim
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=JaCL_N(dd,1:3,0   ,:   ,:   )
    CASE(XI_PLUS)
      tmp=JaCL_N(dd,1:3,Nloc,:   ,:   )
    CASE(ETA_MINUS)
      tmp=JaCL_N(dd,1:3,:   ,0   ,:   )
    CASE(ETA_PLUS)
      tmp=JaCL_N(dd,1:3,:   ,Nloc,:   )
    CASE(ZETA_MINUS)
      tmp=JaCL_N(dd,1:3,:   ,:   ,0   )
    CASE(ZETA_PLUS)
      tmp=JaCL_N(dd,1:3,:   ,:   ,Nloc)
    END SELECT
    CALL ChangeBasisSurf(3,Nloc,Nloc,Vdm_CLN_N,tmp,tmp2)
    ! turn into right hand system of side
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      pq=SideToVol2(Nloc,p,q,0,iLocSide,PP_dim)
      Ja_Face_l(dd,1:3,p,q)=tmp2(:,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! dd
  IF(PRESENT(Ja_Face)) Ja_Face(:,:,:,:,SideID)=Ja_Face_l


  NormalDir=NormalDirs(iLocSide); TangDir=TangDirs(iLocSide); NormalSign=NormalSigns(iLocSide)
  CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face_l,&
                         NormVec(:,:,:,0,SideID),TangVec1(:,:,:,0,SideID),&
                         TangVec2(:,:,:,0,SideID),SurfElem(:,:,0,SideID))


#if PP_dim == 2
  IF (iLocSide.EQ.XI_MINUS) THEN
    nMortars=MERGE(2,0,MortarType(1,sideID).EQ.2 .OR. MortarType(1,sideID).EQ.3)
    IF (nMortars.EQ.2) THEN
      SideID_Mortar = MortarType(2,sideID)
      tmp_MI    = MortarInfo(:,1,SideID_Mortar)
      MortarInfo(:,1,SideID_Mortar) = MortarInfo(:,2,SideID_Mortar)
      MortarInfo(:,2,SideID_Mortar) = tmp_MI
    END IF
  END IF
#endif

  !compute metrics for mortar faces, interpolate Ja_Face to small sides
  IF(MortarType(1,SideID).GT.0)THEN
    CALL Mortar_CalcSurfMetrics(SideID,Nloc,Ja_Face_l,Face_xGP(:,:,:,0,SideID),&
                                            Mortar_Ja,Mortar_xGP,nbSideIDs)
    DO iMortar=1,4
      SideID2=nbSideIDs(iMortar)
      IF(SideID2.LT.1) CYCLE ! for MPI sides some sides are built from the inside and for type 2/3 there are only 2 neighbours
      IF(PRESENT(Ja_Face)) Ja_Face(:,:,:,:,SideID2)=Mortar_Ja(:,:,:,:,iMortar)
      Face_xGP(:,:,:,0,SideID2) = Mortar_xGP(:,:,:,iMortar)
      CALL SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Mortar_Ja(:,:,:,:,iMortar),&
                             NormVec(:,:,:,0,SideID2),TangVec1(:,:,:,0,SideID2),&
                             TangVec2(:,:,:,0,SideID2),SurfElem(:,:,0,SideID2))
    END DO
  END IF
END DO

IF (meshHasMortars) DEALLOCATE(Mortar_Ja)

END SUBROUTINE CalcSurfMetrics

!==================================================================================================================================
!> Computes surface normal and tangential vectors and surface area from surface metrics Ja_Face.
!==================================================================================================================================
SUBROUTINE SurfMetricsFromJa(Nloc,NormalDir,TangDir,NormalSign,Ja_Face,NormVec,TangVec1,TangVec2,SurfElem)
! MODULES
USE MOD_Mathtools,   ONLY: CROSS
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nloc                             !< polynomial degree
INTEGER,INTENT(IN) :: NormalDir                        !< direction of normal vector
INTEGER,INTENT(IN) :: TangDir                          !< direction of 1. tangential vector
REAL,INTENT(IN)    :: NormalSign                       !< sign of normal vector
REAL,INTENT(IN)    :: Ja_Face(3,3,0:Nloc,0:ZDIM(Nloc)) !< face metrics
REAL,INTENT(OUT)   ::   NormVec(3,0:Nloc,0:ZDIM(Nloc)) !< element face normal vectors
REAL,INTENT(OUT)   ::  TangVec1(3,0:Nloc,0:ZDIM(Nloc)) !< element face tangential vectors
REAL,INTENT(OUT)   ::  TangVec2(3,0:Nloc,0:ZDIM(Nloc)) !< element face tangential vectors
REAL,INTENT(OUT)   ::  SurfElem(  0:Nloc,0:ZDIM(Nloc)) !< element face surface area
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q
!==================================================================================================================================
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  ! norm of the non-normalized physical normal vector
  ! |n| = |J dX^-T n_ref| = |J (a^i)^T n_ref|
  SurfElem(  p,q) = SQRT(SUM(Ja_Face(NormalDir,1:PP_dim,p,q)**2))
  NormVec( :,p,q) = NormalSign*Ja_Face(NormalDir,:,p,q)/SurfElem(p,q)
  ! For two-dimensional computations, the normal direction will be 1 or 2. For the tangential direction
  ! we then set 2 or 1 accordingly.
  TangVec1(:,p,q) = Ja_Face(TangDir,:,p,q) - SUM(Ja_Face(TangDir,:,p,q)*NormVec(:,p,q)) &
                    *NormVec(:,p,q)
  TangVec1(:,p,q) = TangVec1(:,p,q)/SQRT(SUM(TangVec1(:,p,q)**2))
#if (PP_dim == 2)
  TangVec2(:,p,q) = 0.
#else
  TangVec2(:,p,q) = CROSS(NormVec(:,p,q),TangVec1(:,p,q))
#endif
END DO; END DO ! p,q
END SUBROUTINE SurfMetricsFromJa

END MODULE MOD_Metrics
