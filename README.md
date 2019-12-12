good luck!


ref: https://github.com/tcew/nodal-dg
Follow the book 'Nodal Discontinuous Galerkin Method' to have a standalone Burgers' Equation Solver
1. Jacobi Polynomial: JacobiP(r,alpha,beta,n)
2. Weights and nodes for Gaussian quadrature: JacobiGQ
3. Legendre-Gauss-Lobatto (LGL) quadrature points: JacobiGL
4. Modal-Nodal transformation, the Vandermonde matrix: Vandermonde1D
5. Gradient of Jacobi Polynomial: GradJacobiP
6. Gradient of Vandermonde at grid points: GradVandermonde1D
7. differentiation matrix (MD_r=S): Dmatrix1D
8. Extract the surface terms: Lift1D
9. End points of elements as vertex coordinates of size N_v: VX
10. EToV
11. Map to the physical coordinates using affine mapping: GeometricFactors1D
12. Local outward pointing normals and inverse Jacobian at the interfaces: Normals1D
13. FToF, EToE, EToF to connect vertices: Connect1D
14. BuildMap1D
15. All in StartUp1D
16. parameters: Globals1D
17. Time evolution, low-storage five-stage fourth-order explicit RK method
18. Driver routine to read grid information
19. MeshGen1D

module Nodal1D {
subroutines:
JacobiP_(x, alpha, beta, N,Nc, P)
  :x(1:Nc),P(0:N,1:Nc)
GradJacobiP_(x,alpha,beta,N,Nc,dP)
  :x(1:Nc),dP(0:N,1:Nc)
JacobiGQ_(alpha,beta,N,P)
  :P(0:N,1:2),1 for nodes, 2 for weights
JacobiGL_(alpha,beta,N,P)
  :P(0:N)
LegendreGL_(N,P)
  :P(0:N,1:2),1 for nodes, 2 for weights
Vandermonde1D_(N,Nc,r,P)
  :r(1:Nc),P(1:Nc,0:N)
GradVandermonde1D_(N,Nc,r,DVr)
  :r(1:Nc),DVr(1:Nc,0:N)
Dmatrix1D_(N,Nc,r,V)
  : Nc=N+1
  : r(1:Nc),V(1:Nc,1:Nc)

Functions:
Lift1D_(Nc,V) output P(Nc,2)
  : V(1:Nc,1:Nc)
}
module Grid1D{

}
