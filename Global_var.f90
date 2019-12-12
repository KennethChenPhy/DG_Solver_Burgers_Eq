module Global_var
  implicit none
  integer::info,opt,opt2
!Jacobi Polynomial class (alpha,beta)
real(8):: alpha=0.0
real(8):: beta=0.0
real(8)::nodetol=1.d-10 !nodal tolorence
!K # of elements
!N order of polynomial approximatiom
!Nc # of nodes, usually Nc=N+1
!Nv # of vertex, usually Nv=K+1i
!Nfp # of nodes on face
!Only changable after commands: make clean and make
real(8),parameter:: Pi=3.14159265358979323846
integer,parameter:: K=100
integer,parameter:: N=20
integer,parameter:: Nc=N+1
integer,parameter:: Nv=K+1
integer,parameter:: Nfp=1
!1D end-point coordinates: xL, xR
real(8):: cL=-1.d0
real(8):: cR=1.d0
real(8):: r(1:Nc,1:2) !LGL nodes and weights
real(8):: V(1:Nc,0:N) !Vandermonde matrix
real(8):: invV(0:N,1:Nc) !Inverse of  V matrix
real(8):: Dr(1:Nc,1:Nc) !Differentiation matrix
real(8):: Lift(1:Nc,1:2) !Suface terms
real(8):: VX(0:K) !Vertices coordinates
integer:: EToV(1:K,1:2) !Element to Vertex mapping
!integer:: FToV(1:2*K,0:K)    !2K Faces mapping to K+1 Vertices
!integer:: FToF(1:2*K,1:2*K)  !Face to face mapping
integer:: EToE(1:K,1:2)
!ith element's jth face connects $(i,j) element
!if $(i,j)=i then it is a boundary
integer:: EToF(1:K,1:2)
!ith element's jst face connects EToE(i,j)th element's $EToF(i,j) face
!if $(i,j)=j then it is a boundary
real(8):: x(1:Nc,1:K) !Coordinates of all nodes
!real(8):: u(Nc,K) !solution
real(8):: bcL,bcR !Boundary Condition
real(8):: xr(1:Nc,1:K) !Local Jacobian of x: dx_i/dr_j = Dr*x
real(8):: rx(1:Nc,1:K) !element-wise inverse of J: rx_ij=1/J_ij
integer::Fmask(2)
real(8):: Fx(1:2,1:K) !Coordinates of the edge nodes
real(8):: nx(1:2,1:K)     !Outward pointing normals at element faces
real(8):: Fscale(1:2,1:K)   !Inverse metric at surface 1./J
integer::vmapM(2*K),vmapP(2*K),mapB(2),vmapB(2)
!Assume the boundary is at face 1 of element 1 and face 2 of element K:
integer::mapI=1,mapO=K*2,vmapI=1,vmapO=K*Nc



endmodule Global_var
