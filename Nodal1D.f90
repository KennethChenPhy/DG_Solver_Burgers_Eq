! DG method 1D element nodes and modes
! Jacobi Polynomials
!==========================================================
        module Nodal1D
                implicit none
                contains
! Normalized Jacobi Polynomials
! JacobiP(x,alpha,beta,N,Nc) x /in [-1,1]
! Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1 and alpha+beta/=-1
!    at points x(1:Nc) for order N and returns P(0:N,1:Nc) of values
! Note   : They are normalized to be orthonormal.
         subroutine JacobiP_(x, alpha, beta, N,Nc, P)
            implicit none
            real(8),intent(in) :: x(:),alpha,beta
            integer,intent(in) :: N,Nc
            real(8) :: P(0:N,1:Nc)

            real(8) :: P0,P1,a,a1,b
            integer :: i,j
!initial values of P_0(x) and P_1(x)
      P0=sqrt(2.d0**(-alpha-beta-1)*gamma(alpha+beta+2)&
              /gamma(alpha+1)/gamma(beta+1))
      P1=0.5*P0*sqrt((alpha+beta+3)/((alpha+1)*(beta+1)))
      if (N==0) then
         do i=1,Nc
                 P(0,i)=P0
         enddo
 elseif(N==1)then
                do i=1,Nc
                P(0,i)=P0
                P(1,i)=P1*((alpha+beta+2)*x(i)+(alpha-beta))
                enddo
        else
!recurrence

                do i=1,Nc
                P(0,i)=P0
                P(1,i)=P1*((alpha+beta+2)*x(i)+(alpha-beta))
                enddo

        a=2.d0/(2+alpha+beta)*sqrt((1+alpha)*(1+beta)/(3+alpha+beta))
        do j=2,N
        a1=2.d0/(2.d0*j+alpha+beta)*sqrt(j*(j+alpha+beta)*(j+alpha)&
        *(j+beta)/((2.d0*j+alpha+beta-1)*(2.d0*j+alpha+beta+1)))
        b=-(alpha**2-beta**2)/(2.d0*(j-1)+alpha+beta)/(2.d0*(j-1)+alpha+beta+2)
                do i=1,Nc
                P(j,i)=(x(i)*P(j-1,i)-a*P(j-2,i)-b*P(j-1,i))/a1
                enddo
                a=a1
        enddo
        endif


endsubroutine JacobiP_
!Gradient of Jacobi Polynomials
!Gives only upto N-1 th order dP_n=P_{n-1}
        subroutine GradJacobiP_(x,alpha,beta,N,Nc,dP)
                            implicit none
            integer,intent(in) :: N,Nc
            real(8),intent(in) :: x(1:Nc),alpha,beta
            real(8) :: dP(0:N,1:Nc),dPtemp(0:N-1,1:Nc)
            integer::i,j
            if(N==0)then
                    do i=1,Nc
                    dP(0,i)=0.d0
                    enddo
            else
                    do i=1,Nc
                    dP(0,i)=0.d0
                    enddo
                    call JacobiP_(x,alpha+1,beta+1,N-1,Nc,dPtemp)
                    do i=1,N
                    do j=1,Nc
                    dP(i,j)=sqrt(i*(i+alpha+beta+1))*dPtemp(i-1,j)
                    enddo
                    enddo
            endif
        endsubroutine GradJacobiP_
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function JacobiGQ(alpha,beta,N)
!Purpose: Compute the N'th order Gaussian quadrature points, x and wights, w
!         associated with the Jacobi polynomial of type (alpha,beta)
! alpha,beta>-1 and alpha+beta/=-1  P(row 1): x  ;   P(row 2): w
! USE Lapack for eigensystem of tridiagonal Jacobi matrix
! CHECK LAPACK STEQR for details
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine JacobiGQ_(alpha,beta,N,P)
    implicit none
        real(8),intent(in) :: alpha,beta
        integer,intent(in) :: N
        real(8) :: P(0:N,1:2)
        ! diagonal elemetents of Jacobi Matrix d, output ascending eigenvalues
        ! off-diagonal elements
        ! orthonormal eigenvectors i-th column corresponds to i-th eigenvalue
        ! info=0 means successful execution
        real(8) :: d(0:N),e(1:N),z(0:N,0:N),work(2*N)
        integer :: info,ldz,j
        ldz=N+1
        if (N==0) then
                P(0,1)=(beta-alpha)/(alpha+beta+2.0)
                P(0,2)=2.0
        else
                do j=1,N
                d(j)=-(alpha**2-beta**2)/(2.d0*j+alpha+beta)/(2.d0*j+alpha+beta+2)
                enddo
                do j=1,N
                e(j)=2.d0/(2.d0*j+alpha+beta)*sqrt(j*(j+alpha+beta)*(j+alpha)*(j+beta)&
                     /((2.d0*j+alpha+beta-1)*(2.d0*j+alpha+beta+1)))
                enddo

                d(0)=(beta-alpha)/(alpha+beta+2.0)

        call dsteqr('I',N+1,d,e,z,ldz,work,info)
              if (info==0) then
!                write(*,*) 'STEQR succesful execuation'
                do j=0,N
                P(j,1)=d(j)
                P(j,2)=z(0,j)**2*2.d0**(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)&
                        /gamma(alpha+beta+2)
                enddo
               else
                write (*,*) 'err: fail to execuate LAPACK:STEQR'
              endif

        endif
  endsubroutine JacobiGQ_

!Jacobi Gauss-Lobatto points
!JacobiGQ(alpha+1,beta+1,N-2)
!Only output nodes but not weights
! see module LegendreGL for both nodes and weights
  subroutine JacobiGL_(alpha,beta,N,P)
     implicit none
         real(8),intent(in) :: alpha,beta
         integer,intent(in) :: N
         real(8) :: P(0:N),Ptemp(0:N-2,1:2)
         integer :: i
         if (N==1) then
             P(0)=-1.0
             P(1)=1.0
         else
             call JacobiGQ_(alpha+1,beta+1,N-2,Ptemp)
             P(0)=-1.0
             do i=1,N-1
                 P(i)=Ptemp(i-1,1)
             enddo
             P(N)=1.0
         endif
  end subroutine JacobiGL_

!Jegendre Gauss-Lobatto nodes and weights
! see module JacobiGL for only nodes but not weights

  subroutine LegendreGL_(N,P)
     implicit none

         integer,intent(in) :: N
         real(8) :: P(0:N,1:2),Ptemp(0:N-2,1:2),Ptemp2(0:N,1:N-1)
         integer :: i
         real(8) :: alpha,beta
         alpha=0.0
         beta=0.0
         if (N==1) then
             P(0,1)=-1.0
             P(1,1)=1.0
             P(0,2)=1.0
             P(1,2)=1.0
         else
             call JacobiGQ_(alpha+1,beta+1,N-2,Ptemp)
             P(0,1)=-1.0
             P(0,2)=2.d0/N/(N+1)
             call JacobiP_(Ptemp(:,1),0.d0,0.d0,N,N-1,Ptemp2)
             do i=1,N-1
                 P(i,1)=Ptemp(i-1,1)
                 P(i,2)=(2.d0*N+1.d0)/N/(N+1.d0)/(Ptemp2(N,i))**2
!Normalized Legendre Polynomial
             enddo
             P(N,1)=1.0
             P(N,2)=2.d0/N/(N+1)
         endif
  end subroutine LegendreGL_

! Vandermonde Matrix to transfer between nodal and modal representations
! Vandermonde1D(N,r,P): V_ij=P_j(r_i)
!N: order of polynomials; r: vector of nodes size of Nc
!V_ij \hat{u}_j = u_i    \hat{u}_i are the modal, u_i are the nodal
        subroutine Vandermonde1D_(N,Nc,r,P)
                implicit none
        integer,intent(in)::N,Nc
        real(8),intent(in)::r(1:Nc)
        real(8)::P(1:Nc,0:N),Ptemp(0:N,1:Nc)
        call JacobiP_(r,0.d0,0.d0,N,Nc,Ptemp)
        P=transpose(Ptemp)
        endsubroutine Vandermonde1D_
! Gradient of Vandermonde Matrix V_r,ij=dP_j/dr at r_i
        subroutine GradVandermonde1D_(N,Nc,r,DVr)
                implicit none
                integer,intent(in)::N,Nc
                real(8),intent(in)::r(1:Nc)
                real(8)::DVr(1:Nc,0:N),DVtemp(0:N,1:Nc)
                call GradJacobiP_(r,0.d0,0.d0,N,Nc,DVtemp)
                DVr=transpose(DVtemp)
        endsubroutine GradVandermonde1D_
!Initialize the differentiation matrix D_r=V_r V^-1 on the interval

        subroutine Dmatrix1D_(N,Nc,r,V)
                use Global_math
                implicit none
                integer,intent(in)::N,Nc
                real(8),intent(in)::r(1:Nc)
                real(8)::V(1:Nc,1:Nc),Vtemp(1:Nc,0:N),Vrtemp(1:Nc,0:N)

                integer::info
                if(Nc==N+1)then
                call Vandermonde1D_(N,Nc,r,Vtemp)
                call GradVandermonde1D_(N,Nc,r,Vrtemp)
                call InverseMatrix(Nc,Vtemp,info)
                if(info==0)then
                ! return Vtemp as Nc*Nc size inverse matrix of the original

                V=matmul(Vrtemp,Vtemp)
                else
                write(*,*) 'err: fail to invert Vandermonde matrix'
        endif
        else
                write(*,*) 'err: Modal and Nodal sizes do not match N_c\=N+1'
        endif
        endsubroutine Dmatrix1D_
! Lift1D function output Nc*2 array for surface term
! Input Nc*Nc Vandermonde Matrix
! Output inverse(mass matrix)*(edge unit vector)
        function Lift1D_(Nc,V) result(P)
                implicit none
                integer,intent(in)::Nc
                real(8),intent(in)::V(1:Nc,1:Nc)
                real(8)::P(1:Nc,1:2),emat(Nc,2)
                integer::i,j
                emat=0.d0
                emat(1,1)=1.d0
                emat(Nc,2)=1.d0
                P=matmul(matmul(V,transpose(V)),emat)
!                do i=1,Nc
!                    do j=1,2
!                        write(*,'(F10.4)',advance='no') P(i,j)
!                            enddo
!                    write(*,*)
!                enddo
!                do i=1,Nc
!                        do j=1,Nc
!                                P(i,1)=V(i,j)*V(1,j)
!                                P(i,2)=V(i,j)*V(Nc,j)
!                        enddo
!                enddo
!                do i=1,Nc
!                    do j=1,2
!                        write(*,'(F10.4)',advance='no') P(i,j)
!                            enddo
!                    write(*,*)
!                enddo

        end function Lift1D_

        endmodule Nodal1D
