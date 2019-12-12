module Global_math
        implicit none
real(8),parameter::rk4a(5)=(/0.d0,&
       -567301805773.d0/1357537059087.d0,&
       -2404267990393.d0/2016746695238.d0,&
       -3550918686646.d0/2091501179385.d0,&
       -1275806237668.d0/842570457699.d0/)
real(8),parameter::rk4b(5)= (/1432997174477.d0/9575080441755.d0,&
         5161836677717.d0/13612068292357.d0,&
         1720146321549.d0/2090206949498.d0,&
         3134564353537.d0/4481467310338.d0,&
         2277821191437.d0/14882151754819.d0/)
real(8),parameter::rk4c(5)= (/0.d0,&
         1432997174477.d0/9575080441755.d0,&
         2526269341429.d0/6820363962896.d0,&
         2006345519317.d0/3224310063776.d0,&
         2802321613138.d0/2924317926251.d0/)

      contains



        subroutine InverseMatrix(N,M,info)
!Call LAPACK to INVERSE a N by N Matrix
        implicit none
        integer,intent(in)::N
        real(8)::M(1:N,1:N),work(N)
        integer:: info,ipiv(N)
        call dgetrf(N,N,M,N,ipiv,info)
        if(info==0)then
        call dgetri(N,M,N,ipiv,work,N,info)
        endif
        end subroutine InverseMatrix

!************************************************************************************
! One Step Classic Runge Kutta 4th order method
!************************************************************************************

SUBROUTINE RK4(IFcn, n, t, h, xi, y)
IMPLICIT NONE
!------------------------------------------------------------------
!       Solving system of ODE by RK4 method
!       input:
!       IFcn    - subroutine of the input fcn
!       n               - number of ODEs
!       t              - initial time
!       h               -time step
!       xi()    - initial values
!       output:
!       y()             - output value
!------------------------------------------------------------------
INTEGER :: i, j, n
DOUBLE PRECISION, DIMENSION(1:n) :: K1, K2, K3, K4, xi, x, y, fcn
REAL(8) :: h, t
EXTERNAL IFcn


        CALL IFcn(n, t, xi, fcn)
        DO j = 1, n
                K1(j) = fcn(j) *(h)
                x(j) = xi(j) + K1(j)/2.D0
        ENDDO

        CALL IFcn(n, t + h/2.D0, x, fcn)
        DO j = 1, n
                K2(j) = fcn(j) *(h)
                x(j) = xi(j) + K2(j)/2.D0
        ENDDO

        CALL IFcn(n, t + h/2.D0, x, fcn)
        DO j = 1, n
                K3(j) = fcn(j) *(h)
                x(j) = xi(j) + K3(j)
        ENDDO

        CALL IFcn(n, t + h, x, fcn)
        DO j = 1, n
                K4(j) = fcn(j) *(h)
                y(j) = xi(j) + (K1(j) + 2.D0*K2(j) + 2.D0*K3(j) + K4(j))/6.D0
        ENDDO
ENDSUBROUTINE RK4

!************************************************************************************
! One Step Low-Storage Runge Kutta 4th order method
!************************************************************************************

SUBROUTINE LSRK4(IFcn, n, t, h, xi, y)
IMPLICIT NONE
!------------------------------------------------------------------
!       Solving system of ODE by LSRK4 method
!       input:
!       IFcn    - subroutine of the input fcn
!       n               - number of ODEs
!       t              - initial time
!       h              -  time step
!       xi()    - initial values
!       output:
!       y()             - output value
!------------------------------------------------------------------
INTEGER :: i, j, n
DOUBLE PRECISION, DIMENSION(1:n) :: k, xi, y, fcn
REAL(8) :: h, t
EXTERNAL IFcn
        y=xi
        do i=1,5
        CALL IFcn(n, t+rk4c(i)*h, y, fcn)
        DO j = 1, n
        k(j)=rk4a(i)*k(j)+h*fcn(j)
        ENDDO
        do j=1,n
        y(j)=y(j)+rk4b(i)*k(j)
        enddo
        enddo
ENDSUBROUTINE LSRK4

SUBROUTINE Write2R8(Un, F,E1, E2)
implicit none
REAL(8) :: E1, E2
integer :: Un
character(len=*) :: F
    write(Un, F, advance = "no") E1
    write(Un, F) E2
ENDSUBROUTINE Write2R8


      end module Global_math
