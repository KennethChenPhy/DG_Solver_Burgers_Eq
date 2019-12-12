module SlopeLimiter
        implicit none
contains
function minmod(v) result(a)
        implicit none
        real(8)::v(:),a
        integer::n,s,i
        a=minval(abs(v))
        n=size(v)
        s=0
        do i=1,n
        s=s+int(sign(1.d0+1d-10,v(i)))
        enddo
        s=int(s/n)
        if (abs(s)==1) then
                a=s*a
        else
                a=0.d0
        endif
endfunction minmod

!Scope Limiter

function SlopeLimit1Element(uk,cord,uave) result(PIuk)
        use Global_var
        implicit none
        !for element u^k with coordinates
        !ukp average of u^k+1
        real(8),intent(in)::uk(Nc),cord(Nc),uave(-1:1)
        real(8)::x0,PIuk(Nc),h,uv(1:3)
        integer::i

        h=cord(Nc)-cord(1)
        x0=cord(1)+h/2.d0
        PIuk=2.d0/h*matmul(Dr,uk)
        uv(2)=(uave(1)-uave(0))/h*2.d0
        uv(3)=(uave(0)-uave(-1))/h*2.d0

        do i=1,Nc

        uv(1)=PIuk(i)
        PIuk(i)=uave(0)+(cord(i)-x0)*minmod(uv)

        enddo
endfunction SlopeLimit1Element

subroutine SlopeLimit1(u)
        use Global_var
        implicit none
        real(8)::u(Nc,K),utilde(0:N,K),uave(-1:1,K),uavg(Nc,K)
        integer::i
utilde=matmul(invV,u)
utilde(1:N,:)=0.d0
uavg=matmul(V,utilde)
uave(0,:)=uavg(1,:)
do i=1,K
uave(-1,i)=uave(0,EToE(i,1))
uave(1,i)=uave(0,EToE(i,2))
enddo
utilde=matmul(invV,u)
utilde(2:N,:)=0.d0
uavg=matmul(V,utilde)

do i=1,K
u(:,i)=SlopeLimit1Element(uavg(:,i),x(:,i),uave(:,i))
enddo

endsubroutine SlopeLimit1

subroutine SlopeLimitN(u)
        use Global_var
        implicit none
        real(8)::u(Nc,K),uavg(Nc,K)
        real(8)::uv(2,K),uave(-1:1,K),m(3),x0(K),utilde(0:N,K)
        integer::i,j
!find average of each element and averages of connected elements
utilde=matmul(invV,u)
utilde(1:N,:)=0.d0
uavg=matmul(V,utilde)
uave(0,:)=uavg(1,:)
do i=1,K
uave(-1,i)=uave(0,EToE(i,1))
uave(1,i)=uave(0,EToE(i,2))
enddo

!interface flux
do i=1,K
m(1)=uave(0,i)-u(1,i)
m(2)=uave(0,i)-uave(-1,i)
m(3)=uave(1,i)-uave(0,i)
uv(1,i)=uave(0,i)-minmod(m)
m(1)=u(Nc,i)-uave(0,i)
m(2)=uave(0,i)-uave(-1,i)
m(3)=uave(1,i)-uave(0,i)
uv(2,i)=uave(0,i)+minmod(m)
enddo
do i=1,K
x0(i)=(x(Nc,i)-x(1,i))/2.d0+x(1,i)
enddo
do i=1,K
if (abs(uv(1,i)-u(1,i))<nodetol .and. abs(uv(2,i)-u(Nc,i)<nodetol)) cycle
        utilde=matmul(invV,u)
        utilde(2:N,:)=0.d0
        uavg=matmul(V,utilde)
!uavg=matmul(Dr,u)*rx
  do j=1,Nc
  uavg(j,i)=uave(0,i)+(x(j,i)-x0(i))*uavg(j,i)
  enddo
u(:,i)=SlopeLimit1Element(uavg(:,i),x(:,i),uave(:,i))

enddo
endsubroutine SlopeLimitN
endmodule SlopeLimiter
