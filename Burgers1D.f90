!Right Hand Side of Burgers' equation
module Burgers1D
  use Global_var
  use Global_math
  implicit none

contains
  function BurgersRHS1D(uwork,time) result(RHS)
    !initial value u of size Nc*K
    !use EToE, EToF
    !page 141 or 154 in Book Jan S. Hesthaven & Tim Warburton
    !f(u)=u^2
    use Global_var
    implicit none
    real(8),intent(in)::uwork(1:Nc,1:K),time
    real(8)::RHS(Nc,K)
    real(8)::C,fast(2,K),du(2,K),du2(2,K)
    integer::row,col,i,j
    !Set boundary condition
    if (cL-3*time<=-0.5) then
      bcL=2.d0
    else
      bcL=1.d0
    endif
    if (cR-3*time<=-0.5) then
      bcR=2.d0
    else
      bcR=1.d0
    endif
    !Lax-Friedrichs Flux
    C=2*maxval(abs(uwork))

  do i=1,K
      do j=1,2
      col=EToE(i,j) !id of exterior element
      row=int(EToF(i,j)/2)*(Nc-1)+1 !id of exterior element's face
      du(j,i)=uwork(int(j/2)*(Nc-1)+1,i)-uwork(row,col)
      du2(j,i)=uwork(int(j/2)*(Nc-1)+1,i)**2-uwork(row,col)**2
    enddo
    enddo
    du(1,1)=2.d0*(uwork(1,1)-bcL)
    du(2,K)=2.d0*(uwork(Nc,K)-bcR)
    du2(1,1)=2.d0*(uwork(1,1)**2-bcL**2)
    du2(2,K)=2.d0*(uwork(Nc,K)**2-bcR**2)
    fast=du2/2.d0-nx*C/2.d0*du
    RHS=-rx*matmul(Dr,(uwork*uwork))+matmul(Lift,(Fscale*nx*fast))

  endfunction BurgersRHS1D
!============================================================================

  subroutine Burgers1D_(u,ftime)
    use Global_var
    use Global_math
    use SlopeLimiter
    implicit none
    real(8)::ftime,u(Nc,K)
    real(8)::time,C,CFL,dt,xmin,work(Nc,K)
    real(8)::v1(Nc,K),v2(Nc,K),v3(Nc,K),v4(Nc,K)
    integer::Nsteps,i,j
  integer::ti,tj
    !set inital time step
    time=0.d0
    CFL=0.25
    xmin=minval(abs(Fx(1,1:K)-Fx(2,1:K)))
    C=2.d0*maxval(abs(u))
    dt=CFL*xmin/C
    Nsteps=int(ftime/dt)
    dt=ftime/Nsteps
if(opt2==2)then
call SlopeLimitN(u)
else 
        call SlopeLimit1(u)
endif
    do j=1,Nsteps
!Classic Explicit RK4 method
    if(opt==1)then
        v1=BurgersRHS1D(u,time)*dt
        work=u+v1/2.d0
if(opt2==2)then
call SlopeLimitN(work)
else 
        call SlopeLimit1(work)
endif
        v2=BurgersRHS1D(work,time+dt/2.d0)*dt
        work=u+v2/2.d0
if(opt2==2)then
call SlopeLimitN(work)
else 
        call SlopeLimit1(work)
endif
        v3=BurgersRHS1D(work,time+dt/2.d0)*dt
        work=u+v3
if(opt2==2)then
call SlopeLimitN(work)
else 
        call SlopeLimit1(work)
endif
        v4=BurgersRHS1D(work,time+dt)*dt
        u=u+(v1+2.d0*v2+2.d0*v3+v4)/6.d0
if(opt2==2)then
call SlopeLimitN(u)
else 
        call SlopeLimit1(u)
endif
    elseif(opt==2)then
!LS-RK4 Method:
      do i=1,5
        work=rk4a(i)*work+dt*BurgersRHS1D(u,time+rk4c(i)*dt)
        u=u+rk4b(i)*work

if(opt2==2)then
call SlopeLimitN(u)
else 
        call SlopeLimit1(u)
endif

      enddo
      elseif(opt==3) then
!  write(*,*) time
!SSP-RK4:
v1=u+0.39175222700392*dt*BurgersRHS1D(u,time)
if(opt2==2)then
call SlopeLimitN(v1)
else 
        call SlopeLimit1(v1)
endif
v2=0.44437049406734*u+0.55562950593266*v1&
  +0.36841059262959*dt*BurgersRHS1D(v1,time+0.39175222700392*dt)
if(opt2==2)then
call SlopeLimitN(v2)
else 
        call SlopeLimit1(v2)
endif
v3=0.62010185138540*u+0.37989814861460*v2&
  +0.25189177424738*dt*BurgersRHS1D(v2,time+0.58607968896780*dt)
if(opt2==2)then
call SlopeLimitN(v3)
else 
        call SlopeLimit1(v3)
endif
v4=0.17807995410773*u+0.82192004589227*v3&
  +0.54497475021237*dt*BurgersRHS1D(v3,time+0.47454236302687*dt)
if(opt2==2)then
call SlopeLimitN(v4)
else 
        call SlopeLimit1(v4)
endif
u=0.00683325884039*u+0.51723167208978*v2&
  +0.12759831133288*v3+0.34833675773694*v4&
  +0.08460416338212*dt*BurgersRHS1D(v3,time+0.47454236302687*dt)&
  +0.22600748319395*dt*BurgersRHS1D(v4,time+0.93501063100924*dt)
if(opt2==2)then
call SlopeLimitN(u)
else 
        call SlopeLimit1(u)
endif

else
write(*,*) 'err: try again.'
exit
endif
      time=time+dt
    enddo
  endsubroutine Burgers1D_



endmodule Burgers1D
