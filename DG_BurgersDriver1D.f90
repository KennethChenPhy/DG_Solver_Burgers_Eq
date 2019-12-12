program main
        use Global_var
        use Global_math
        use Grid1D
        use SlopeLimiter
        use Burgers1D
  implicit none
real(8)::ftime,u(Nc,K),RHS(Nc,K)
integer::i,j
    call Initial1D

do i=1,Nc
    do j=1,K
    if (x(i,j)<=-0.5) then
    u(i,j)=2.d0
    else
    u(i,j)=1.d0
    endif
enddo
enddo



!u=sin(x)

write(*,*) 'Input final time'
read(*,*) ftime
write(*,*) 'Option 1: SlopeLimit1'
write(*,*) 'Option 2: SlopeLimitN'
read(*,*) opt2
write(*,*) 'Option 1: ERK4 method'
write(*,*) 'Option 2: LS-RK4 method'
write(*,*) 'Option 3: SSP-RK4 method'
read(*,*) opt

call Burgers1D_(u,ftime)


open(201,file='sln.txt',status='replace')
do j=1,K
  do i=1,Nc
    call Write2R8(201,'(F15.8)',x(i,j),u(i,j))
enddo
enddo
close(201)


endprogram main
