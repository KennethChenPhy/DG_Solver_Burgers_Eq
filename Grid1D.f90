!With 1D element in place, this module assembles 1D grid.
!E: element; F: face; V: vertex

      module Grid1D
      implicit none
      contains
      subroutine Uniform1D_(xL,xR,K,VX,EToV)
!initialize default uniform 1D grid of K elements
              implicit none
              real(8)::xL,xR
              integer::K,i
              real(8)::VX(0:K),h
              integer::EToV(1:K,1:2)
              h=(xR-xL)/K
              do i=0,K
                VX(i)=xL+h*i
              enddo
              do i=1,K
                EToV(i,1)=i-1
                EToV(i,2)=i
              enddo
      endsubroutine Uniform1D_
      subroutine NodalCoord_(K,Nc,r,VX,EToV,x)
              implicit none
              integer,intent(in)::K,Nc,EToV(1:K,1:2)
              real(8),intent(in)::r(1:Nc),VX(0:K)
              real(8)::x(1:Nc,1:K)
              integer::i,j
              do i=1,Nc
              do j=1,K
         x(i,j)=VX(EToV(j,1))+0.5*(r(i)+1)*(VX(EToV(j,2))-VX(EToV(j,1)))
              enddo
              enddo
      endsubroutine NodalCoord_

      subroutine Connect1D_(K,EToV,EToE,EToF)
     !Build global connectivity arrays from EtoV input
              implicit none
              integer,intent(in)::K,EToV(1:K,1:2)
              integer::EToE(1:K,1:2),EToF(1:K,1:2)
              integer::f2f(2*K-2,2),ele1(2*K-2,2),ele2(2*K-2,2)
        integer::FToV(1:2*K,0:K),FToF(1:2*K,1:2*K),Ide(1:2*K,1:2*K)
        integer::i,j,p
              !Global #face to #face mapping
              !Local # face 1: left end
              !Local # face 2: right end
                FToV=0
                FToF=0
                Ide=0

              do i=1,K
              do j=1,2
              FToV((i-1)*2+j,EToV(i,j))=1
              enddo
              enddo


              do i=1,2*K
              Ide(i,i)=1
              enddo
              FToF=matmul(FToV,transpose(FToV))-Ide
              p=1
              do i=1,2*K
              do j=1,2*K
              IF (FToF(i,j)==1) then
                if(p>2*K-2)then
                        exit
                endif
                      f2f(p,1)=i
                      f2f(p,2)=j
                      p=p+1
              endif
              enddo
              enddo
              do i=1,2*K-2
                ele1(i,1)=int((f2f(i,1)-1)/2)+1
                ele1(i,2)=mod(f2f(i,1)-1,2)+1
                ele2(i,1)=int((f2f(i,2)-1)/2)+1
                ele2(i,2)=mod(f2f(i,2)-1,2)+1
              enddo
              do i=1,K
              EToE(i,1)=i
              EToE(i,2)=i
              EToF(i,1)=1
              EToF(i,2)=2
              enddo
              do i=1,2*K-2
              EToE(ele1(i,1),ele1(i,2))=ele2(i,1)
              EToF(ele1(i,1),ele1(i,2))=ele2(i,2)
              enddo
      endsubroutine Connect1D_
      subroutine BuildMaps1D
              !Assuming first node of element 1 and last node of element
              !K to be boundaries
              !global node # for interior and exterior at faces, each size 2*K
              use Global_var
              implicit none
              integer::nodeids(Nc,K),k2,f2,vidM,vidP,row,col
              integer::vM(2,K),vP(2,K),i,j
              real(8)::x1,x2,d
                do i=1,Nc
                do j=1,K
                nodeids(i,j)=i+(j-1)*Nc
                enddo
                enddo
                vM=0*vM
                vP=0*vP
                do i=1,K
                vM(1,i)=nodeids(1,i)
                vM(2,i)=nodeids(Nc,i)
                enddo
                do i=1,K
                do j=1,2
                k2=EToE(i,j)
                f2=EToF(i,j)
                !element i 's face j connects element k2 's face f2
                !they are entries of vM
                vidM=vM(j,i)
                vidP=vM(f2,k2)
                !from their global index id, get the coordinates
                col=i
                row=int(j/2)*(Nc-1)+1
                x1=x(row,col)
                col=k2
                row=int(f2/2)*(Nc-1)+1
                x2=x(row,col)
                d=abs(x1-x2)
                if (d<nodetol)then
                        vP(j,i)=vidP
                endif
                enddo
                enddo
                do i=1,K
                do j=1,2
                vmapM(j+(i-1)*2)=vM(j,i)
                vmapP(j+(i-1)*2)=vP(j,i)
                enddo
                enddo
                j=1
                do i=1,2*K
                if(vmapM(i)==vmapP(i))then
                        if(j>2)then
                                exit
                        endif
                        mapB(j)=i
                        j=j+1
                endif
                enddo
                do i=1,2
                vmapB(i)=vmapM(mapB(i))
                enddo
      endsubroutine BuildMaps1D
      subroutine Initial1D
              ! call global variables
              ! Initialization
              use Global_var
              use Global_math
              use Nodal1D
              implicit none
              integer::i,j
              call LegendreGL_(N,r)
              call Vandermonde1D_(N,Nc,r(:,1),V)
              invV=V
              call InverseMatrix(Nc,invV,info)
              call Dmatrix1D_(N,Nc,r(:,1),Dr)
              Lift=Lift1D_(Nc,V)
              call Uniform1D_(cL,cR,K,VX,EToV)
              call NodalCoord_(K,Nc,r(:,1),VX,EToV,x)
              xr=matmul(Dr,x)
              do i=1,Nc
              do j=1,K
              rx(i,j)=1.d0/xr(i,j)
              enddo
              enddo
              nx=0.d0
              Fmask(1)=1
              Fmask(2)=Nc
              do i=1,K
              Fx(1,i)=x(1,i)
              Fx(2,i)=x(Nc,i)
              nx(1,i)=-1.d0
              nx(2,i)=1.d0
              Fscale(1,i)=rx(1,i)
              Fscale(2,i)=rx(Nc,i)
              enddo
              call Connect1D_(K,EToV,EToE,EToF)
              !call BuildMaps1D
      end subroutine Initial1D
      end module Grid1D
