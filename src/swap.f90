!###################################################################################
!#####  Purpose : Sort of Re of eigenvalues with same order change in eigenfunctions
!#####  created : 30.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################
subroutine SwapReal(aX,pX,ni)!aX=Array(Länge ni) mit Eigenwerten, pX=Matrix(ni x ni, erster Index: k, zweiter Index: Eigenwert), ni=Arraylänge
     implicit NONE
     integer ni
     complex*16 pX(ni,ni),aX(ni),swap
     integer i,j,kk
     do i=1,ni-1
        do j=i+1,ni
           if (dreal(aX(i)).GE.dreal(aX(j))) then
              swap=aX(i)
              aX(i)=aX(j)
              aX(j)=swap
              do kk=1,ni
                 swap=pX(kk,i)
                 pX(kk,i)=pX(kk,j)
                 pX(kk,j)=swap
              enddo
           endif
        enddo
     enddo
    return
end subroutine SwapReal
!###################################################################################
!#####  Purpose : Sort of Im of eigenvalues with same order change in eigenfunctions
!#####  created : 30.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################
subroutine SwapImag(aX,pX,ni)!aX=Array(Länge ni) mit Eigenwerten, pX=Matrix(ni x ni, erster Index: k, zweiter Index: Eigenwert), ni=Arraylänge
     implicit NONE
     integer ni
     complex*16 pX(ni,ni),aX(ni),swap
     integer i,j,kk
     do i=1,ni-1
        do j=i+1,ni
           if (aimag(aX(i)).GE.aimag(aX(j))) then
              swap=aX(i)
              aX(i)=aX(j)
              aX(j)=swap
              do kk=1,ni
                 swap=pX(kk,i)
                 pX(kk,i)=pX(kk,j)
                 pX(kk,j)=swap
              enddo
           endif
        enddo
     enddo
    return
end subroutine SwapImag


!###################################################################################
!#####  Purpose : Sort of RE of eigenvalues with same order change in left AND rigth eigenfunctions
!#####  created : 30.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################
subroutine SwapRealLR(aX,pX,gX,ni)!aX=Array(Länge ni) mit Eigenwerten, pX=Matrix(ni x ni, erster Index: k, zweiter Index: Eigenwert), ni=Arraylänge
     implicit NONE
     integer ni
     complex*16 pX(ni,ni),gX(ni,ni),aX(ni),swap
     integer i,j,kk
     do i=1,ni-1
        do j=i+1,ni
           if (dreal(aX(i)).GE.dreal(aX(j))) then
              swap=aX(i)
              aX(i)=aX(j)
              aX(j)=swap
              do kk=1,ni
                 swap=pX(kk,i)
                 pX(kk,i)=pX(kk,j)
                 pX(kk,j)=swap
                 swap=gX(kk,i)
                 gX(kk,i)=gX(kk,j)
                 gX(kk,j)=swap
              enddo
           endif
        enddo
     enddo
    return
end subroutine SwapRealLR

!###################################################################################
!#####  Purpose : Sort of Im of eigenvalues with same order change in left AND rigth eigenfunctions
!#####  created : 30.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################
subroutine SwapImagLR(aX,pX,gX,ni)!aX=Array(Länge ni) mit Eigenwerten, pX=Matrix(ni x ni, erster Index: k, zweiter Index: Eigenwert), ni=Arraylänge
     implicit NONE
     integer ni
     complex*16 pX(ni,ni),gX(ni,ni),aX(ni),swap
     integer i,j,kk
     do i=1,ni-1
        do j=i+1,ni
           if (aimag(aX(i)).GE.aimag(aX(j))) then
              swap=aX(i)
              aX(i)=aX(j)
              aX(j)=swap
              do kk=1,ni
                 swap=pX(kk,i)
                 pX(kk,i)=pX(kk,j)
                 pX(kk,j)=swap
                 swap=gX(kk,i)
                 gX(kk,i)=gX(kk,j)
                 gX(kk,j)=swap
              enddo
           endif
        enddo
     enddo
    return
end subroutine SwapImagLR





!######Real Valued
subroutine Swap2DR(aX,pX,ni)!aX=Array(Länge ni) mit Eigenwerten, pX=Matrix(ni x ni, erster Index: k, zweiter Index: Eigenwert), ni=Arraylänge
     implicit NONE
     integer ni
     real*8 pX(ni,ni),aX(ni),swap
     integer i,j,kk
     do i=1,ni-1
        do j=i+1,ni
           if (aX(i).GE.aX(j)) then
              swap=aX(i)
              aX(i)=aX(j)
              aX(j)=swap
              do kk=1,ni
                 swap=pX(kk,i)
                 pX(kk,i)=pX(kk,j)
                 pX(kk,j)=swap
              enddo
           endif
        enddo
     enddo
    return
end subroutine Swap2DR
!######Real Valued + 2nd array
subroutine Swap2DR2(aX,pX,ni,bX)!aX=Array(Länge ni) mit Eigenwerten, pX=Matrix(ni x ni, erster Index: k, zweiter Index: Eigenwert), ni=Arraylänge , bX = array(ni) 
     implicit NONE
     integer ni
     real*8 pX(ni,ni),aX(ni),bX(ni),swap
     integer i,j,kk
     do i=1,ni-1
        do j=i+1,ni
           if (aX(i).GE.aX(j)) then
              swap=aX(i)
              aX(i)=aX(j)
              aX(j)=swap
              swap=bX(i)
              bX(i)=bX(j)
              bX(j)=swap
              do kk=1,ni
                 swap=pX(kk,i)
                 pX(kk,i)=pX(kk,j)
                 pX(kk,j)=swap
              enddo
           endif
        enddo
     enddo
    return
end subroutine Swap2DR2
