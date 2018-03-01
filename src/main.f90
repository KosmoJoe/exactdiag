
!############   Main for the calculation of Three body problem curves #################
!#####  Purpose : Solve Relative Part of three body problem
!#####  created : 02.7.2015 Johannes
!#############################################################################
program main
  use grids        
  implicit none
  complex*16,dimension(:,:),allocatable         :: Hphi,H1B		!Matrix for diagonalization
  real*8,dimension(:),allocatable	        :: Vpot
  real*8,dimension(:,:),allocatable             :: norm
  integer					:: ii,jj
  character(50)                                 :: xstring

  complex*16,dimension(:),allocatable	:: eigen,epsi !eigenvalues 
  complex*16,dimension(:,:),allocatable :: vL,vR
  complex*16,dimension(:),allocatable	:: work
  real*8,dimension(:),allocatable	:: rwork
  integer				:: info




  
!####################  Load from input.txt 
  call get_data()  
  call system('mkdir -p '//trim(path))
  call save_data()

!#################### Allocate
  write (*,'(a,$)') 'Allocate'
  allocate(x(nx),dx(nx),H1B(nx,nx),epsi(nx),Vpot(nx))
  allocate(Hphi(m**2,m**2),eigen(m**2))
  allocate(vL(nx,nx),vR(nx,nx),work(4*nx),rwork(2*nx),norm(nx,nx))
  write (*,'(16X,a)') 'Done'
  
!#################### Create Grid
  write (*,'(a,$)') 'Create Grid'
  call create_eqgrid2(nx, -xmax, xmax, x, dx(1))
  write (*,'(13X,a)') 'Done'

!#################### Create Matrix 1Body
   write (*,'(a,$)') 'Prepare One Body'
   Vpot = V1d(nx,x)
   call create_matrix1B(nx,H1B,Vpot,dx)
 !#################### Diagonalize Matrix
   call zgeev('V','V',nx,H1B,nx,epsi,vL,nx,vR,nx,work,4*nx,rwork,info )
   !write(*,*) work(1),info
!#################### Sort
   call SwapRealLR(epsi,vL,vR,nx)	
!#################### Normalize
   do ii=1,nx
	vR(:,ii)=vR(:,ii)/dsqrt(dx(:))
        !vL(:,ii)=vL(:,ii)/dsqrt(dx(:))  
   end do
   do ii=1,nx
	do jj=1,nx
	    !norm(ii,jj) = sum(conjg(vL(:,ii))*vR(:,jj)*dx(:))
            norm(ii,jj) = sum(conjg(vR(:,ii))*vR(:,jj)*dx(:))
	end do	
	vR(:,ii)=vR(:,ii)/dsqrt(abs(norm(ii,ii)))
	!vL(:,ii)=vL(:,ii)/dsqrt(abs(norm(ii,ii)))
   end do
   call saveIN_matrixR(norm,nx, nx,trim(path)//'norm')

   call saveIN(x,dx, trim(path)//'grid',nx)
   call saveIN(x,Vpot,trim(path)//'potN',nx)
   call saveIN(dreal(epsi),aimag(epsi), trim(path)//'epsi',nx)

   call saveIN_matrixR(vR,nx, nx,trim(path)//'eigenV1B')
   !do ii = 1,5
   !  write( xstring, '(i10)' )  ii
   !  call saveIN(dreal(vR(:,ii)),aimag(vR(:,ii)), trim(path)//'eigenV1B_'//trim(adjustl(xstring)),nx)
   !end do
   write (*,'(8X,a)') 'Done'
  

!#################### Create Matrix 2Body
   write (*,'(a,$)') 'Prepare Two Body'
   call create_matrixPhi(epsi,vR,Hphi)
   deallocate(vR,vL,work,rwork,norm)
   allocate(vL(m**2,m**2),vR(m**2,m**2),work(4*m**2),rwork(2*m**2))
   call saveIN_matrixR(Hphi,m**2, m**2,trim(path)//'Hr')
   write (*,'(8X,a)') 'Done'
 !#################### Diagonalize Matrix
   write (*,'(a,$)') 'Main'
   call zgeev('V','V',m**2,Hphi,m**2,eigen,vL,m**2,vR,m**2,work,4*m**2,rwork,info )
   !write(*,*) work(1),info
   write (*,'(20X,a)') 'Done'
!#################### Sort
   write (*,'(a,$)') 'Sort'
   call SwapRealLR(eigen,vL,vR,m**2)
   write (*,'(20X,a)') 'Done'


!#################### Save
  write (*,'(a,$)') 'Save'

  !call saveIN_matrixR(vR,m**2, m**2,trim(path)//'eigenV')
  do ii = 1,5
     write( xstring, '(i10)' )  ii
     !call saveIN(dreal(vR(:,ii)),aimag(vR(:,ii)), trim(path)//'eigenV_'//trim(adjustl(xstring)),m**2)
     call saveIN_matrixR(reshape(vR(:,ii),(/m, m/)),m,m, trim(path)//'eigenV_'//trim(adjustl(xstring)))
     call saveIN_matrixI(reshape(vR(:,ii),(/m, m/)),m,m, trim(path)//'eigenVI_'//trim(adjustl(xstring)))
  end do
  call saveIN(dreal(eigen),aimag(eigen), trim(path)//'eigen',m**2)
  write (*,'(20X,a)') 'Done'  
  
!#################### Clean
  write (*,'(a,$)') 'Clean'
  deallocate(x,dx,Hphi,eigen,Vpot,H1B,epsi)
  deallocate(vR,vL,work,rwork)
  write (*,'(19X,a)') 'Done'  


end program main

!#############################################################################
!#####  Purpose : Derive Hamilton for two reltiave Particle on x
!#####  created : 07.7.2015 Johannes
!#############################################################################
subroutine create_matrix(nx,H,Vpot,dx)
	use paramsKonst
	implicit none
	integer,intent(in)					:: nx	! number of grid points
	real*8 ,intent(in)					:: dx	
	real*8,dimension(nx**2),intent(in)  			:: Vpot	! vector of potential
	complex*16,dimension(nx**2,nx**2),intent(inout)  	:: H	! Matrix for diagonalization
	real*8,dimension(nx,nx)                	:: help	
	integer:: jj,ii

	H = 0d0
!#################### kinetic I
	do jj=1,nx**2
		do ii=1,nx**2
                      if (jj == ii) then
		      	H(jj,ii) = -2d0
                      elseif ( (abs(jj-ii)==1) .AND.  (.NOT.(mod(ii, nx) == 0 .AND. jj/nx > ii/nx)) .AND. (.NOT.(mod(jj, nx) == 0 .AND. ii/nx > jj/nx)) ) then
			H(jj,ii) = 1d0
                      endif
		end do		
	end do
        H = H*(-1d0/(2d0*dx**2))
!#################### kinetic II
	do jj=1,nx**2
		do ii=1,nx**2
		      if (jj == ii) then
		      	H(jj,ii) = H(jj,ii) + 1d0/dx**2  !-2*(-0.5/dx^2)
                      elseif (abs(jj-ii)==nx) then
			H(jj,ii) = H(jj,ii) - 1d0/(2d0*dx**2) !1*(-0.5/dx^2)
                      endif
		end do		
	end do
!#################### kinetic III
        help =0d0
	do jj=1,nx
		do ii=1,nx
		      if (jj-ii==1) then
		      	help(jj,ii) = 1d0
                      elseif (jj-ii== -1) then
			help(jj,ii) = -1d0
                      endif
		end do		
	end do
        help = help*(-1d0/(8d0*dx**2))

	do jj=1,nx
		do ii=1,nx
		      if (jj-ii==1) then
		      	H(nx*(jj-1)+1:nx*jj,nx*(ii-1)+1:nx*ii) =  H(nx*(jj-1)+1:nx*jj,nx*(ii-1)+1:nx*ii) +  help
                      elseif (jj-ii== -1) then
			H(nx*(jj-1)+1:nx*jj,nx*(ii-1)+1:nx*ii) =  H(nx*(jj-1)+1:nx*jj,nx*(ii-1)+1:nx*ii) + transpose(help)
                      endif
		end do		
	end do
!#################### potential
	do jj=1,nx**2
		H(jj,jj) = H(jj,jj) + Vpot(jj)
	end do
	
	return
end subroutine create_matrix

!#############################################################################
!#####  Purpose : Derive Hamilton for one reltiave Particle on x
!#####  created : 07.7.2015 Johannes
!#############################################################################
subroutine create_matrix1B(nx,H,Vpot,dx)
	use paramsKonst
	implicit none
	integer,intent(in)					:: nx	! number of grid points
	real*8 ,intent(in)					:: dx	
	real*8,dimension(nx),intent(in)  			:: Vpot	! vector of potential
	complex*16,dimension(nx,nx),intent(inout)  	:: H	! Matrix for diagonalization
	integer:: jj,ii

	H = 0d0
!#################### kinetic I
	do jj=1,nx
		do ii=1,nx
                      if (jj == ii) then
		      	H(jj,ii) = -2d0
                      elseif (abs(jj-ii)==1) then
			H(jj,ii) = 1d0
                      endif
		end do		
	end do
        H = H*(-1d0/(2d0*dx**2))
!#################### potential
	do jj=1,nx
		H(jj,jj) = H(jj,jj) + Vpot(jj)
	end do
	
	return
end subroutine create_matrix1B
	

!#############################################################################
!#####  Purpose : Derive Hamilton for two reltiave Particle on phi_1b
!#####  created : 07.7.2015 Johannes
!#############################################################################
subroutine create_matrixPhi(epsi,phi,H)
	use grids
	implicit none
	complex*16,dimension(nx),intent(in)  	        :: epsi	! single particle eigenenergies
	complex*16,dimension(nx,nx),intent(in)          :: phi	! single particle eigenfunctions
	complex*16,dimension(m**2,m**2),intent(inout)  	:: H	! Matrix for diagonalization

	complex*16,dimension(:,:),allocatable          :: dM	!
	complex*16,dimension(:,:,:,:),allocatable      :: dD	!
	complex*16,dimension(:,:),allocatable          :: dphi	! derivative of single particle eigenfunctions
        complex*16                          :: temp
        integer,dimension(:,:),allocatable  :: tab
	integer:: jj,ii,k,q,kp,qp,pairI,iip,jjp

        allocate(dM(m,m),dphi(nx,nx),tab(m**2,2),dD(m,m,m,m))
	H = 0d0
!#################### Derive dM(k,kp)= \int \phi_k^* \dphi_kp \dx
	do jj=1,nx
            dphi(:,jj) = dfx2(nx, dx, phi(:,jj))
	end do

	!write(*,*) dphi(:,1)
        !write(*,*) sum(conjg(phi(:,1))*dphi(:,1)*dx(:))
        !write(*,*) sum(phi(:,2)*phi(:,1)*dx(:))
        !write(*,*) sum(dphi(:,1)*dx(:))

        call saveIN_matrixR(dphi,nx, nx,trim(path)//'dphi')
        dM = 0d0
	do jj=1,m
	      do ii=jj,m
                temp = sum(conjg(phi(:,jj))*dphi(:,ii)*dx(:))
		dM(jj,ii) = temp
		if (.NOT.(jj==ii)) then
                	dM(ii,jj) = -temp
                end if
              end do
	end do
        call saveIN_matrixR(dM,m, m,trim(path)//'dM')

!#################### Derive D(k,q,qp,kq)= \int \phi_k^* \phi_q^* \phi_kp \phi_qp \dx

        dD = 0d0
	do jj=1,m
	      do ii=1,m
		  do iip=1,m
		      do jjp=1,m
                	  dD(jj,ii,iip,jjp) = sum(conjg(phi(:,jj))*conjg(phi(:,ii))*phi(:,iip)*phi(:,jjp)*dx(:))
		      end do
		  end do
              end do
	end do
	dD = dD*g/dsqrt(2d0)
!#################### pairI Table
        tab = 0
	do jj=1,m**2
                k = pairI(jj,1,m) 
                q = pairI(jj,2,m)
                tab(jj,1) = k
                tab(jj,2) = q 	
	end do
        call saveINint(tab(:,1),tab(:,2),trim(path)//'tab',m**2)
!#################### kinetic I	
	do jj=1,m**2
                k = tab(jj,1) 
                q = tab(jj,2) 
		do ii=1,m**2
		      kp = tab(ii,1) 
                      qp = tab(ii,2) 
                      if (jj == ii) then
		      	H(jj,jj) = H(jj,ii) + epsi(k)+epsi(q)
                      endif
                      H(jj,ii) = H(jj,ii) -1d0/2d0*dM(k,kp)*dM(q,qp) + dD(k,q,qp,kp)
		end do		
	end do
	
	deallocate(dM,dphi,tab,dD)

	return
end subroutine create_matrixPhi



!#############################################################################
!#####  Purpose : Derives the index of a pair (k,q)
!#####  created : 07.7.2015 Johannes
!#############################################################################
    function pair(k,q,M) result(j)
        implicit none 
        integer, intent(in)  :: k,q,M
        integer              :: j
	j = (k-1)*M + q
    end function pair

!#############################################################################
!#####  Purpose : Derives the pair (k,q) from the index j
!#####  created : 07.7.2015 Johannes
!#############################################################################
    function pairI(j,n,M) result(k)
        implicit none 
        integer, intent(in)  :: j,n,M
        integer              :: tmp,k,ind
	
        tmp = j
        ind = 1
	do while (tmp>0)
          tmp = tmp - M
          ind = ind + 1
	end do
        k = ind-1   ! Find k

        if(n==2) k = j - (k-1)*M  !Derive q from k and return (as k)

    end function pairI

