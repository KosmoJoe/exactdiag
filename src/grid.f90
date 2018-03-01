!###################################################################################
!#####  Purpose : Creates equidistant grid
!#####  created : 26.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################


  subroutine create_eqgrid(n_x, x_min, x_max, x, dx)
	implicit none
	integer,intent(in)::n_x		!number of intervals/grid points
	real*8,intent(in)::x_min	!Lower border for grid
	real*8,intent(in)::x_max	!Upper border for grid
	real*8,dimension(n_x),intent(out)::x		!@out : grid points
	real*8,dimension(n_x),intent(out)::dx		!@out : diffenrential
	integer:: j

	dx=(x_max-x_min)/dble(n_x)
	do j=1, n_x
		x(j)=x_min+(j-0.5d0)*dx(j)
	end do
	return 
  end subroutine create_eqgrid

  subroutine create_eqgrid2(n_x, x_min, x_max, x, dx)
	implicit none
	integer,intent(in)::n_x		!number of intervals/grid points
	real*8,intent(in)::x_min	!Lower border for grid
	real*8,intent(in)::x_max	!Upper border for grid
	real*8,dimension(n_x),intent(out)::x		!@out : grid points
	real*8,dimension(n_x),intent(out)::dx		!@out : diffenrential
	integer:: j

	dx=(x_max-x_min)/dble(n_x-1)
	do j=1, n_x
		x(j)=x_min+(j-1)*dx(j)
	end do
	return 
  end subroutine create_eqgrid2

!###################################################################################
!#####  Purpose : Creates piecewise equidistant grid
!#####  created : 26.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################

subroutine create_pwgrid (n_x,brd,n_int,n_subint,x,dx) 
     implicit NONE
	integer, intent(in)	::n_x,n_int  			!number of intervals /grid points
	real*8 , dimension(n_int+1),intent(in)	 :: brd		!borders for subgrids
	integer, dimension(n_int),intent(in)	 :: n_subint	!number of grid points for each subgrid 
	real*8 , dimension(n_x),intent(out)	 :: x		!@out : grid points
	real*8 , dimension(n_x),intent(out) 	 :: dx		!@out : diffenrential

	integer		:: j ,jj, js
	real*8 dx_hlp
	
	if (n_x /= sum(n_subint)) then
	    write(*,*) "n_tot not equals sum of single subintervals" 
	    stop 
	end if
	js = 0
        do j= 1,n_int
	    dx_hlp = (brd(j+1) - brd(j))/n_subint(j)
	    do jj = 1,n_subint(j)
		dx(js+jj) 	= dx_hlp
		x(js+jj)	= brd(j) + (jj-0.5d0)*dx_hlp
	    end do
	    js = sum(n_subint(1:j))
	end do
	return
end subroutine create_pwgrid

!###################################################################################
!#####  Purpose : Creates general piecewise grid
!#####  created : 26.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################

subroutine create_gpwgrid (n_x,brd,n_int,n_subint,typ,x,dx) 
     implicit NONE
	integer, intent(in)	::n_x,n_int  			!number of intervals /grid points
	real*8 , dimension(n_int+1),intent(in)	 :: brd		!borders for subgrids
	integer, dimension(n_int),intent(in)	 :: n_subint	!number of grid points for each subgrid 
	integer, dimension(n_int),intent(in)	 :: typ		!Type for subgrid
	real*8 , dimension(n_x),intent(out)	 :: x		!@out : grid points
	real*8 , dimension(n_x),intent(out) 	 :: dx		!@out : diffenrential

	integer		:: j , temp
	real*8 dx_hlp
	
	if (n_x /= sum(n_subint)) then
	    write(*,*) "n_tot not equals sum of single subintervals" 
	    stop 
	end if
	temp = 0
	do j = 1,n_int
		 select case (typ(j))
			case (1) !Equidistant
			    call create_eqgrid(n_subint(j), brd(j), brd(j+1), x(temp+1:temp+n_subint(j)), dx(temp+1:temp+n_subint(j)))
			case (2) !Gauss
			    call create_gcgrid(brd(j), brd(j+1), x(temp+1:temp+n_subint(j)), dx(temp+1:temp+n_subint(j)),n_subint(j))
			case default 
			    write (*,*) 'This should not occur!!!'
		end select 
		temp = sum(n_subint(1:j))
	end do

	return
end subroutine create_gpwgrid

!###################################################################################
!#####  Purpose : Creates gauss-chebyshev grid
!#####  created : 26.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################

  subroutine create_gcgrid(x_min,x_max,x,dx,n_x)
     !INTRINSIC dble
	implicit none
	integer,intent(in)::n_x		!number of intervals/grid points
	real*8,intent(in)::x_min	!Lower border for grid
	real*8,intent(in)::x_max	!Upper border for grid
	real*8,dimension(n_x),intent(out)::x		!@out : grid points
	real*8,dimension(n_x),intent(out)::dx		!@out : diffenrential
	integer:: j

     real*8 tpi_hlp,apu_hlp,ai_hlp
     tpi_hlp =8.0d0*atan(1.0d0)
     apu_hlp =(x_max-x_min)/dble(n_x+1)
     do j=1,n_x
        ai_hlp =dble(j)/dble(n_x + 1)
        x(j)=(ai_hlp - dsin(tpi_hlp*ai_hlp)/tpi_hlp)*(x_max - x_min) + x_min
        dx(j)=(1.0d0-dcos(tpi_hlp*ai_hlp))*apu_hlp
     enddo
     return
  end subroutine create_gcgrid



!###################################################################################
!#####  Purpose : Interpolates a real valued function at position x
!#####  created : 26.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################

   subroutine lin_interpol(xi,yi,n_xi,x,y)
      implicit none
      integer, intent(in) :: n_xi  			!length of xi and yi
      real*8, dimension(n_xi), intent(in) :: xi 	!array of x-values  
      real*8, dimension(n_xi), intent(in) :: yi 	!the corresponding array of y-values
      real*8, intent(in) :: x  				!Position to interpolate
      real*8, intent(out) :: y				!@out : interpolated y-value at position x (in)

      integer :: im,ip,mid
      real*8 :: a

      ! find neighboring grid points of x using the bisection method
      im=1; ip=n_xi
      do
         if (x<=xi(im+1)) then
            ip=im+1
            exit
         end if
         if (x>=xi(ip-1)) then
            im=ip-1
            exit
         end if
         mid=(im+ip)/2
         if (x>=xi(mid)) im=mid
         if (x<xi(mid)) ip=mid
      end do

      ! interpolate
      if (xi(ip)-xi(im)==0d0) then
         write(*,*) 'Two identical x-values in interpolation routine. Program stopped!'
         stop
      end if
      a = (xi(ip)-x)/(xi(ip)-xi(im))
      y = a*yi(im) + (1d0-a)*yi(ip)
      return
   end subroutine lin_interpol


!###################################################################################
!#####  Purpose : Interpolates a real valued function on a new grid
!#####  created : 26.7.2013 Johannes
!#####  changed : xxxx
!###################################################################################

   subroutine lin_interpolate(xi,yi,n_xi,x,y,n_x)
      implicit none
      integer, intent(in) :: n_xi,n_x  			!n_xi: length of xi and yi ; 	n_x:length of x and y
      real*8, dimension(n_xi), intent(in) :: xi 	!array of x-values  
      real*8, dimension(n_xi), intent(in) :: yi 	!the corresponding array of y-values
      real*8, dimension(n_x), intent(in) :: x  		!Grid to interpolate
      real*8, dimension(n_x), intent(out) :: y		!@out : interpolated y-value at positions x (in)
	  integer:: j
      y = 0d0
      do j = 1,n_xi
	call lin_interpol(xi,yi,n_xi,x(j),y(j))
      end do
      return
   end subroutine lin_interpolate

