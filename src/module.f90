
!###################################################################################
!#####  Purpose : physical parameter and parameter of the model
!#####  created : 02.7.2015 Johannes Schurer
!#####  changed : xxxx
!###################################################################################
module paramsKonst !Physik
     implicit none
     save 
	real*8		 :: mAOmI   	! mass ratio between m_atom and m_ion
	real*8	 	 :: g 	        ! interatomic interaction strength

	real*8	 	 :: v0 		! barrier height
	real*8	 	 :: omega	! short range cut of energy
	real*8	 	 :: gamma 	! barrier width
 	real*8	 	 :: l 		! trap length

  contains

	!#############################################################################
	!#####  Purpose : Potential 
	!#####  created : 02.7.2015 Johannes
	!#############################################################################
	function V(nx, x)
		implicit none
		integer,intent(in)				:: nx	!number used basis functions
		real*8,dimension(nx),intent(in)			:: x
		real*8,dimension(nx**2)			        :: V
		integer						:: jj,ii
	
                gamma = gamma*4d0*dsqrt(10*omega)
		V = 0d0
		do jj = 1,nx
			do ii = 1,nx
		            V( (jj-1)*nx + ii) = - 1d0/(4d0*(x(jj))**4 + 1d0/omega) + v0*dexp(-gamma*2d0*x(jj)**2)
                            V( (jj-1)*nx + ii) = V( (jj-1)*nx + ii) - 1d0/(4d0*(x(ii))**4 + 1d0/omega) + v0*dexp(-gamma*2d0*x(ii)**2)
                            !if (x(jj) == 0d0) V( (jj-1)*nx + ii) = V( (jj-1)*nx + ii) + g/sqrt(2d0)
			end do
		end do
	end function

	function V1d(nx, x)
		implicit none
		integer,intent(in)				:: nx	!number used basis functions
		real*8,dimension(nx),intent(in)			:: x
		real*8,dimension(nx)			        :: V1d
		integer						:: jj,ii
	
                gamma = gamma*4d0*dsqrt(10*omega)
		V1d = 0d0
		do jj = 1,nx
		   V1d(jj) = - 1d0/(4d0*(x(jj))**4 + 1d0/omega) + v0*dexp(-gamma*2d0*x(jj)**2)
		end do
	end function
	


	function V_harm(nx, x)
		implicit none
		integer,intent(in)				:: nx	!number used basis functions
		real*8,dimension(nx),intent(in)			:: x
		real*8,dimension(nx)			        :: V_harm
		integer						:: jj,ii
	
		V_harm = 0d0
		do jj = 1,nx
		   V_harm(jj) = 0.5d0/l**4 * x(jj)**2
		end do
	end function

	function V_ww(nx, x)
		implicit none
		integer,intent(in)				:: nx	!number used basis functions
		real*8,dimension(nx),intent(in)		:: x
		real*8,dimension(nx,nx)			        :: V_ww
		integer						:: jj,ii
	
                gamma = gamma*4d0*dsqrt(10d0*omega)
		V_ww = 0d0
		do jj = 1,nx
		    do ii = 1,nx
		            V_ww(jj,ii) = - 1d0/((x(jj) - x(ii))**4 + 1d0/omega) + v0*dexp(-gamma*(x(jj)-x(ii))**2)
                    end do
		end do
	end function
end module paramsKonst

!###################################################################################
!#####  Purpose : numerical parameter and const. and Routines for getting parameters of input file
!#####  created : 02.7.2015 Johannes Schurer
!#####  changed : xxxx
!###################################################################################
module grids 
     use paramsKonst
     implicit none
     save 
	complex*16, parameter :: i = (0.0d0,1.0d0)	! complex number
	integer 	      :: j		        ! counter

	integer         		     :: nx   !Number of grid points

	real*8				     :: xMax ! maximal x value
	real*8, dimension(:),allocatable     :: x    ! xgrid
	real*8, dimension(:),allocatable     :: dx   

        integer                              :: m ! number of basis functions to use

	character(len=200)  :: path

 contains

!###################################################################################
!#####  Purpose : Read in from input.txt file
!#####  created : 02.7.2015 Johannes Schurer
!#####  changed : xxxx
!###################################################################################
subroutine get_data
  implicit NONE

!------------| pysical parameters |------------------------------------------------------
  NAMELIST /params/ mAOmI,g,v0,omega,gamma,l

!------------| numerical parameters |------------------------------------------------------
  NAMELIST /main_settings/ path,nx,xMax,m

  OPEN(100,FILE='input.txt')
     READ(UNIT=100,NML=params)
  CLOSE(100)
  OPEN(100,FILE='input.txt')
     READ(UNIT=100,NML=main_settings)
  CLOSE(100)

end subroutine get_data

!###################################################################################
!#####  Purpose : Store data in data folder
!#####  created : 28.7.2015 Johannes Schurer
!#####  changed : xxxx
!###################################################################################
subroutine save_data
  implicit NONE
  integer :: io_error

!------------| pysical parameters |------------------------------------------------------
  NAMELIST /params/ mAOmI,g,v0,omega,gamma,l

!------------| numerical parameters |------------------------------------------------------
  NAMELIST /main_settings/ path,nx,xMax,m

  OPEN(unit=20,file= './'//trim(trim(path)//'params')//'.txt',status='replace',action='write', iostat=io_error)
     write(20,NML=params)
     write(20,NML=main_settings)
  CLOSE(20)

end subroutine save_data

!#############################################################################
!#####  Purpose : Derive 1st derivative
!#####  created : 29.7.2013 Johannes
!#############################################################################
function dfx(nx, dx, fx)
	implicit none
	integer,intent(in)				::nx	!number used basis functions
	real*8,dimension(nx),intent(in)			::dx
	complex*16,dimension(nx),intent(in)		::fx
	complex*16,dimension(nx)			::dfx
	integer						::j
	
	dfx = 0d0
	dfx(1) = 1d0/dx(1)*(fx(2) - fx(1))
	do j = 2,nx-1
		dfx(j) = 1d0/(dx(j))  *(fx(j+1) - fx(j))
	end do
	dfx(nx) = dfx(nx-1)
end function

function dfx2(nx, dx, fx)
	implicit none
	integer,intent(in)				::nx	!number used basis functions
	real*8,dimension(nx),intent(in)			::dx
	complex*16,dimension(nx),intent(in)		::fx
	complex*16,dimension(nx)			::dfx2
	integer						::j
	
	dfx2 = 0d0
	dfx2(1) = 1d0/(2d0*dx(1))*fx(2)
	do j = 2,nx-1
		dfx2(j) = 1d0/(dx(j)+dx(j-1))  *(fx(j+1) - fx(j-1))
	end do
	dfx2(nx) = 1d0/(2d0*dx(nx))  *(fx(nx-1))
end function

!#############################################################################
!#####  Purpose : Derive 2nd derivative
!#####  created : 29.7.2013 Johannes
!#############################################################################
function d2fx(nx, dx, fx)
	implicit none
	integer,intent(in)				::nx	!number used basis functions
	real*8,dimension(nx),intent(in)			::dx
	complex*16,dimension(nx),intent(in)		::fx
	complex*16,dimension(nx)			::d2fx
	integer						::j
		
	d2fx = 0d0
	d2fx(1) = 1d0/dx(1)*(  (fx(3)-fx(2))/dx(2) - (fx(2) - fx(1))/dx(1))
	d2fx(2) = 1d0/dx(2)*(  (fx(4)-fx(3))/dx(3) - (fx(3) - fx(2))/dx(2))
	do j = 3,nx-2
		d2fx(j) = 1d0/(dx(j)+dx(j-1))  *(  (fx(j+2)-fx(j))/(dx(j+1)+dx(j)) - (fx(j) - fx(j-2))/(dx(j-1)+dx(j-2)))
	end do
	d2fx(nx-1) = d2fx(nx-2)
	d2fx(nx) = d2fx(nx-2)
end function

end module grids

