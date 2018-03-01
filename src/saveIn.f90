!# Routine for saveing functions y(x)
!# 02.7.2015 Johannes Schurer

subroutine saveIN (xdata, ydata ,Dname, Dsize)!, string) ! Speicher Routine
implicit none
  integer , intent(in)  :: Dsize			! array size
  real*8,dimension(Dsize), intent(in) :: xdata		! xDATA
  real*8,dimension(Dsize), intent(in) :: ydata		! yDATA
  !character(len=200)  ,intent(in)  :: string		! Comment line
  integer :: io_error,counter1
  character(len=*), intent(in) :: Dname			! Filename without ending
	open(unit=20,file= './'//trim(Dname)//'.txt',status='replace',action='write', iostat=io_error)
		!write(20,'(a,$)') '%#'
		!write(20,*) trim(string)
    	  if ( io_error == 0) then
		do counter1=1,Dsize
  			write (20,'(1X,ES30.22E3,$)') xdata(counter1) 
			write(20,'(1X,ES30.22E3)')   ydata(counter1)
		enddo	
	  else 
         	 write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', io_error,' aufgetreten'
          end if
   close(unit=20)
    return
end subroutine saveIN

subroutine saveINint (xdata, ydata ,Dname, Dsize)!, string) ! Speicher Routine
implicit none
  integer , intent(in)  :: Dsize			! array size
  integer,dimension(Dsize), intent(in) :: xdata		! xDATA
  integer,dimension(Dsize), intent(in) :: ydata		! yDATA
  !character(len=200)  ,intent(in)  :: string		! Comment line
  integer :: io_error,counter1
  character(len=*), intent(in) :: Dname			! Filename without ending
	open(unit=20,file= './'//trim(Dname)//'.txt',status='replace',action='write', iostat=io_error)
		!write(20,'(a,$)') '%#'
		!write(20,*) trim(string)
    	  if ( io_error == 0) then
		do counter1=1,Dsize
  			write (20,*) xdata(counter1),  ydata(counter1)
		enddo	
	  else 
         	 write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', io_error,' aufgetreten'
          end if
   close(unit=20)
    return
end subroutine saveINint

subroutine saveIN_matrix(matrixdata,xSize, ySize,Dname) ! Speicher Routine
implicit none
  integer , intent(in)  :: xSize,ySize
  real*8,dimension(xSize,ySize), intent(in) :: matrixdata
  integer :: io_error,counter1,counter2
  character(len=*), intent(in) :: Dname
     open(unit=20,file= './'//trim(Dname)//'.txt',status='replace',action='write', iostat=io_error)
    	  if ( io_error == 0) then
		do counter2 = 1,ySize
			write(20,'(*(F14.7))') (  matrixdata(counter1,counter2) , counter1=1,xSize) 
		enddo	
	  else 
         	 write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', io_error,' aufgetreten'
          end if
    close(unit=20)
    return
end subroutine saveIN_matrix


subroutine saveIN_matrixR(matrixdata,xSize, ySize,Dname) ! Speicher Routine
implicit none
  integer , intent(in)  :: xSize,ySize
  complex*16,dimension(xSize,ySize), intent(in) :: matrixdata
  integer :: io_error,counter1,counter2
  character(len=*), intent(in) :: Dname
    open(unit=20,file= './'//trim(Dname)//'.txt',status='replace',action='write', iostat=io_error)
    	  if ( io_error == 0) then
		do counter2=1,ySize
			write(20,'(*(F14.7))')  ( dreal(matrixdata(counter1,counter2)),  counter1=1,xSize) 
		enddo	
	  else 
         	 write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', io_error,' aufgetreten'
          end if
    close(unit=20)
    return
end subroutine saveIN_matrixR

subroutine saveIN_matrixI(matrixdata,xSize, ySize,Dname) ! Speicher Routine
implicit none
  integer , intent(in)  :: xSize,ySize
  complex*16,dimension(xSize,ySize), intent(in) :: matrixdata
  integer :: io_error,counter1,counter2
  character(len=*), intent(in) :: Dname
   open(unit=20,file= './'//trim(Dname)//'.txt',status='replace',action='write', iostat=io_error)
    	  if ( io_error == 0) then
		do counter2=1,ySize
			write(20,'(*(F14.7))')   ( aimag(matrixdata(counter1,counter2)),  counter1=1,xSize) 
		enddo	
	  else 
         	 write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', io_error,' aufgetreten'
          end if
    close(unit=20)
    return
end subroutine saveIN_matrixI
