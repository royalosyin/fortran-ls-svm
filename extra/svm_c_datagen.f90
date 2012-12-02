!******************************************************************************
!     svm_c_datagen.f90 - Version 0.1 
!	  Copyright (C) Jose Colmenares
! 	  Developed by:	Jose Colmenares (jbcolmenares@gmail.com)
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.     
!*******************************************************************************

! this program generates (a) a number of points (the number given by console imput)
! of which those inside a circle of radio 0.2 are labeled with a one, and
! those outside are labeled with cero. (b) an equal number of random points to 
! test the training.

program genera_data

implicit none

integer			::			ns, nsn, nd, i, j, io, cac
real			::			tmp,tmp2
real, dimension(:), allocatable		::		x,y
character (len=10)		::		tmpc

cac = command_argument_count()
if (cac.NE.1) then
	call get_command_argument(0,tmpc)
	print *,''
	print *,'svm_c_datagen. To use:'
	print *,tmpc,'number_of_points'
	print *,'for example:'
	print *,tmpc,' 1000'
	print *,'will generate 1000 points'
	stop
end if

call get_command_argument(1,tmpc)
tmpc = trim(tmpc)
read(tmpc,'(I5)') ns

nsn = 1*ns
nd = 2

call random

allocate(x(ns),y(ns))
call random_number(x)
call random_number(y)

open(unit=3, file='training', form='formatted', action='write')
open(unit=4, file='training.gnuplot',action='write')

write(unit=3,fmt='(I4)') ns
write(unit=3,fmt='(I4)') nd

do i=1,ns
	x(i) = x(i)*2.0 - 1.0
	y(i) = y(i)*2.0 - 1.0
end do

do i=1,ns
	if (sqrt(x(i)**2+y(i)**2).LT.0.2) then
		tmp = 1.0
	else
		tmp = -1.0
	end if
	write(unit=3,fmt='(3f16.4)') x(i),y(i),tmp
	write(unit=4,fmt='(3f16.4)') x(i),y(i),tmp
end do

close(unit=3)
close(unit=4)
open(unit=2, file='predict', form='formatted', action='write')

deallocate(x,y)
allocate(x(nsn),y(nsn))

write(2,fmt=*) nsn
call random_number(x)
call random_number(y)

do i=1,ns
	x(i) = x(i)*2.0 - 1.0
	y(i) = y(i)*2.0 - 1.0
end do

do i=1,nsn
	write (unit=2,fmt='(2f16.4)') x(i),y(i)
end do

close(2)

deallocate(x,y)

contains

subroutine random

  implicit none

  character(len=10)                   ::  date, time, zone
  integer, dimension(:), allocatable  ::  seed
  integer                             ::    n
  integer                             ::    i,io
  integer, dimension(8)               ::  values
  
  n = 10
  
  allocate(seed(n),stat = io)
    if (io.NE.0) then
        print *,'Memory error while runing genetic algorithms'
        stop
    end if

  do i = 1,n
    call date_and_time(date,time,zone,values)
    seed(i) = real(values(8))
  end do
  
  call random_seed(size=n)

  call random_seed(put=seed)

  deallocate(seed)
  
end subroutine random


end program genera_data
