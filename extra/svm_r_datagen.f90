!******************************************************************************
!     svm_r_datagen.f90 - Version 0.1 
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
! of a normalized sinc pulse with a uniform distribution noise. And (b) an equal number of random points to 
! test the training.

program genera_data

implicit none

integer			::			ns, nsn, nd, i, j, io, cac
real			::			tmp,tmp2,pi,alpha1,alpha2
real, dimension(:), allocatable		::		x,noise
character (len=10)		::		tmpc

! get the number of points from the argument line
cac = command_argument_count()
if (cac.NE.1) then
	call get_command_argument(0,tmpc)
	print *,''
	print *,'svm_r_datagen. To use:'
	print *,tmpc,'number_of_points'
	print *,'for example:'
	print *,tmpc,' 1000'
	print *,'will generate 1000 points'
	stop
end if

call get_command_argument(1,tmpc)
tmpc = trim(tmpc)
read(tmpc,'(I5)') ns

nsn = ns
pi = 3.1415
nd = 1
alpha1 = pi/4.0
alpha2 = 4.0

call random

allocate(x(ns),noise(ns))
call random_number(x)
call random_number(noise)
open(unit=3, file='training', form='formatted', action='write')
open(unit=4, file='training.gnuplot',action='write')

write(unit=3,fmt='(I4)') ns
write(unit=3,fmt='(I4)') nd
write(unit=3,fmt='(I4)') nd

do i=1,ns
	tmp2 = x(i)*6.0 - 3.0
	tmp = sin(abs(pi*tmp2))/abs(pi*tmp2) + noise(i)/5.0
	write(unit=3,fmt='(2f16.4)') tmp2, tmp
	write(unit=4,fmt='(2f16.4)') tmp2, tmp
end do

close(unit=3)
close(unit=4)
open(unit=2, file='predict', form='formatted', action='write')

deallocate(x)
allocate(x(nsn))

write(2,fmt=*) nsn
call random_number(x)

do i=1,nsn
	tmp = x(i)*6.0 - 3.0
	write (unit=2,fmt='(f16.4)') tmp
end do

close(2)

deallocate(x)

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
