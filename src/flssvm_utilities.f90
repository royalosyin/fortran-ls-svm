!******************************************************************************
!     flssvm_utilities.f90 - Version 0.1 
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

module flssvm_globals

implicit none

!!!!!!!!!!!!!!!  SVM variables !!!!!!!!!

character (len=50)							::		new_data, out_data, input_data, ga_par
integer										::		kernel_type, ndx, ndy, ns, nsn, io, merror
real(kind=8)								::		s1, s2, gamma, mean, variance, tol_cj
logical										::		regression,train, ga_use
real(kind=8), dimension(:,:), allocatable	::		xin, yin, xnew, ynew, H,alpha
real(kind=8), dimension(:), allocatable		::		b, omega
integer, dimension(:), allocatable			::		position

end module flssvm_globals

module flssvm_utilities

use flssvm_globals

implicit none

contains

subroutine load_input_data

	implicit none
	
	integer			::			i,j

	open(unit = 1, file = input_data, action = 'read', status='old', iostat=io)
	if (io.NE.0) stop 'could not open data file.'

	if (train) then
		read(unit = 1,fmt = *) ns
		read(unit = 1,fmt = *) ndx
		if (regression) then
			read(unit = 1,fmt = *) ndy
		else
			ndy = 1
		end if
		
		call svm_memory_load()
		
		do i=1,ns
				read(unit = 1,fmt = *) xin(i,:),yin(i,:)
		end do

	else
	
		read(unit = 1,fmt = *) ns
		if (regression) then
			read(unit = 1,fmt = *) ndy
		else
			ndy = 1
		end if
		read(unit = 1,fmt = *) ndx
		read(unit = 1,fmt = *) kernel_type
		read(unit = 1,fmt = *) s1
		read(unit = 1,fmt = *) s2
		read(unit = 1,fmt = *) gamma

		call svm_memory_load()
		
		do j = 1,ndy
			do i = 1,ns
				read(unit = 1,fmt='(E16.4)') alpha(i,j)
			end do
		end do
		
		do i = 1,ndy
			read(unit = 1,fmt='(E16.4)') b(i)
		end do
		
		do i = 1,ns
			read(unit = 1,fmt='(E16.4)') xin(i,:)
		end do
		
	end if
	
	close(unit = 1)

end subroutine load_input_data


subroutine svm_memory_load

	if (train) then
		allocate(H(ns,ns),stat=merror)
		if (merror.NE.0) then
			print *,"Memory error. Training matrix could not be allocated"; call svm_end; end if
	end if

	allocate(xin(ns,ndx),stat=merror)
	if (merror.NE.0) then
		print *,"Memory error. data matrix xin could not be allocated"; call svm_end; end if

	allocate(yin(ns,ndy),stat=merror)
	if (merror.NE.0) then
		print *,"Memory error. data matrix yin could not be allocated"; call svm_end; end if

	allocate(b(ndy),stat=merror)
	
	allocate(alpha(ns,ndy),stat=merror)

end subroutine svm_memory_load

subroutine svm_end

	if (allocated(H)) deallocate(H)
	if (allocated(xin)) deallocate(xin)
	if (allocated(yin)) deallocate(yin)
	if (allocated(b)) deallocate(b)
	if (allocated(alpha)) deallocate(alpha)
	
	stop
		
end subroutine svm_end

end module flssvm_utilities
