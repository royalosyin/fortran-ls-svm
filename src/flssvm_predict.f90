!******************************************************************************
!     flssvm_predict.f90 - Version 0.1 
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

module	flssvm_predict

use flssvm_globals
use flssvm_utilities

implicit none

contains

subroutine predict

	real(kind=8)		::		tmp1, tmp2, tmp3
	integer				::		i,j,k,l
	logical, dimension(:), allocatable			::		tmp4

	open(unit = 15, file = new_data, action = 'read', status='old', form='formatted', iostat=io)
	if (io.NE.0) then
		print *,'could not open data file.'; call svm_end;
	end if
	
	read (unit = 15,fmt=*) nsn

	allocate(xnew(nsn,ndx),ynew(nsn,ndy),omega(ns),stat=merror)
	if (merror.NE.0) then
		print *,"Memory error. New data could not be allocated"; call svm_end; end if

	do i=1,nsn
		do j=1,ndx
			read(unit = 15,fmt='(F16.4)', advance = 'no') xnew(i,j)
		end do
		read(unit = 15,fmt='(x)')
	end do

	close(unit = 15)


	! ******************************* REGRESSION
	if (regression) then
		select case (kernel_type)
		case(1)
			do k = 1,ndy
				do i = 1,nsn
					do j = 1,ns
						omega(j) = (dot_product(xnew(i,:),xin(j,:))+s1)**s2
					end do
					ynew(i,k) = dot_product(omega,alpha(:,k))+b(k)
				end do
			end do
		case(2)
			do k = 1,ndy
				do i = 1,nsn
					do j = 1,ns
						tmp1 = dot_product(xnew(i,:),xnew(i,:))
						tmp2 = dot_product(xnew(i,:),xin(j,:))
						tmp3 = dot_product(xin(j,:),xin(j,:))
						omega(j) = dexp(-(tmp1-2.d0*tmp2+tmp3)/(s1*s1))
					end do
					ynew(i,k) = dot_product(omega(1:ns),alpha(:,k))+b(k)
				end do
			end do
		case(3)
			do k = 1,ndy
				do i = 1,nsn
					do j = 1,ns
						tmp1 = 1.d0
						do l = 1,ndx
							tmp2 = (1.d0-s1*s1)/(2.d0*(1.d0-2.d0*s1*dcos(xin(i,l)-xin(j,l))+s1*s1))
							tmp1 = tmp1*tmp2
						end do
						omega(j) = tmp1
					end do
					ynew(i,k) = dot_product(omega,alpha(:,k))+b(k)
				end do
			end do
		end select
	else ! *************************************** CLASSIFICATION
		
		select case (kernel_type)
		case(1)

			do i = 1,nsn
				do j = 1,ns
					omega(j) = (dot_product(xnew(i,:),xin(j,:))+s1)**s2
				end do
				ynew(i,1) = sign(1.d0,dot_product(omega,alpha(:,1))+b(1))
			end do

		case(2)

			do i = 1,nsn
				do j = 1,ns
					tmp1 = dot_product(xnew(i,:),xnew(i,:))
					tmp2 = dot_product(xnew(i,:),xin(j,:))
					tmp3 = dot_product(xin(j,:),xin(j,:))
					omega(j) = dexp(-(tmp1-2.d0*tmp2+tmp3)/(s1*s1))
				end do
				ynew(i,1) = sign(1.d0,dot_product(omega,alpha(:,1))+b(1))
			end do

		case(3)

			do i = 1,nsn
				do j = 1,ns
					tmp1 = 1.d0
					do k = 1,ndx
						tmp2 = (1.d0-s1*s1)/(2.d0*(1.d0-2.d0*s1*dcos(xin(i,k)-xin(j,k))+s1*s1))
						tmp1 = tmp1*tmp2
					end do
					omega(j) = tmp1
				end do
				ynew(i,1) = sign(1.d0,dot_product(omega,alpha(:,1))+b(1))
			end do

		end select
	end if

	if (ga_use.EQV..FALSE.) then
		open(unit = 3, file = out_data, form = 'formatted', action = 'write', status = 'replace', iostat=io)
		
		do i = 1,nsn
			do j = 1,ndx
				write (unit = 3, fmt = '(E16.4,x)', advance = 'no') xnew(i,j)
			end do
			do j = 1,ndy
				write (unit = 3, fmt = '(E16.4,x)', advance = 'no') ynew(i,j)
			end do
			write (unit = 3, fmt = '(x)') 
		end do
		
		close(unit = 3)
		
		deallocate(xnew,ynew,omega)
		
		print *,''
		print *,'Prediction completed.'
		print *,''
		
		call svm_end
	else
	
	
	end if
	
end subroutine predict

end module flssvm_predict
