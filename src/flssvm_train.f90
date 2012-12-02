!******************************************************************************
!     flssvm_train.f90 - Version 0.1 
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

module		flssvm_train

use flssvm_globals
use flssvm_utilities

implicit none

contains

subroutine training

	real (kind=8)								::		rtmp1, rtmp2, rtmp3, s
	real (kind=8), dimension(:), allocatable	::		n, v, ones
	integer										::		i,j,k

	allocate(n(ns),v(ns),ones(ns),stat=merror)
	if (merror.NE.0) then
		print *,"Memory error. Training could be carryed out"; call svm_end; end if	
	
	if (regression) then
	! **************************************** REGRESSION
		select case (kernel_type)
			case(1)
				do j = 1,ns
					do i = 1,j
						H(i,j) =  (dot_product(xin(i,:),xin(j,:)) + s1)**s2
						H(j,i) = H(i,j) 
					end do
				end do
			case(2)
				do j = 1,ns
					do i = 1,j
						rtmp1 = dot_product(xin(i,:),xin(i,:))
						rtmp2 = dot_product(xin(i,:),xin(j,:))
						rtmp3 = dot_product(xin(j,:),xin(j,:))
						H(i,j) = dexp(-(rtmp1-2.d0*rtmp2+rtmp3)/(s1*s1))
						H(j,i) = H(i,j)
					end do
				end do
			case(3)
				do j = 1,ns
					do i = i,j
						rtmp1 = 1.d0
						do k = 1,ndx
							rtmp2 = (1.d0-s1*s1)/(2.d0*(1.d0-2.d0*s1*dcos(xin(i,k)-xin(j,k))+s1*s1))
							rtmp1 = rtmp1*rtmp2
						end do
						H(i,j) = rtmp1
						H(j,i) = rtmp1
					end do
				end do
		end select
	else ! ************************************ CLASSIFICATION
		select case (kernel_type)
			case(1)
				do j = 1,ns
					do i = 1,j
						H(i,j) =  yin(j,1)*yin(i,1)*(dot_product(xin(i,:),xin(j,:)) + s1)**s2
						H(j,i) = H(i,j) 
					end do
				end do
			case(2)
				do j = 1,ns
					do i = 1,j
						rtmp1 = dot_product(xin(i,:),xin(i,:))
						rtmp2 = dot_product(xin(i,:),xin(j,:))
						rtmp3 = dot_product(xin(j,:),xin(j,:))
						H(i,j) = yin(j,1)*yin(i,1)*dexp(-(rtmp1-2.d0*rtmp2+rtmp3)/(s1*s1))
						H(j,i) = H(i,j)
					end do
				end do
			case(3)
				do j = 1,ns
					do i = i,j
						rtmp1 = 1.d0
						do k = 1,ndx
							rtmp2 = (1.d0-s1*s1)/(2.d0*(1.d0-2.d0*s1*dcos(xin(i,k)-xin(j,k))+s1*s1))
							rtmp1 = rtmp1*rtmp2
						end do
						H(i,j) = yin(j,1)*yin(i,1)*rtmp1
						H(j,i) = H(i,j)
					end do
				end do
		end select
	end if

	do i = 1,ns
		H(i,i) = 1/gamma + H(i,i)
	end  do

    ones = 1.0

	if (regression) then	
		do j = 1,ndy
			n = conjgrad(H,ones)
			v = conjgrad(H,yin(1:ns,j))
			s = dot_product(ones,n)
			b(j) = dot_product(n,yin(1:ns,j))/s

			do i = 1,ns
				alpha(i,j) = v(i) - n(i)*b(j)
			end do
		end do
	else
		n = conjgrad(H,yin(:,1))
		v = conjgrad(H,ones)
		s = dot_product(yin(:,1),n)
		b(1) = dot_product(n,ones)/s
		
		do i = 1,ns
			alpha(i,1) = v(i) - n(i)*b(1)
		end do
	end if

	if (ga_use.EQV..FALSE.) then
		open(unit = 10, file = out_data, status = 'replace')
		write(10,fmt=*) ns
		if (regression) write(10,fmt=*) ndy
		write(10,fmt=*) ndx
		write(10,fmt=*) kernel_type
		write(10,fmt=*) s1
		write(10,fmt=*) s2
		write(10,fmt=*) gamma
		
		do j = 1,ndy
			do i = 1,ns
				write(10,fmt='(E16.4)') alpha(i,j)
			end do
		end do
		
		do i = 1,ndy
			write(10,fmt='(E16.4)') b(i)
		end do
		
		do i = 1,ns
			write(10,fmt='(E16.4)') xin(i,:)
		end do
		
		close(10)
	
	
		deallocate(n,v,ones,stat=merror)
		
		print *,''
		print *,'Training completed.'
		print *,''
		
		call svm_end
			if (merror.NE.0) then
				print *,'Memory error finishing training'; call svm_end; end if
	else
		deallocate(n,v,ones,stat=merror)
			if (merror.NE.0) then
				print *,'Memory error finishing training'; call svm_end; end if
	end if
end subroutine training

function conjgrad(A,b)

	real(kind=8), dimension(:), intent(in)		::		b
    real(kind=8), dimension(size(b))			::		conjgrad
	real(kind=8), dimension(:), allocatable		::		x,p,r1,r2
	real(kind=8), dimension(:,:), intent(in)	::		A
	real(kind=8)								::		beta, lambda, tmp1, tmp2
	integer										::		i,k

	allocate(x(ns),p(ns),r1(ns),r2(ns),stat=merror)
	if (merror.NE.0) then
		print *,"Memory error. Conjugate gradient method could not be runned"; call svm_end; end if
	
	i = 0
	x = 0.d0
	r1 = b
	r2 = 0.d0
	p = 0.d0

	do while ((dot_product(r1,r1).GT.tol_cj).OR.(i.GT.10000))

		i = i + 1
		
		if (i==1) then
			p = b
		else
			beta = dot_product(r1,r1)/dot_product(r2,r2)
			do k = 1,ns
				p(k) = r1(k) + beta*p(k)
			end do
		end  if
	
		tmp1 = dot_product(r1,r1)

		tmp2 = dot_product(matmul(p,A),p)

		lambda = tmp1 / tmp2 !dot_product(r1,r1)/dot_product(matmul(p,A),p)
			
		do k = 1,ns
			x(k) = x(k) + lambda*p(k)
		end do

		r2 = r1
		r1 = r1 - lambda*matmul(A,p)
		
	end  do
		
	conjgrad = x

	deallocate(x,p,r1,r2)

end function conjgrad



end module flssvm_train
