!******************************************************************************
!     flssvm_main.f90 - Version 0.1 
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

program flssvm_main

use flssvm_globals
use flssvm_utilities
use flssvm_train
use flssvm_predict

implicit none
	
	call get_arguments
	
	call load_input_data

	if (train) then
		call training
	else
		call predict
	end if

contains

subroutine get_arguments

	integer					::		i, n, cac, op_type
	character(len=50)		::		argument,tmpc
	logical					::		pass

	regression = .true.
	pass = .false.
	train = .true.
	out_data = "output"
	input_data = "training"
	kernel_type = 2
	s1 = 1.0
	s2 = 1.0
	new_data = "predict"
	gamma = 0.5
	tol_cj = 1.d-5

	call get_command_argument(1,argument)
	if (trim(argument).EQ.'--help') then

		print *,''
		print *,'FLSSVM - Fortran Least Squares Support Vector Machine' 
		print *,'         currently only in regression mode'
		print *,''
		print *,'flssvm [parameters ]'
		print *,''
		print *,'   -v              version'
		print *,'   help            this menu'
		print *,''
		print *,'Parameters:'
		print *,''
		print *,'   type=1          type of svm: 0 = classification, 1=regression'
		print *,'   act=0           action. 0 = training, 1 = prediction'
		print *,'   in=training     training data file'
		print *,'   out=output      output file'
		print *,'   g=0.5           gamma'
		print *,'   k=2             Kernel type.'	
		print *,'                   k = 1 for polinomial kernel K(x,x'') = (<x,x''> + s1)^s2' 
		print *,'                   k = 2 for RBF kernel, K(x,x'') = exp(-1/s1*||x-x''||^2)'
		print *,'                   k = 3 for fourier kernel, '
		print *,'                   K(x,x'') = II (1-s1^2)/(2(1-2scos(xi-xj))+s^2)'
		print *,'   s1=1.0          first kernel parameter'
		print *,'   s2=1.0          second kernel parameter'
		print *,'   new=predict     data to predict'
		print *,'   tol_cj=1e-10    conjugate gradient tolerance'
		print *,''
		print *,''
		print *,''
		print *,'example: '
		print *,'./flssvm k=2 s1=0.2 g=10  '
		print *,'./flssvm act=1 in=output out=result '
		print *,''
		print *,'This software is free software, distribuided under the GNU license'
		print *,'Hope it''s usefull. Please feel free to contruibute'
		
		stop
				
	else
	
		cac = command_argument_count()

		if (cac.EQ.0) then
			print *,''
			print *,'Assuming default parameters....'
			print *,'press q to quit, any other key to continue'
			print *,'for help tipe ./flssvm --help'
			print *,''
			read *,tmpc
			if (tmpc.EQ.'q') stop
		end if
		
		do i=1,cac

			call get_command_argument(i,argument)
			n = scan(argument,'=')

			select case(trim(argument(1:n-1)))
				case('-v')
					print *,'flssvm version 0.1'
					print *,''
					stop
				case('type')
					argument = trim(argument(n+1:50)) 
					read(argument,'(I1)') op_type
					if (op_type.EQ.0) then
						regression = .FALSE.
					else if (op_type.EQ.1) then
						regression = .TRUE.
					else
						print *,'wrong type argument'
						stop
					end if
				case('act')
					argument = trim(argument(n+1:50)) 
					read(argument,'(I1)') op_type
					if (op_type.EQ.0) then
						train = .TRUE.
					else if (op_type.EQ.1) then
						train = .FALSE.
					else
						print *,'wrong act argument'
						stop
					end if
				case('out')
					out_data = trim(argument(n+1:50))
				case('in')
					input_data = trim(argument(n+1:50))
				case('k')
					argument = trim(argument(n+1:50)) 
					read(argument,'(I1)') kernel_type
					if ((kernel_type.NE.1).AND.(kernel_type.NE.2).AND.(kernel_type.NE.3)) then
						print *,'wrong kernel type'
						stop
					end if
				case('new')
					new_data = trim(argument(n+1:50))
				case('s1')
					argument = trim(argument(n+1:50))
					read(argument,'(f16.4)') s1
				case('s2')
					argument = trim(argument(n+1:50))
					read(argument,'(f16.4)') s2
				case('g')
					argument = trim(argument(n+1:50))
					read(argument,'(f16.4)') gamma
					if (gamma.LE.0.0) then
						print *,'gamma cannot be equal or less than cero'
						stop
					end if
				case('tol_cj')
					argument = trim(argument(n+1:50))
					read(argument,'(E16.4)') tol_cj
					if (tol_cj.LE.0.d0) then
						print *,'tolerance must be greater than cero'
						stop
					end if
			end select
			
		end do
		
		pass = .true.
		
			if ((kernel_type.EQ.1).AND.((s1.LT.0.0).OR.(s2.LT.0.d0).OR.(s2-nint(s2).GT.1.d-5))) then
				print *,'first kernel parameter must be greater than cero, '
				print *,'and second kernel parameter must be an integer greater or equal'
				print *,'to cero'
				pass = .false.
			end if
			
			if ((kernel_type.EQ.2).AND.(s1.LE.0.d0)) then
				print *,'kernel parameter must be greater than cero'
				pass = .false.
			end if
			
			if ((kernel_type.EQ.3).AND.((s1.LE.0.d0).OR.(s1.GT.1.0))) then
				print *,'kernel parameter must be greater than cero'
				print *,'and less than one'
				pass = .false.
			end if
			
			if (.not.pass) then
				print *,''
				print *,'At least one wrong argument. Quiting'
				print *,''
				print *,'for help about the usage tipe --help'
				print *,''
				stop
			end if
	end if
		

end subroutine get_arguments


end program flssvm_main
