subroutine impl_step_forward(dt, time, s)
!-------------------------------
! Take 1 step forward
! Call LAPACK DGESV to solve the
! the linear system resulted
! from implicity timestepping
!-------------------------------
	use field
	implicit none
	! Fine or coarse - which would
	! determine the lin sys size
	real(dp), intent(in)	:: dt
	real(dp), intent(inout)	:: time
	integer, intent(inout)	:: s
	! local variables
	real(dp)				:: left_out, right_in, old, mass_c, mass_f, mass_reflux
	integer					:: i, j, w

	time = time+dt
	right_in = 0
	left_out = 0
	w = xor(s,1)

	! Fill the ghost cells with prd bc
	call fill_ghost_cells(s)

	! Take an step forward on the coarse grid
	call coarse_impl(s)

	! Fill the ghost cells on the fine grid
	! s array has the ghost values at t_{n+1/2}
	! Take 1 out of 2 half steps forward on the coarse grid
	call quad_interpolate(s,1)
	call fine_impl(s)

	! Get the fluxes
	left_out =-(fine_soln(2*flo,w) - fine_soln(2*flo-1,s)) 
	right_in =-(fine_soln(2*fhi+2,s) - fine_soln(2*fhi+1,w)) 
	
	! Take second half step forward
	call quad_interpolate(w,0)
	call fine_impl(w)

	! Get the rest of fluxes
	left_out =left_out-(fine_soln(2*flo,s) - fine_soln(2*flo-1, w))
	right_in =right_in-(fine_soln(2*fhi+2,w) - fine_soln(2*fhi+1, s))

	! Transfer the solution over to the correct array
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = fine_soln(i,s)
	end do

	! Compute coarse fluxes before averaging
	old = -(soln(flo,w)-soln(flo-1,w))
	left_out = (old-left_out)*cfl

	old = -(soln(fhi+1,w) - soln(fhi,w))
	right_in = (-old+right_in)*cfl

	! Average solution onto coarse grid
	do i = flo, fhi
		soln(i, w) = 0.5_dp*(fine_soln(2*i, w) + fine_soln(2*i+1, w))
	end do

	! Reflux: replace coarse flux with fine flux
	! Remember in implicit solve, don't overwrite solution
	! before taking derivative to get old fluxes
	soln(flo-1, w) = soln(flo-1,w) + left_out
	soln(fhi+1, w) = soln(fhi+1,w) + right_in

	s = w
end subroutine impl_step_forward


subroutine coarse_impl(s)
!-------------------------------
! Implicit solve on coarse level
! Periodic boundaries
!-------------------------------
	use field
	use coarse_lp
	implicit none
	integer, intent(in)		:: s
	integer					:: i,j

	! Note the index start at 1
	do i=1, n
		do j=1, n
			ca(i,j)=0
		end do

		ca(i,i) = 1._dp+ 2._dp*cfl
		if(i>1 .and. i<n) then
			ca(i,i-1) = -cfl
			ca(i,i+1) = -cfl
		else if(i==1) then
			ca(i,i+1) = -cfl
			ca(i,n) = -cfl
		else
			ca(i,i-1) = -cfl
			ca(i,1) = -cfl
		end if

	end do

	! Populate the source term with current soln
	do i =1, n
		cb(i,1) = soln(i-1, s)
	end do

	!external dgesv
	call dgesv(n,cnrhs,ca,n,cipiv,cb,n,cinfo)

	if(cinfo .gt. 0) then
		write(*,*) "dgesv return non-zero info. Something went wrong."
	end if

	! Copy the solution over to the current solution array
	do i=1, n
		soln(i-1, xor(s,1)) = cb(i,1)
	end do
end subroutine coarse_impl


subroutine  fine_impl(s)
!-------------------------------
! Implicit solve on fine level
! Boundaries values treated
! as part of the source term
!-------------------------------
	use field
	use fine_lp
	implicit none
	integer, intent(in)	:: s
	real(dp), parameter	:: fine_cfl = 1._dp*cfl
	integer				:: i, j, nf

	nf = 2*(fhi-flo+1)
	
	! Note the index start at 1
	do i=1, nf
		do j=1, nf
			fa(i,j)=0
		end do

		fa(i,i) = 1._dp+ 2._dp*fine_cfl
		if(i>1 .and. i<nf) then
			fa(i,i-1) = -fine_cfl
			fa(i,i+1) = -fine_cfl
		else if(i==1) then
			fa(i,i+1) = -fine_cfl
		else
			fa(i,i-1) = -fine_cfl
		end if
	end do

	! Populate the source term with current soln
	! At the boundary, ghost cell values become
	! part of the source term
	fb(1,1) = fine_soln(2*flo, s) + fine_soln(2*flo-1,s) * fine_cfl
	j = 2*flo-1
	do i =2, nf-1
		fb(i,1) = fine_soln(j+i, s)
	end do
	fb(nf,1) = fine_soln(2*fhi+1,s) + fine_soln(2*fhi+2,s) * fine_cfl

	call dgesv(nf,fnrhs,fa,nf,fipiv,fb,nf,finfo)

	if(finfo .gt. 0) then
		write(*,*) "dgesv return non-zero info. Something went wrong."
	end if

	! Copy the solution over to the current solution array
	do i=1, nf
		fine_soln(j+i, xor(s,1)) = fb(i,1)
	end do
end subroutine fine_impl
