subroutine step_forward (dt, time, s)
!-------------------------------------------------
! Carry out two steps forward, so the most
! up to date soluiton is again in the old array
!-------------------------------------------------
	use field
	implicit none
	real, intent(in)		:: dt
	real, intent(inout)		:: time
	integer, intent (inout)	:: s
	! local varaibles
	real 					:: udt_by_dx, left_out, right_in
	integer					:: i, j, w

	udt_by_dx = u*dt/dx
	right_in = 0
	left_out = 0
	w = xor(s,1)

	call fill_ghost_cells(s)
	! Update the coarse level cell averages using
	! one sided upwind method
	do j=pl, pr
		soln(j, w) =udt_by_dx*soln(j-1, s) + (1-udt_by_dx)*soln(j, s)
	end do
	! Refinement: fill ghost cell and update fine level
	! Then average the fine level values onto coarse level
	call fine_step(dt, s, left_out, right_in)
	! Reflux: replace coarse flux with fine flux
	! Due to upwinding, only right edge need reflux
	soln(fhi+1, w) = soln(fhi+1, w) + udt_by_dx*(-soln(fhi, s) + right_in)

	s = w
	time = time+dt

end subroutine step_forward

subroutine fill_ghost_cells (s)
!-------------------------------------------------
! Fill the ghost cells on either end
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in) :: s

	! fill the ghost nodes to enforce bc
	! Dirichlet on the left,
	! value = riemann problem amplitude
	! Homogeneous Neumann on the right
	soln(lo, s) = riemann_val
	soln(hi, s) = soln(hi-1, s)
end subroutine fill_ghost_cells

subroutine fill_fine_ghosts (s)
!-------------------------------------------------
! Fill the ghost cells on the fine grid
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in)	:: s

	! Even though we only need 
	fine_soln(2*flo-1, s) = soln(flo-1, s)
	fine_soln(2*fhi+2, s) = soln(fhi+1, s)
end subroutine fill_fine_ghosts

subroutine fine_step (dt, s, left_out, right_in)
!-------------------------------------------------
! Take a fine step on the given region
! (flo, fhi) is the low and high indices of cells.
! The computation in this subroutine would also
! return the flux out of the left coarse cell and
! the right coarse cell bordering the fine cells.
!-------------------------------------------------
	use field
	implicit none
	real, intent(in)		:: dt
	integer, intent(in)		:: s
	real, intent(out)		:: left_out, right_in
	! local variables
	integer					:: i,w
	! a factor of two comes from refinement
	real					:: udt_by_dx

	udt_by_dx = u*dt*2/dx
	call fill_fine_ghosts(s)

	! Now take a step forward on the fine grid
	w = xor(s,1)
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = udt_by_dx*fine_soln(i-1, s)&
						+ (1-udt_by_dx)*fine_soln(i,s)
	end do

	do i = flo, fhi
		soln(i, w) = 0.5*(fine_soln(2*i, w) + fine_soln(2*i+1, w))
	end do

	! Since we use upwinding scheme and assume u>0
	! The flux is the upwind cell average
	left_out = fine_soln(2*flo-1, s)
	right_in = fine_soln(2*fhi+1, s)

end subroutine fine_step

subroutine interpolate (more, less, val)
!-------------------------------------------------
! Linear interpolation, assume interpolation
! point is half at 1/4 or 3/4 between two coarser
! points
!-------------------------------------------------
	real, intent(in)	:: more,less
	real, intent(out)	:: val
	val = 0.25*less + 0.75*more
end subroutine interpolate
