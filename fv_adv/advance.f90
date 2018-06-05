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
	real 					:: left_out, right_in
	integer					:: i, j, w

	time = time+dt
	right_in = 0
	left_out = 0
	w = xor(s,1)

	call fill_ghost_cells(s)
	! Update the coarse level cell averages using
	! one sided upwind method
	!!! write(*,*) "Time: ", time
	do j=pl, pr
		soln(j, w) =cfl*soln(j-1, s) + (1-cfl)*soln(j, s)
	!!!	write(*,*) j, &
	!!!			"coarse", soln(j-1,s), soln(j, s), &
	!!!			"updated", soln(j,w)

	end do
	! Refinement: fill ghost cell and update fine level
	! Then average the fine level values onto coarse level
	call fine_step(dt, s, left_out, right_in)
	! Reflux: replace coarse flux with fine flux
	! Due to upwinding, only right edge need reflux
	soln(fhi+1, w) = soln(fhi+1, w) + cfl*(-soln(fhi, s) + right_in)

	s = w

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

	! Keeping cfl number constant, we need to sub cycle
	! assume dt_fine=0.5*dt_coarse
	w = xor(s,1)
	call fill_fine_ghosts(s)
	call fill_fine_ghosts(w)

	! Now take a step forward on the fine grid
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = cfl*fine_soln(i-1, s)&
						+ (1-cfl)*fine_soln(i,s)
	end do

	do i = 2*flo, 2*fhi+1
		fine_soln(i,s) = cfl*fine_soln(i-1, w)&
						+ (1-cfl)*fine_soln(i,w)
	end do

	! Transfer the solution over to the correct array
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = fine_soln(i,s)
	end do

	do i = flo, fhi
		!!!write(*,*) i, "fine ", fine_soln(2*i, w), fine_soln(2*i+1, w),&
		!!!	"coarse", soln(i, w), 0.5*(fine_soln(2*i, w) + fine_soln(2*i+1, w))
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
