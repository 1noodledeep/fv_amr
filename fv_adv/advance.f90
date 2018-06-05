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

	! Update the coarse level cell averages using
	! one sided upwind method
	do j=pl, pr
		soln(j, w) =cfl*soln(j-1, s) + (1-cfl)*soln(j, s)
	end do
	! Refinement: fill ghost cell and update fine level
	! Then average the fine level values onto coarse level
	call fine_step(dt, s, left_out, right_in)
	! Reflux: replace coarse flux with fine flux
	! Due to upwinding, only right edge need reflux
	soln(fhi+1, w) = soln(fhi+1, w) + cfl*(-soln(fhi, s) + right_in)

	s = w
end subroutine step_forward

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

	! Fill the current fine solution with coarse soln
	! at t_n
	call fill_fine_ghosts(s, 0)

	! Fill the next fine solution step with average
	! between t_n and t_{n+1} on coarse grid
	call fill_fine_ghosts(w, 1)

	! Now take 1 half step forward on the fine grid
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = cfl*fine_soln(i-1, s)&
						+ (1-cfl)*fine_soln(i,s)
	end do

	! Since we use upwinding scheme and assume u>0
	! The flux is the upwind cell average
	! We also must capture the flux before the array gets
	! reused by second half step update
	left_out = 0.5*(fine_soln(2*flo-1, s) + fine_soln(2*flo-1, w))
	right_in = 0.5*(fine_soln(2*fhi+1, s) + fine_soln(2*fhi+1, w))

	! Now take second half step forward
	do i = 2*flo, 2*fhi+1
		fine_soln(i,s) = cfl*fine_soln(i-1, w)&
						+ (1-cfl)*fine_soln(i,w)
	end do


	! Transfer the solution over to the correct array
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = fine_soln(i,s)
	end do

	! Average solution onto coarse grid
	do i = flo, fhi
		soln(i, w) = 0.5*(fine_soln(2*i, w) + fine_soln(2*i+1, w))
	end do
end subroutine fine_step

subroutine fill_fine_ghosts (s, ts)
!-------------------------------------------------
! Fill the ghost cells on the fine grid
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in)	:: s, ts

	! Even though we only need the left edge
	! We fill both ends anyway
	! In between two fine steps, we use
	! average of the coarse grid solution
	if(ts == 0) then
		fine_soln(2*flo-1, s) = soln(flo-1, s)
		fine_soln(2*fhi+2, s) = soln(fhi+1, s)
	else
		fine_soln(2*flo-1, s) = 0.5*(soln(flo-1,0) + soln(flo-1,1))
		fine_soln(2*fhi+2, s) = 0.5*(soln(fhi+1, 0) + soln(fhi+1, 1))
	end if
end subroutine fill_fine_ghosts
