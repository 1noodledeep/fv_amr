subroutine step_forward (dt, time, s)
!-------------------------------------------------
! Carry 1 step forward
!-------------------------------------------------
	use field
	implicit none
	real(dp), intent(in)		:: dt
	real(dp), intent(inout)		:: time
	integer, intent (inout)		:: s
	! local varaibles
	real(dp) 					:: left_out, right_in, old
	integer						:: i, j, w

	time = time+dt
	right_in = 0
	left_out = 0
	w = xor(s,1)

	! Since we have prd bc, we need to sync the
	! boundary before each time step
	call fill_ghost_cells(s)

	! Update the coarse level cell averages using
	! 2nd order centered difference 2nd derivative
	do j=pl, pr
		soln(j, w) = (1._dp - 2_dp*cfl) * soln(j, s)&
					+ cfl*(soln(j+1, s)+soln(j-1,s))
	end do

	! Refinement: fill ghost cell and update fine level
	! Then average the fine level values onto coarse level
	call fine_step(s, left_out, right_in)

	! Reflux: replace coarse flux with fine flux
	! Remember flux is negative first derivative
	old = -(soln(flo,s)-soln(flo-1,s))
	soln(flo-1, w) = soln(flo-1,w) + cfl*(old - left_out)
	old = -(soln(fhi+1,s) - soln(fhi,s))
	soln(fhi+1, w) = soln(fhi+1,w) + cfl*(-old + right_in)

	s = w
end subroutine step_forward

subroutine fine_step (s, left_out, right_in)
!-------------------------------------------------
! Take a fine step on the given region
! (flo, fhi) is the low and high indices of cells.
! The computation in this subroutine would also
! return the flux out of the left coarse cell and
! the right coarse cell bordering the fine cells.
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in)			:: s
	real(dp), intent(out)		:: left_out, right_in
	real(dp)					:: fine_cfl
	! local variables
	integer						:: i,w

	! Keeping cfl number constant, we need to sub cycle
	! assume dt_fine=0.5*dt_coarse
	w = xor(s,1)

	! Fill the current fine solution with coarse soln
	! at t_n
	! Fill the next fine solution step with average
	! between t_n and t_{n+1} on coarse grid
	!call fill_fine_ghosts(s, 0)
	!call fill_fine_ghosts(w, 1)

	!call lin_interpolate(s,0)
	!call lin_interpolate(w,1)

	call quad_interpolate(s,0)
	call quad_interpolate(w,1)

	! Since the factor here is dt/dx**2
	! dt and dx both decrease by a factor of 2
	! hence the ratio increase by a factor of 2
	fine_cfl = 2*cfl

	! Now take 1 half step forward on the fine grid
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = (1._dp - 2*fine_cfl) * fine_soln(i, s) &
					+ fine_cfl*(fine_soln(i+1, s)+fine_soln(i-1,s))
	end do

	! We also must capture the flux before the array gets
	! reused by second half step update

	! Flux in this case is -D phi_x, but D is absorbed in CFL#
	left_out =-(fine_soln(2*flo,s) - fine_soln(2*flo-1, s) &
				+fine_soln(2*flo,w) - fine_soln(2*flo-1,w)) 
	right_in =-(fine_soln(2*fhi+2,s) - fine_soln(2*fhi+1, s) &
				+fine_soln(2*fhi+2,w) - fine_soln(2*fhi+1,w)) 

	! Now take second half step forward
	do i = 2*flo, 2*fhi+1
		fine_soln(i,s) = (1. - 2*fine_cfl) * fine_soln(i, w) &
					+ fine_cfl*(fine_soln(i+1, w)+fine_soln(i-1,w))
	end do


	! Transfer the solution over to the correct array
	do i = 2*flo, 2*fhi+1
		fine_soln(i,w) = fine_soln(i,s)
	end do

	! Average solution onto coarse grid
	do i = flo, fhi
		soln(i, w) = 0.5_dp*(fine_soln(2*i, w) + fine_soln(2*i+1, w))
	end do
end subroutine fine_step

subroutine fill_fine_ghosts (s, ts)
!-------------------------------------------------
! Fill the ghost cells on the fine grid
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in)	:: s, ts

	! We fill ghost cells at the end of the fine grid
	! In between two fine steps, we use
	! average of the coarse grid solution

	! 0 for beginning of t_n
	if(ts == 0) then
		fine_soln(2*flo-1, s) = soln(flo-1, s)
		fine_soln(2*fhi+2, s) = soln(fhi+1, s)
	! non-zero for in between time steps
	else
		fine_soln(2*flo-1, s) = 0.5_dp*(soln(flo-1,0) + soln(flo-1,1))
		fine_soln(2*fhi+2, s) = 0.5_dp*(soln(fhi+1, 0) + soln(fhi+1, 1))
	end if
end subroutine fill_fine_ghosts

subroutine lin_interpolate(s, ts)
!-------------------------------------------------
! Fill the ghost cells by linear interpolation
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in) :: s, ts

	
	! 0 for beginning of t_n
	if(ts == 0) then
		fine_soln(2*flo-1, s) = 0.75_dp*soln(flo-1, s)+0.25_dp*soln(flo,s)
		fine_soln(2*fhi+2, s) = 0.25_dp*soln(fhi, s)+0.75_dp*soln(fhi+1,s)
	! non-zero for in between time steps
	else
		fine_soln(2*flo-1, s) = 0.5_dp*(0.75_dp*soln(flo-1, 1)+0.25_dp*soln(flo,1)&
								+  0.75_dp*soln(flo-1, 0)+0.25_dp*soln(flo,0))
		fine_soln(2*fhi+2, s) = 0.5_dp*(0.25_dp*soln(fhi, 1)+0.75_dp*soln(fhi+1,1)&
								+ 0.25_dp*soln(fhi, 0)+0.75_dp*soln(fhi+1,0))
	end if

end subroutine lin_interpolate

subroutine quad_interpolate(s, ts)
!-------------------------------------------------
! Fill the ghost cells by quadratic interpolation
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in) 	:: s, ts
	real(dp)				:: alo,blo,clo,ahi,bhi,chi
	real(dp)				:: sixteenth_dxsp, quarter_dxsp

	sixteenth_dxsp = 0.0625_dp*dx*dx
	quarter_dxsp = 0.25_dp*dx
	! 0 for beginning of t_n
	if(ts == 0) then
		call fit_quadratic( soln(fhi, s), soln(fhi+1, s), soln(fhi+2, s),&
							ahi, bhi, chi)
		call fit_quadratic( soln(flo-2, s), soln(flo-1, s), soln(flo, s),&
							alo, blo, clo)
		fine_soln(2*flo-1, s) = alo*sixteenth_dxsp + blo*quarter_dxsp + clo
		fine_soln(2*fhi+2, s) = ahi*sixteenth_dxsp - bhi*quarter_dxsp +chi
	! non-zero for in between time steps
	else
		! Average in time at the left edge
		call fit_quadratic( soln(flo-2, 0), soln(flo-1, 0), soln(flo, 0),&
							alo, blo, clo)
		call fit_quadratic( soln(flo-2, 1), soln(flo-1, 1), soln(flo, 1),&
							ahi, bhi, chi)
		fine_soln(2*flo-1, s) = 0.5_dp*(alo*sixteenth_dxsp + blo*quarter_dxsp + clo&
								+  ahi*sixteenth_dxsp + bhi*quarter_dxsp +chi)

		! Now average in time at the right edge
		call fit_quadratic( soln(fhi, 0), soln(fhi+1, 0), soln(fhi+2, 0),&
							alo, blo, clo)
		call fit_quadratic( soln(fhi, 1), soln(fhi+1, 1), soln(fhi+2, 1),&
							ahi, bhi, chi)
		fine_soln(2*fhi+2, s) = 0.5_dp*(alo*sixteenth_dxsp - blo*quarter_dxsp + clo&
								+ ahi*sixteenth_dxsp - bhi*quarter_dxsp +chi)
	end if

end subroutine quad_interpolate

subroutine fit_quadratic(left, cent, right, a, b, c)
!-------------------------------------------------
! A helper function that takes in 3 data points
! and fits a quadratic to it
! Assume data points are dx apart
! and the center data point is at x=0
!-------------------------------------------------
	use field
	implicit none
	real(dp), intent(in)	:: left, cent, right
	real(dp), intent(out)	:: a, b, c

	c = cent
	b = 0.5_dp*(right-left)/dx
	a = 0.5_dp*(right+left - 2*cent)/dx/dx

end subroutine fit_quadratic

subroutine integrate_mass(s, m)
!-------------------------------------------------
! Integrate the current mass of the solution
!-------------------------------------------------
	use field
	implicit none
	integer, intent(in)				:: s
	real(dp), intent(out)	:: m
	integer							:: i
	real(dp)				:: mass

	mass = 0
	call fill_ghost_cells(0)
	call fill_ghost_cells(1)

	do i = pl,pr
		mass = mass+soln(i,s)
	end do

	m = mass*dx
end subroutine integrate_mass
