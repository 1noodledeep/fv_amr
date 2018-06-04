subroutine step_forward (dt, time, s)
!-------------------------------------------------
! Carry out two steps forward, so the most
! up to date soluiton is again in the old array
!-------------------------------------------------
	use field
	implicit none
	real, intent(in)		:: dt
	real, intent(inout)		:: time
	real 					:: udt_by_dx
	integer					:: i, j
	integer, intent (inout)	:: s

	udt_by_dx = u*dt/dx

	call fill_ghost_nodes(s)
	! Update the cell averages using
	! one sided upwind method
	do j=pl, pr
		soln(j, xor(s, 1)) =udt_by_dx*soln(j-1, s) + (1-udt_by_dx)*soln(j, s)
	end do

	s = xor(s,1)
	time = time+dt

end subroutine step_forward

subroutine fill_ghost_nodes (s)
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
end subroutine fill_ghost_nodes
