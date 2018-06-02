subroutine step_forward (dt, time)
!-------------------------------------------------
! Carry out two steps forward, so the most
! up to date soluiton is again in the old array
!-------------------------------------------------
	use field
	implicit none
	real, intent(in)	:: dt
	real, intent(inout)	:: time
	real 				:: udt_by_dx
	integer				:: i, j

	udt_by_dx = u*dt/dx

	do i=0,1
		call fill_ghost_nodes(i)
		! Update the cell averages using
		! one sided upwind method
		do j=pl, pr
			soln(j, xor(i, 1)) =udt_by_dx*soln(j-1, i) + (1-udt_by_dx)*soln(j, i)
		end do
	end do
	call fill_ghost_nodes(0)

	time = time+2*dt

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
