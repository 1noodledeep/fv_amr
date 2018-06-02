module field
!------------------------------------------------
!  Module containing dynamically allocated
!  arrays that holds the solutions,
!  one old soln array, and one new soln array
!------------------------------------------------
	implicit none
	! length of the array
	integer				:: n, n2
	! starting and ending index of the array
	integer				:: lo, hi
	! the first and last non-ghost node index
	integer				:: pl, pr
	! physical coordinates of lower and higher bounds
	real				:: rlo=0., rhi=1., dx, max_dt, u
	real, parameter		:: riemann_val = 1.
	real, parameter		:: riemann_boundary = 0.1
	! the soln array
	real, allocatable	:: soln(:, :)
end module field

