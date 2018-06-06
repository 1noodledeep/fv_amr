module field
!------------------------------------------------
!  Module containing dynamically allocated
!  arrays that holds the solutions,
!  one old soln array, and one new soln array
!------------------------------------------------
	implicit none
	! length of the array
	integer				:: n, n2
	! starting and ending index of
	! the coarse and fine array
	integer				:: lo, hi, flo, fhi
	! the first and last non-ghost node index
	integer				:: pl, pr
	! physical coordinates of lower and higher bounds
	real				:: rlo=0., rhi=1., dx, max_dt
	real, parameter		:: source = 1.
	real, parameter		:: src_bnry_left = 0.4, src_bnry_right = 0.6
	real, parameter		:: cfl=0.2
	! Diffusion constant, rescale to be 1
	real, parameter		:: D=1.0
	! the soln array
	real, allocatable	:: soln(:, :)
	! the refined solution array
	real, allocatable	:: fine_soln(:,:)
end module field
