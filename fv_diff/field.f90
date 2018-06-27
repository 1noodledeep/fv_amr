module field
!------------------------------------------------
!  Module containing dynamically allocated
!  arrays that holds the solutions,
!  one old soln array, and one new soln array
!------------------------------------------------
	implicit none
	! length of the array
	integer				:: n, n2
	integer, parameter	:: dp=kind(0.d0)
	! starting and ending index of
	! the coarse and fine array
	integer				:: lo, hi, flo, fhi
	! the first and last non-ghost node index
	integer				:: pl, pr
	! physical coordinates of lower and higher bounds
	real(dp)				:: rlo=0._dp, rhi=1._dp, dx, max_dt
	real(dp), parameter		:: source = 1._dp
	real(dp), parameter		:: src_bnry_left = 0.4_dp, &
									   src_bnry_right = 0.5_dp
	real(dp), parameter		:: cfl=40.0_dp
	! Diffusion constant, rescale to be 1
	real(dp), parameter		:: D=1.0_dp
	! the soln array
	real(dp), allocatable	:: soln(:, :)
	! the refined solution array
	real(dp), allocatable	:: fine_soln(:,:)
end module field

module coarse_lp
!--------------------------------------------
!  Module containing constants and arrays
!  needed for DGESV call on coarse grid
!--------------------------------------------
	use field
	implicit none
	integer, parameter		:: cnrhs =1
	real(dp), allocatable	:: ca(:,:), cb(:,:)
	integer, allocatable	:: cipiv(:)
	integer					:: cinfo
end module coarse_lp

module fine_lp
!--------------------------------------------
!  Module containing constants and arrays
!  needed for DGESV call on coarse grid
!--------------------------------------------
	use field
	implicit none
	integer, parameter		:: fnrhs =1
	real(dp), allocatable	:: fa(:,:), fb(:,:)
	integer, allocatable	:: fipiv(:)
	integer					:: finfo
end module fine_lp
