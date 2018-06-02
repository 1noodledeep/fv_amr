subroutine coarse_setup
!-------------------------------------------------
! Read in the system parameters and
! allocate solution array
!-------------------------------------------------
	use field
	implicit none

	! local variables
	integer :: a_stat
	! number of ghost cells on each side
	integer, parameter :: ngc=1
	write(*,'(A)',advance="NO") "The number of valid cells to use : "
	read(*, *) n

	! add ghost cells
	n2 = n+2
	! figure out the indexing scheme, going from -1 to n
	! -1 and nth cell are ghost cells
	! doesn't have to be -1 to n, though
	lo = -1
	hi = n

	pl = lo+ngc
	pr = hi-ngc

	u = 0.1
	dx = (rhi-rlo) / real(n)
	max_dt = 0.5* dx / u

	allocate(soln(lo:hi, 0:1), stat=a_stat)
	if (a_stat /=0) stop "coarse_setup: allocation failed. "
end subroutine coarse_setup

subroutine dealloc
!-------------------------------------------------
! Deallocate the dynamically allocated arrays
!-------------------------------------------------
	use field
	implicit none

	integer :: da_stat
	deallocate(soln, stat=da_stat)
	if (da_stat /=0) stop "coarse_setup: deallocation failed. "

end subroutine dealloc

subroutine initialize
!-------------------------------------------------
! Initialize the old solution array with a left
! Riemann problem
!-------------------------------------------------
	use field
	implicit none
	real 			 	:: x
	integer 		 	:: i

	do i = lo,hi
		x = dx*(real(i)+0.5)
		soln(i, 1)=0
		if (x<riemann_boundary) then
			soln(i, 0) = riemann_val
		else
			soln(i, 0) = 0
		end if
	end do
end subroutine initialize
