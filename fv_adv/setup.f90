subroutine setup
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
	!n = 200
	! add ghost cells
	n2 = n+2
	! figure out the indexing scheme, going from -1 to n
	! -1 and nth cell are ghost cells
	! doesn't have to be -1 to n, though
	lo = -1
	hi = n

	pl = lo+ngc
	pr = hi-ngc

	dx = (rhi-rlo) / dble(n)
	max_dt = cfl* dx / u

	allocate(soln(lo:hi, 0:1), stat=a_stat)
	if (a_stat/=0) stop "setup: allocation of coarse foln failed. "

	! Allocate the refined region solution
	flo = pl + int(n*0.40)
	fhi = pl + int(n*0.60)

	if(fhi>pr) then
		fhi = pr
	end if

	if(flo < pl) then
		flo = pl
	end if

	! Pad each side with 1 ghost cell
	allocate(fine_soln(2*flo-ngc:2*fhi+1+ngc, 0:1), stat=a_stat)
	if(a_stat/=0) stop "setup: allocation of fine soln failed. "

end subroutine setup

subroutine dealloc
!-------------------------------------------------
! Deallocate the dynamically allocated arrays
!-------------------------------------------------
	use field
	implicit none

	integer :: da_stat
	deallocate(soln, stat=da_stat)
	if (da_stat /=0) stop "dealloc: deallocation failed. "
	deallocate(fine_soln, stat=da_stat)
	if (da_stat /=0) stop "dealloc: deallocation failed. "

end subroutine dealloc

subroutine initialize
!-------------------------------------------------
! Initialize the old solution array with a left
! Riemann problem
!-------------------------------------------------
	use field
	implicit none
	real(dp) 		 	:: x
	integer 		 	:: i

	do i = lo,hi
		x = dx*(dble(i)+0.5)
		soln(i, 1)=0
		if (x<riemann_boundary) then
			soln(i, 0) = riemann_val
		else
			soln(i, 0) = 0
		end if
	end do

	do i = 2*flo-1, 2*fhi+2
		x = 0.5*dx*(dble(i)+0.5)
		fine_soln(i,1)=0
		if(x<riemann_boundary) then
			fine_soln(i,0) = riemann_val
		else
			fine_soln(i,0) = 0
		end if
	end do

	call fill_ghost_cells(0)
	call fill_ghost_cells(1)
end subroutine initialize

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

