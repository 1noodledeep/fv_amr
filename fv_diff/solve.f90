program solve
	use field
	implicit none
	real(dp)				:: T, dt, last_dt, time=0._dp
	real(dp), allocatable	:: output(:,:), mass(:)
	character(len=1024) 	:: filename, fmt_str
	integer					:: i, j, s=0, nt, frames,&
							   incr, last_incr, co, impl
	! Set up the solution array and the problem parameters
	call setup

	! Initialize the solution
	call initialize

	write(*, '(a)', advance="no") "Implicit enter nonzero #: "
	read(*,*) impl
	write (*,'(a)', advance="no") "Minimal duration to solve to: "
	read (*, *) T
	write (*,'(a)', advance="no") "dt to use: "
	read (*, *) dt
	write(*, '(a)', advance="no") "Number of frames to write: "
	read (*, *) frames

	!T = 2*max_dt
	!dt = max_dt
	!frames = 2
	! Compute the number of timesteps to take
	if(dt> max_dt) dt = max_dt
	if(T<dt) T=dt
	nt = int(T/dt)
	last_dt = (T-nt*dt)

	! Write at least the beginning and the end frame
	! At most 1 frame/timestep
	if (frames < 2) then
		frames =2
	else if (frames > nt+1) then
		frames = nt +1
	end if

	! get the number of timesteps per output frame
	! last frame might have more steps
	incr = int(nt/ (frames-1))

	if(incr * (frames-1) < nt) then
		last_incr = incr + nt - incr * (frames-1)
	else
		last_incr =  incr
	end if

	! Store the initial snap shot
	allocate(output(lo:hi, 0:frames-1))
	allocate(mass(0:frames-1))
	do i=lo, hi
		output(i, 0) = soln(i, 0)
	end do
	call integrate_mass(0, mass(0))
	write(*,*) "[time = ", time," m= ", mass(0),"]"

	! Initialize everything else in output to zero
	do j=1, frames-1
		do i=lo, hi
			output(i, j) = 0
		end do
	end do

	open(unit = 7, file="soln.dat")
	write(7, *) "# ================================================"
	write(7, *) "# Total time = ", T
	write(7, *) "# Frames to write = ", frames
	write(7, *) "# dt = ", dt
	write(7,*)  "# last dt = ", last_dt
	write(7, *) "# Number of steps per frame = ", incr, "(", last_incr, ")"
	write(7, *) "# Current time = ", time
	write(7, *) "# ================================================"


	write(7, *) "# Using ", n, "cells"
	write(7, *) "# Non-ghost cell indices ", pl, pr
	write(7, *) "# dx = ", dx
	write(7, *) "# ================================================"

	write(7, '(a)', advance = "no") " #   x   "
	! step forward
	write(7, '(f8.5 x)', advance="no") time
	do j = 1, frames-1

		if (j<frames-1) then
			co = incr
		else
			co = last_incr
		end if

		if(impl>0) then
			do i=1,co
				call impl_step_forward(dt, time, s)
			end do

			if (j==frames-1) then
				call impl_step_forward(last_dt, time, s)
			end if
		else
			do i=1, co
				call step_forward(dt, time, s)
			end do

			if (j==frames-1) then
				call step_forward(last_dt, time, s)
			end if

		end if
		do i=lo, hi
			output(i, j) =  soln(i, s)
		end do
		call integrate_mass(s, mass(j))
		write(*,*) "[time = ", time," m= ", mass(j),"]"

		write(7, '(f8.5 x)', advance="no") time
		
	end do
	write(7,*)

	! print out the results
	write(7, '(a)', advance = "no") " #  mass "
	do i=0, frames-1
		write(7, '(f8.5 x)', advance="no") mass(i)
	end do
	write(7,*)

	do i = lo, hi
		write(7, '(f8.5 x)', advance="no") dx*(i+0.5)
		do j = 0, frames-1
			write(7,'(f8.5 x)', advance="no") output(i, j)
		end do
		write(7,*)
	end do
	close(7)

	! free the dynamically allocated memory
	deallocate(output, mass)
	call dealloc
	stop
end program solve
