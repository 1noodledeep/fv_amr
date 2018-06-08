import numpy as np
from matplotlib import pyplot as plt

def interpolate(fine):
	r,c = fine.shape
	nr = int((r-2)/2)

	# Create an array but leave it empty
	tmp = np.empty([nr, c])
	# Populate each row, one row at a time
	# By taking averages of fine cells onto coarse grid
	for i in range(nr):
		tmp[i] = 0.5*(fine[1+2*i] + fine[2+2*i])
	return tmp

# The next section only gets run when
# the script is executed as a main function
# Won't run if loaded as module
if __name__=="__main__":
	fine = np.loadtxt("soln_fine.dat")
	reflux = np.loadtxt("soln_reflux_const.dat")
	no_reflux = np.loadtxt("soln_no_reflux.dat")
	lin_interp = np.loadtxt("soln_lin_interp.dat")
	quad_interp = np.loadtxt("soln_quad_interp.dat")
	quad_impl =np.loadtxt("soln_impl.dat")

	# interpolate the fine solution (linear) onto the coarse grid
	# for comparison
	interp_fine = interpolate(fine)

	# start plotting
	plt.plot(reflux[1:-1,0], reflux[1:-1,1], 'b-', lw=2, label="Initial condition")
	plt.legend(loc='best')
	plt.ylim([0,1.5])
	plt.show()
	
	plt.plot(reflux[1:-1,0], interp_fine[:,2], 'g-', lw=2, label='Fine ncells=200')
	plt.plot(reflux[1:-1,0], reflux[1:-1,2], 'b--', lw=2, label='Locally refined [.3,.7]')
	plt.plot(reflux[1:-1,0], no_reflux[1:-1,2], 'c-', lw=2, label='Locally refined, no reflux')
	plt.plot(reflux[1:-1,0], quad_impl[1:-1,2], 'r.', label="Implicit, locally refine")

	# other plotting configs
	plt.ylim([-0.1, 1.])
	plt.legend(loc='best')
	plt.title("1D diffusion equation comparison")
	plt.show()


	# get the different
	plt.title("Difference compared to fine solution")
	diff = interp_fine - no_reflux[1:-1,:]
	plt.plot(reflux[1:-1,0], diff[:,2], 'r-', lw=2, label="Diff w no reflux")

	diff = interp_fine - reflux[1:-1,:]
	plt.plot(reflux[1:-1,0], diff[:,2], 'k-', lw=2, label="Diff w const interp")

	diff = interp_fine - lin_interp[1:-1,:]
	plt.plot(reflux[1:-1,0], diff[:,2], 'g-', lw=2, label="Diff w linear interpolation")

	diff = interp_fine - quad_interp[1:-1,:]
	plt.plot(reflux[1:-1,0], diff[:,2], 'c--', lw=2, label="Diff w quadratic interpolation")

	diff = interp_fine - quad_impl[1:-1,:]
	plt.plot(reflux[1:-1,0], diff[:,2], 'r--', lw=2, label="Diff w quad impl")

	plt.legend(loc='best')
	plt.show()
