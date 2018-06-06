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
	reflux = np.loadtxt("soln_reflux.dat")
	no_reflux = np.loadtxt("soln_no_reflux.dat")

	# interpolate the fine solution (linear) onto the coarse grid
	# for comparison
	interp_fine = interpolate(fine)

	# start plotting
	plt.plot(reflux[1:-1,0], interp_fine[:,2], 'g-', lw=2, label='Fine ncells=200')
	plt.plot(reflux[1:-1,0], reflux[1:-1,2], 'b--', lw=2, label='Locally refined [.3,.7]')
	plt.plot(reflux[1:-1,0], no_reflux[1:-1,2], 'c-', lw=2, label='Locally refined, no reflux')

	# get the different
	diff = interp_fine - reflux[1:-1,:]
	plt.plot(reflux[1:-1,0], diff[:,2], 'k-', lw=2, label="Diff between local refined and fine")
	diff = reflux[1:-1,:] - no_reflux[1:-1,:]
	plt.plot(reflux[1:-1,0], diff[:,2], 'r-', lw=2, label="Diff between reflux and no reflux")

	# other plotting configs
	plt.ylim([-0.1, 1.])
	plt.legend(loc='best')
	plt.title("1D diffusion equation comparison")
	plt.show()
