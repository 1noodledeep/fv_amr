import numpy as np
from matplotlib import pyplot as plt

def max_first_deriv(soln, dx):
	fd = -1000000
	row, col = soln.shape
	# since there are ghost cells
	# on either side of the domain
	# we are computing the valid range
	for i in range(1, r-1):
		tmp = soln[i+1] - soln[i-1]
		if(abs(tmp)>fd):fd=abs(tmp)
	return fd/(2*dx)

def blip(soln,dx):
	pdf = (np.roll(soln, -1) - np.roll(soln,1))[1:-1]/dx*0.5
	wt = np.sum(pdf)*dx
	
	return pdf/wt

def width(pdf, x, dx):
	fmmt = np.sum(x*pdf)*dx
	smmt = np.sum(x*x*pdf)*dx

	return fmmt, np.sqrt(smmt - fmmt*fmmt)

if __name__ == "__main__":
	# solution on uniform coarse grid
	coarse = np.loadtxt("coarse.dat")
	# solution on uniform fine grid
	fine = np.loadtxt("fine.dat")
	# solution on coarse grid, with refinement in [.4,.6]
	center = np.loadtxt("refine_center.dat")
	# solution with mesh as above, but no reflux
	no_reflux = np.loadtxt("no_reflux.dat")

	dx = coarse[2,0] - coarse[1,0]
	dt = 0.07		# hard code dt
	x = coarse[1:-1, 0]

	N = 100
	pcoarse = blip(coarse[:,N], dx)
	pfine = blip(fine[:,N], dx)
	pcenter = blip(center[:,N], dx)
	pnreflux = blip(no_reflux[:,N], dx)

	print "At t = ",(N-2)*dt,  " width of front on"
	print "Coarse  grid (ncell=100) mean, width: ", width(pcoarse,x,dx)
	print "Fine    grid (ncell=200) mean, width: ", width(pfine,x,dx)
	print "Refined grid [.4,.6] w   reflux mean, width: ", width(pcenter,x,dx)
	print "Refined grid [.4,.6] w/o reflux mean, width: ", width(pnreflux,x,dx)

	plt.plot(x, pcoarse, label='Coarse')
	plt.plot(x, pfine, label='Fine')
	plt.plot(x, pcenter, label='Refined, reflux')
	plt.plot(x, pnreflux, 'c--', label='Refined, no reflux')
	plt.legend(loc='best')
	plt.show()
