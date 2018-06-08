gf=gfortran
cflags=
iflags=-I.
lflags=-L.
lp_flags=-llapack -lblas
# define the object files
objs=field.o setup.o advance.o impl_advance.o
execs=solve

all:
	$(MAKE) executables

executables: $(execs)

# compile the object files
%.o: %.f90
	$(gf) $(cflags) $(iflags) -c $<

solve: solve.f90 $(objs)
	$(gf) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_flags)

test: test.f90
	 $(gf) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_flags)

clean:
	rm -rf $(objs) $(execs) field.mod coarse_lp.mod fine_lp.mod

.PHONY: clean executables
