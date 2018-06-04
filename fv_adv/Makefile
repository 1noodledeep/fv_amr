gf=gfortran
cflags=
iflags=-I.
lflags=-L.
# define the object files
objs=field.o setup.o advance.o refine.o
execs=solve

all:
	$(MAKE) executables

executables: $(execs)

# compile the object files
%.o: %.f90
	$(gf) $(cflags) $(iflags) -c $<

solve: solve.f90 $(objs)
	$(gf) $(cflags) $(iflags) -o $@ $^ $(lflags)

test: test.f90
	 $(gf) $(cflags) $(iflags) -o $@ $^ $(lflags)

clean:
	rm -rf $(objs) $(execs) field.mod

.PHONY: clean executables
