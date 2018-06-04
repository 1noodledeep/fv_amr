# fv_amr
A small code base to solve 1D advection equation with constant velocity (flow left to right).
We use finite volume method, and an implementation of mesh refinement in a fixed region.
The code is primarily written in Fortran 90, source code assume free form.


## Prerequisites

- GNU make

- GNU GCC compiler (using gfortran)

- GNUPLOT (visualization)

## Solve
To run the code, cd into fv_adv/ and do

> make

> ./solve

This would prompt user for

(1) the number of cells to use;

(2) the total time to solve for;

(3) the timestep size (maximum is limited by CFL condition);

(4) and the number of snapshots to output.

## Visualization
The snapshots include ghost cells, and it's output as a table in a text file.
Each column is the solution at a particular time step. 
The time is determined by the total simulation time and the number of snapshots user chooses to output.

In GNUPLOT, do

> plot for [i=2:num_frames+1] u 1:i with line

to visualize all the solution at once. num_frames is the number of snapshots.
