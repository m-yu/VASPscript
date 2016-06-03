VASPscript
==========

Scripts/codes help analyzing VASP output data

An algorithm for Bader Charge Integration
- J. Chem. Phys. 134, 064111 (2011) 
We propose an efficient, accurate method to integrate the basins of attraction of a smooth function defined on a general discrete grid, and apply it to the Bader charge partitioning for the electron charge density. Starting with the evolution of trajectories in space following the gradient of charge density, we derive an expression for the fraction of space neighboring each grid point that flows to its neighbors. This serves as the basis to compute the fraction of each grid volume that belongs to a basin (Bader volume), and as a weight for the discrete integration of functions over the Bader volume. Compared with other grid-based algorithms, our approach is robust, more computationally efficient with linear computational effort, accurate, and has quadratic convergence. Moreover, it is straightforward to extend to non-uniform grids, such as from a mesh-refinement approach, and can be used to both identify basins of attraction of fixed points and integrate functions over the basins.

Source codes written in C and Fortran
- C code ./weight_C
- Fortran code ./weight_F
