2D Drift Diffusion 

This solves the semiconductor drift-diffusion equations in 2D using finite differences.
The "Two-charge-carriers" versions of the models currently solve for a solar cell under illumination. The "Single-charge-carrier" versions solve for the current-voltage curve of a material which only has holes as the free carrier and is under a varying applied voltage in the dark.  All of the models can be modified to solve other systems (i.e. through changing the boundary conditions, adding recombination rates, and modifying the generation rate). 

The equations are Poisson eqn, 
continuity equation, and drift-diffusion equation which are solved in a decoupled iterative method (Gummel method). Scharfetter-Gummel
discretization as well as linear mixing of old and new solutions is used to maintain stability.


--------------------------------------------------------------------------------------------------------------
Requirements for C++ implementations:
In addition to a C++11 compiler, the linear algebra library Eigen is needed. The library is 
composed of only header files, so no linking is necessary. One just needs to specify the directory of the library files in the include. Eigen can be downloaded at: http://eigen.tuxfamily.org/index.php?title=Main_Page

Also the input files: "parameters.inp" and "gen_rate.inp" (for the 2 carrier versions) need to be located in the same directory as the source code.

Also include the openmp compiler flag to allow for Eigen to parallelize the matrix solving.

------------------------------------------------
The Matlab implementations only require Matlab.
