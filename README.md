# Drift-Diffusion_models

Here are 1D, 2D, and 3D models which solve the semiconductor Poisson-Drift-Diffusion equations using finite-differences. These models can be used to model most semiconductor devices. The "Two-charge-carriers" versions of the models currently solve for a solar cell under illumination. The "Single-charge-carrier" versions solve for the current-voltage curve of a material which only has holes as the free carrier and is under a varying applied voltage in the dark.  All of the models can be modified to solve other systems (i.e. through changing the boundary conditions, adding recombination rates, and modifying the generation rate). 

The equations are solved using the self-consistent iterative approach called the Gummel method. In order to ensure numerical stability for the continuity equations, Scharfetter Gummel discretization as well as linear mixing of old and new solutions is used. The 1D/Drift-Diffusion/Single-charge-carrier/src folder also contains an implementation using Slotboom variables, which is an alternative way to achieve stability without using Scharfetter Gummel discretization.


--------------------------------------------------------------------------------------------------------------
Requirements for C++ implementations:
1D version: A C++11 compiler. A make file for the g++ compiler is included as well as a .pro file which can be used for compiling through the IDE QT Creator. Also the input files: "parameters.inp" and "gen_rate.inp" (for the 2 carrier versions) need to be located in the same directory as the source code.

2D and 3D versions: In addition to a C++11 compiler and input files, the linear algebra library Eigen is needed. The library is 
composed of only header files, so no linking is necessary. One just needs to specify the directory of the library files in the include. Eigen can be downloaded at: http://eigen.tuxfamily.org/index.php?title=Main_Page
Also include the openmp compiler flag to allow for Eigen to parallelize the matrix solving.

------------------------------------------------
The Matlab implementations only require Matlab.
