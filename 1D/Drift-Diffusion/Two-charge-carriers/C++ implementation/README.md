# Drift-Diffusion_models

This solves the semiconductor drift-diffusion equations in 1D for both electrons and holes using finite differences. The equations are Poisson eqn, 
continuity equation, and drift-diffusion equation which are solved in a decoupled iterative method (Gummel method). Scharfetter-Gummel
discretization as well as linear mixing of old and new solutions is used to maintain stability.

The gen_rate.txt input file is needed for photogeneration.cpp, or can comment that section out and use an analytic expression or set it to zero if are not 
studying a device under illumination. The gen_rate.txt file should contain just 1 column of generation rates corresponding to each mesh point in the device (except the end points, x = 0 and x=L). An example generation rate file for device of 300nm thickness is included.

The code was originally developed for solar cells, but can be used for any semiconductor devices.

Code can be compiled using the makefile or the QT Creator .pro project file (just open that and run within QT).
