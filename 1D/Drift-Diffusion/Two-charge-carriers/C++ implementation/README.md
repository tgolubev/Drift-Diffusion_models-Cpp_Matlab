# Drift-Diffusion_models

1D Drift Diffusion 

This solves the semiconductor drift-diffusion equations in 1D for both electrons and holes using finite differences. The equations are Poisson eqn, 
continuity equation, and drift-diffusion equation which are solved in a decoupled iterative method (Gummel method). Scharfetter-Gummel
discretization as well as linear mixing of old and new solutions is used to maintain stability.

The code was originally developed for solar cells, but can be used for any semiconductor devices.

------------------------------------------------------------------------------------------------------

Input parameters are specified in parameters.inp.

Other input files:
Photogeneration rate: The generation rate file (filename can be specified in parameters.inp) should contain just 1 column of generation rates corresponding to each mesh point in the device (except the end points, x = 0 and x=L). An example generation rate file for device of 300nm thickness is included.

Experimental JV curve (optional): File name can be specified in parameters.inp. File should contain 2 columns: 1st column are voltage values and 2nd are current values (in A/m^3).

------------------------------------------------------------------------------------------------------

Code can be compiled using the makefile or the QT Creator .pro project file (just open that and run within QT).

------------------------------------------------------------------------------------------------------

Version History:

Current version 3.0 (under development):
   - Add optional automatic fitting algorithm to fit the model to an experimental JV curve.

Past versions:
   v2.0     Object orientation. Parameter input from file. Speed improvements.
   v1.0     First working version. Parameters are in parameters.h. No object orientation.
