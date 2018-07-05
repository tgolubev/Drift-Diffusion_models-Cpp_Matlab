# Drift-Diffusion_models

3D Drift Diffusion 

This solves the semiconductor drift-diffusion equations in 3D for both electrons and holes using finite differences. The equations are Poisson eqn, 
continuity equation, and drift-diffusion equation which are solved in a decoupled iterative method (Gummel method). Scharfetter-Gummel
discretization as well as linear mixing of old and new solutions is used to maintain stability.

The code was originally developed for solar cells, but can be used for any semiconductor devices.



