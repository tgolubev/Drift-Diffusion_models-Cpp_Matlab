# Drift-Diffusion_models

Here are 1D, 2D, and 3D models which solve the semiconductor Poisson-Drift-Diffusion equations using finite-differences. These models can be used to model most semiconductor devices. The models currently solve for the current-voltage curve and charge densities for a solar cell under illumination but they can be easily modified to solve other systems (i.e. through changing the boundary conditions). 

The equations are solved using the self-consistent iterative approach called the Gummel method. In order to ensure numerical stability for the continuity equations, Scharfetter Gummel discretization is used.

