# Drift-Diffusion_models

Here are 1D, 2D, and 3D models which solve the semiconductor Poisson-Drift-Diffusion equations using finite-differences. These models can be used to model most semiconductor devices. The "Two-charge-carriers" versions of the models currently solve for a solar cell under illumination. The "Single-charge-carrier" versions solve for the current-voltage curve of a material which only has holes as the free carrier and is under a varying applied voltage in the dark.  All of the models can be modified to solve other systems (i.e. through changing the boundary conditions, adding recombination rates, and modifying the generation rate). 

The equations are solved using the self-consistent iterative approach called the Gummel method. In order to ensure numerical stability for the continuity equations, Scharfetter Gummel discretization as well as linear mixing of old and new solutions is used. The 1D/Drift-Diffusion/Single-charge-carrier/src folder also contains an implementation using Slotboom variables, which is an alternative way to achieve stability without using Scharfetter Gummel discretization.


