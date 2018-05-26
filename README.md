# Drift-Diffusion_models

The Drift-Diffusion folder contains codes which solve the semiconductor Poisson-Drift-Diffusion equations for both 1 carrier and 2 carriers (electrons and holes).
Scharfetter Gummel discretization is used.

The Mott-Gurney_law folder solves the same system in the high field limit where the diffusion term is neglected using the 
5th order weighted essentially nonoscillatory method (WENO5). The WENO5 code was adapted from codes by Manuel Diaz for solving 1D wave equation.
