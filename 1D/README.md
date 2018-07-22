# Drift-Diffusion_models

1D Drift Diffusion 

The Drift-Diffusion folder contains 1 and 2 carrier models. The 2 carrier version solves for the current-voltage curve of a solar cell under illumination. The 1 carrier version solves for the current-voltage curve for a material whose only free carrier are holes under an applied voltage in the dark.

The 1D/Drift-Diffusion/Single-charge-carrier/src folder also contains an implementation using Slotboom variables, which is an alternative way to achieve stability without using Scharfetter Gummel discretization.

-----------------------------------------------------------------------------------------------------------------------
The Mott-Gurney_law folder uses the weighted essentially non-oscillatory method (WENO) to solve for for the current-voltage curve for a material whose only free carrier are holes under an applied voltage in the dark in the Mott-Gurney limit (no diffusion current).

The WENO implementation is based on original 1D wave eqn code by Manuel Diaz.

WENO reference: Jiang & Shu; Efficient Implementation of Weighted ENO Schemes JCP. vol 126, 202-228 (1996)

-----------------------------------------------------------------------------------------------------------------------
Sample results are located in the Benchmarks folders.
