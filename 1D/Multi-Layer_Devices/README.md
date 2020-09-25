# 1D Drift Diffusion for perovskite device

* This solves the drift-diffusion equations for a perovskite device with structure:
PEDOT:PSS/perovskite/C60/BCP. 
* It includes the effects of Langevin and Shockley-Read-Hall (SRH) trap-assisted recombination
* Photogeneration rate is specified through an input file: gen_rate.txt. The generation rate
  can be calculated using an optical transfer matrix approach. A sample file is included.

This code was used to model experimental perovskite devices and results were published in:
T. Golubev, D. Liu, R. Lunt, P. Duxbury. Understanding the impact of C60 at the interface of perovskite solar cells via drift-diffusion modeling. AIP Advances 9, 035026 (2019)

