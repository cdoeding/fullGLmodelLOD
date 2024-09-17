# fullGLmodelLOD
Code used in the paper "A multiscale approach to the stationary Ginzburg-Landau equations of superconductivity"

by Christian Döding, Benjamin Dörich, and Patrick Henning

This MATLAB implementation computes minimizers (u,A) of the full Ginzburg-Landau energy with a LOD approximation in the order parameter u and a P2 Lagrange FEM approximation in the vector potential A. The energy is minimized by a discretization of the L^2 gradient flow as described in the paper. All numerical results of the paper can be reproduced using the following procedure, which computes a minimizer of the GL energy in a desired setting:

1. in the preamble of the file "main_compute_GL_minimizer.m", adjust and set
	- the model parameters of the desired Ginzburg-Landau model
	- the numerical discretization parameter for LOD, P2-FEM and L^2 gradient flow discretizations
	- the preferences for saving and plotting the results after the computation
	
2. set the initial values for the order parameter and the vector potential in the file "main_compute_GL_minimizer.m". It is possible to use a best approximation of a previously computed minimizer as the initial value. 

3. run the file "main_compute_GL_minimizer.m".

For further details see the comments in "main_compute_GL_minimizer.m". The implementation is created and tested for MATLAB version R2023b. Parallel computing (optional) and plotting (optional) require the MATLAB Add-On-Toolboxes "Parallel Computing Toolbox" (version 23.2) and "Partial Differential Equation Toolbox" (version 23.2).
