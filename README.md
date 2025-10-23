# fullGLmodelLOD
Code used in the paper "A multiscale approach to the stationary Ginzburg-Landau equations of superconductivity"

by Christian Döding, Benjamin Dörich, and Patrick Henning

This MATLAB implementation computes minimizers (u,A) of the full Ginzburg-Landau energy with a LOD approximation in the order parameter u and a P2 Lagrange FEM approximation in the vector potential A. The energy is minimized by a energy adaptive discretization of the Sobolev gradient flow as described in the paper. All of the paper's numerical results can be reproduced using the following procedure. First, a reference minimizer of the GL energy is computed (Section 7.2). Then, minimizers on the desired meshes are computed for convergence verification (Section 7.3 and 7.4).

Compute reference minimizer (Section 7.2):

1. Compute A_star via "main_compute_Astar.m"

2. in the preamble of the file "main_compute_minimizer_level1.m", adjust and set the parameter for the desired experiment:
	- the model parameters of the Ginzburg-Landau model
	- the numerical discretization parameter for LOD, P2-FEM and L^2 gradient flow discretizations
	- the preferences for saving and plotting the results after the computation
	
3. run the file "main_compute_minimizer_level1.m"

4. repeat the procedure for further levels, i.e., "main_compute_minimizer_level2.m", "main_compute_minimizer_level3.m" and "main_compute_minimizer_level4.m".

Note:  In general, the file "main_compute_minimizer_level1.m" can be used to compute a minimizer in an arbitrary setting.

Compute minimizer for convergence verification (Section 7.3 and 7.4) from reference minimizer:

1. in the preamble of the file "main_compute_minimizer_from_reference", adjust and set the parameter for the desired experiment:
	- the model parameters of the Ginzburg-Landau model
	- the numerical discretization parameter for LOD, P2-FEM and L^2 gradient flow discretizations
	- usage of A_star (Section 7.3)
	- the preferences for saving the results after the computation
	
Note:  A reference minimizer from "main_compute_minimizer_level4.m" is required

2. run the file "main_compute_minimizer_from_reference".

The implementation is created and tested for MATLAB version R2023b. Parallel computing (optional) and plotting (optional) require the MATLAB Add-On-Toolboxes "Parallel Computing Toolbox" (version 23.2) and "Partial Differential Equation Toolbox" (version 23.2).
