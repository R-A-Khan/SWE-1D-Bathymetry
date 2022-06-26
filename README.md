# SWE-1D-Bathymetry
 Reconstruction of ocean floor bathymetry for the 1-D shallow water equations using data assimilation of surface wave observations, with applications to tsunami modelling. Full description of model derivation and evaluation can be found at:

R. A. Khan & N. K.-R. Kevlahan (2021) Variational assimilation of surface wave data for bathymetry reconstruction. Part I: algorithm and test cases, Tellus A: Dynamic Meteorology and Oceanography, 73:1, 1-25, DOI: 10.1080/16000870.2021.1976907

## Usage
Adjust model parameters in 
* __run_data_assimil_bath.m__  
or
* __run_data_assimil_bath_CG_tikhonov_multnorm.m__ (if using Conjugate Gradient descent and Tikhonov Regularisation). 
 
Utilises the data assimilation algorithm in __data_assimil_bath.m__ or __data_assimil_bath_CG_tikhonov_multnorm.m__.

Verification of numerical solvers used in DA scheme done using a kappa convergence test in 
* __Kappa_test_bath.m__

Visualize results using 
* __plots.m__ 

Built using Matlab 2018a, including Optimization toolbox.

