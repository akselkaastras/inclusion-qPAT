# Bayesian inclusion detection for quantitative photoacoustic tomography

 should run everything from root folder

## driver.m file

## Parametrizations
Main functions are:
`prior_init.m`
`prior_init2d.m`

returns a structure containing KL (Karhunen-Loeve) eigenfunctions $`\{e_k\}_{k\in \mathbb{Z}^m}`$ and eigenvalues $`\{\lambda_k\}_{k\in \mathbb{Z}^m}`$ truncated corresponding to a user specified max frequency `max_freq` for the Matern covariance.
Given an i.i.d Gaussian vector of same length as the number of eigenfunctions, priorsample.m computes the KL expansion
$$\theta = \sum_{k} \xi_k \sqrt{\lambda_k} e_k.$$


Parametrizations:
`push_forward_levelset2D_smooth.m`
`push_forward_star2D_interp.m`

the first uses continuous approximation of the heaviside function to filter $`\theta`$ .
the second uses `inpoly2.m`

## MCMC

## FEM
Builds on finite element routines in first order Lagrange basis and discontinuous Galerkin basis. Many of these routines build on code by Niko Hänninen, Tanja Tarvainen and Ville Kolehmainen (University of Eastern Finland).

## Driverscripts folder
### Functions
levelPCN
starDGPCN


driver scripts are scripts are scripts that submit batch jobs at DTUs distributed computing service. these will return erros unless on DTU gbar. Outputs and errors are served in the folder 'Output' and 'Error', respectively.

## Figuremaking
Builds on export_fig [link]
Uses tight_subplot
