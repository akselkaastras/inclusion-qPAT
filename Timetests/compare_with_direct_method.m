%% Compare with closed form solution without regularization

clc; clear;
addpath(genpath(pwd))
%% Make mesh
fine_hmax = 0.01;
hmax = 0.0175;

meshpar = mesh_comp(hmax);

%% Load data
load('/work3/akara/qPAT-level/Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat')
N = 13;
trunc = (2*N+1)*N;
v = datapar.bq'*E_coarse(:,1:trunc)';
plot_from_gamma(v,meshpar)

%% Initialize forward model for sampling DG
wfun = @(x1,x2) 0*x1+1;
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');

%% Solve
% System matrix A 
A = fmdl.K;

% Build rhs
Hmatrix = spdiags(v', 0, length(v), length(v));

Q2 = fmdl.L2*Hmatrix;
Q2 = sum(Q2,2);
%fmdl.Q1 = fmdl.Q1(fmdl.phi);
%Q2 = Q2(fmdl.phi);

Q =  - fmdl.Q1 - Q2;

% Solve
R = chol(A);
u = R\(R'\Q);
u = u(fmdl.r);

U = zeros(length(v),1);
U(datapar.meshpar.NZ) = u;
u = U + datapar.W;

gamma_recon = v'./u;

plot_from_gamma(gamma_recon,meshpar)


