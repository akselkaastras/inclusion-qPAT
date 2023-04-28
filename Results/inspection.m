%% mesh
clc;
clear;
addpath(genpath(pwd))
%% Make mesh
hmax = 0.02;
hmax2 = 0.015;
meshpar1 = mesh_comp(hmax);
meshpar2 = mesh_comp(hmax2);

%% Load datapar
load('/work3/akara/qPAT-level/Data/data/data_noiseseed_1_0.01_0.02_star_0.04.mat')
%% Prior
q = 10;
alpha = 1;
tau = 0.5;
maxfreq = 7;
xq = linspace(0,2*pi,256);

priorpar = prior_init(xq,alpha,tau,q,maxfreq);

% For star-shaped inclusion parametrization
priorpar.ninclusions = 2; % number of interfaces
priorpar.v = [0.2 0.4]; % values at each interface

priorpar.background = 0.1;
priorpar.angles = linspace(0,2*pi,256)';
priorpar.mean = -2;
priorpar.dim = [priorpar.M,priorpar.ninclusions];
priorpar.type = 'star';
priorpar.center = [0.37,-0.43;-0.44,0.36];
priorpar.std = 0.2;

%% forward model
fmdl1 = precomputeFEM(meshpar1);
fmdl1 = precomputeRHS(meshpar1,fmdl1,datapar.wfun,datapar.wfungrad);
fmdl1 = fixingD(meshpar1,fmdl1,datapar.D_coarse');
%%
fmdl2 = precomputeFEM(meshpar2);
fmdl2 = precomputeRHS(meshpar2,fmdl2,datapar.wfun,datapar.wfungrad);
%%
D = ones(length(meshpar2.p),1);
fmdl2 = fixingD(meshpar2,fmdl2,D');

%%
xr = results.XR(end,:);
xr = reshape(xr,priorpar.dim);
plot_from_coef_star(xr,priorpar)
tic;
for i = 1:1000
    ll = compute_log_likelihood_pcn_star(xr, datapar, priorpar, fmdl1);
    i
end
T_1 = toc;
%%
tic;
for i = 1:100
    xi_center = priorpar.center;
    
    % Make KL expansion
    theta = priorsample(xr,priorpar);
    theta = priorpar.mean+theta;
    
    % push-forward
    gamma = push_forward_star2D(xi_center,theta,meshpar2.p(1,:)',meshpar2.p(2,:)',priorpar);
    %plot_from_gamma(gamma,datapar.meshpar)
    
    % evaluate forward model
    u = evalFowardModel(fmdl2,meshpar2,gamma);
    i
end
T_2 = toc;
%%
z1 = [9580,5.5];
z2 = [38127,31.8];