clc; clear;
addpath(genpath(pwd))
%% Make mesh
fine_hmax = 0.01;
hmax = 0.0175;

meshpar = mesh_comp(hmax);
%% Make prior for levelset with matern covariance
q = 5^2;
alpha = 1.2;
tau = 10;
maxfreq = 3;


std = 1;
delta = 0.1;

priorpar = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% For levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -1 1 inf]; % contour levels (ninterface+1)
priorpar.v = [0.3 0.1 0.5]; % values at each interface
priorpar.delta = delta; % smoothing factor
priorpar.dim = [priorpar.M,1];
priorpar.type = 'level';
priorpar.std = std;

% Sample prior
xi = priorpar.std * randn(priorpar.M,1);
xi = priorpar.lambdahalf.*xi;
theta = priorsample(xi,priorpar);

% Push-forward
%gamma = push_forward_levelset2D(theta,priorpar);
%gamma = push_forward_levelset2D_smooth(theta,priorpar);

figure(1);
plot_from_coef_2D(xi,priorpar)
%plot_from_gamma(gamma,meshpar);

figure(2);
plot_from_coef_level_smooth(xi,priorpar);



%% Make prior for star-shaped set with matern covariance in DG basis
q = 1e3;
alpha = 2.5;
tau = 4;
maxfreq = 12;
xq = linspace(0,2*pi,512);

priorpar = prior_init(xq,alpha,tau,q,maxfreq);

% For star-shaped inclusion parametrization
priorpar.ninclusions = 2; % number of interfaces
priorpar.v = [0.2 0.4]; % values at each interface

priorpar.background = 0.1;
priorpar.angles = linspace(0,2*pi,512)';
priorpar.mean = -2;
priorpar.dim = [priorpar.M,priorpar.ninclusions];
priorpar.type = 'starDG';
priorpar.center = [0.37,-0.43;-0.44,0.36];
priorpar.std = 0.2;

%% Sample 

% We use M Fourier coefficients and 2 numbers for each inclusion
xi = priorpar.std*randn(priorpar.dim);
xi = priorpar.lambdahalf.*xi;
xi_center = priorpar.center;


% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;

% Push-forward
gamma = push_forward_star2D(xi_center,theta,meshpar.p(1,:)',meshpar.p(2,:)',priorpar);

figure(1);
%plot_from_gamma(gamma,meshpar);

plot_from_coef_star(xi,priorpar)

%% Initialize forward model for sampling DG
fmdl = precomputeFEM_DG(meshpar);
fmdl = precomputeRHS_DG(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');
fmdl = computeProjectionMatrices_coarse(fmdl,meshpar,priorpar);
%% Test sampler
rng(1);
N_iter = 10000;
jump_size = 0.01;
x0 = zeros(priorpar.dim);
x0(1,:) = [1 1];
plot_from_coef_star(x0,priorpar);
ll = compute_log_likelihood_pcn_starDG(x0, datapar, priorpar, fmdl);
%%
samplerpar.N_iter = N_iter;
samplerpar.jump_size = jump_size;
tic;
[LL, N_reject, XR] = pCNsampler(datapar, samplerpar, priorpar, fmdl,x0);
T= toc;
%%
figure;
%plot_from_coef_star(reshape(XR(end-50,:),priorpar.M,priorpar.ninclusions),priorpar);
figure;
plot_from_coef_star(reshape(XR(end,:),priorpar.M,priorpar.ninclusions),priorpar);
%plot(LL)

%% Adaptive sampling
samplerpar.N_iter = N_iter;
samplerpar.jump_size = jump_size;
tic;
[LL, N_reject, XR] = ApCNsampler(datapar, samplerpar, priorpar, fmdl,x0);
T= toc;