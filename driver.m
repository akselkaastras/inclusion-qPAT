clc; clear;
addpath(genpath(pwd))
%% Make mesh
fine_hmax = 0.02;
hmax = 0.0175;

meshpar = mesh_comp(hmax);
%% Make prior for levelset with matern covariance
q = 3;
alpha = 1;
tau = 1;
maxfreq = 4;
delta = 0.2;

priorpar = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% For levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
priorpar.v = [0.3 0.1 0.5]; % values at each interface
priorpar.delta = delta; % smoothing factor
priorpar.dim = [priorpar.M,1];
priorpar.type = 'level';
priorpar.std = 1;

% Sample prior
xi = randn(priorpar.M,1);
theta = priorsample(xi,priorpar);

% Push-forward
%gamma = push_forward_levelset2D(theta,priorpar);
gamma = push_forward_levelset2D_smooth(theta,priorpar);

figure(1);
plot_from_gamma(gamma,meshpar);

figure(2);
plot_from_coef_level_smooth(xi,priorpar);

%% Make prior for star-shaped set with matern covariance
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
%%

% We use M Fourier coefficients and 2 numbers for each inclusion
xi = priorpar.std*randn(priorpar.dim);
%xi_center = 0.5*xi_total(priorpar.M+1:end,:);
xi_center = priorpar.center;


% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;

% Push-forward
gamma = push_forward_star2D(xi_center,theta,meshpar.p(1,:)',meshpar.p(2,:)',priorpar);

figure(1);
plot_from_gamma(gamma,meshpar);

plot_from_coef_star(xi,priorpar)

%% Simulate data

% Kite and cushion parametrizations
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.4]; 0.12*cushioncurve+[-0.4;0.4]];
values = [0.2,0.4,0.1];
noiselevel = 0.01;
seed = 0;

datapar = make_data(curves,values,noiselevel,fine_hmax,meshpar,seed); 

%% Compute approximation error as diagonal normal distribution
N = 100;
datapar = computeApproxError(datapar,N);
%datapar.sigmasq = 0;
%datapar.epssq = 1e-8;
datapar.epssq_approx = datapar.epssq;
%% Initialize forward model for sampling
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');

%% Test sampler
rng(1);
N_iter = 10000;
jump_size = 0.02;
x0 = priorpar.std*randn(priorpar.dim);
plot_from_coef_star(x0,priorpar);
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