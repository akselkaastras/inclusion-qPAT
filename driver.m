%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% DRIVER SCRIPT FOR %%%%%%%%%
%%% qPAT inclusion detection code %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Aksel Kaastrup Rasmussen
% Technical University of Denmark
% 2023

clc; clear;
addpath(genpath(pwd))
%% Make mesh

% Fine mesh for data simulation
fine_hmax = 0.04;

% Coarse mesh for forward map computations
hmax = 0.05;

% Make coarse mesh
meshpar = mesh_comp(hmax);
%% Make prior for level set parametrization with Matern covariance
% Matern parameters:
q = 5^2;
alpha = 1.2;
tau = 10;
maxfreq = 3;

% Prior rescaling factor
std = 1;

% Initialize basis functions and eigenvalues on 2D square domain 
% and evaluated in mesh points of coarse mesh
priorpar = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% Level set parametrization parameters
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -1 1 inf]; % contour levels (ninterface+1)
priorpar.v = [0.3 0.1 0.5]; % values at each interface
priorpar.delta = 0.1; % smoothing factor of level set parametrization
priorpar.dim = [priorpar.M,1];
priorpar.type = 'level';
priorpar.std = std;

% Sample prior example
% Sample random coefficients
xi = priorpar.std * randn(priorpar.M,1);
xi = priorpar.lambdahalf.*xi;

% Make expansion
theta = priorsample(xi,priorpar);

% Push-forward through smoothened and non-smooth level set parametrization
gamma = push_forward_levelset2D(theta,priorpar);
gamma_smooth = push_forward_levelset2D_smooth(theta,priorpar);


% Make trisurf plots based on values of gamma and gamma_smooth in mesh
% nodes:
figure(1);
plot_from_gamma(gamma,meshpar)

figure(2);
plot_from_gamma(gamma_smooth,meshpar)



%% Make prior for star-shaped set parametrization with Matern covariance
% Since the star-shaped set parametrization gives discontinuous samples we 
% consider a discontinuous Galerkin basis (indicator functions over each
% element)

% Matern parameters:
q = 1e3;
alpha = 2.5;
tau = 4;
maxfreq = 12;

% Mesh points on [0,2pi] in which distance function is evaluated
xq = linspace(0,2*pi,512);

% Initialize basis functions and eigenvalues on [0,2pi] and evaluated in
% mesh points
priorpar = prior_init(xq,alpha,tau,q,maxfreq);

% Parameters for star-shaped inclusion parametrization
priorpar.ninclusions = 2; % number of inclusions
priorpar.v = [0.2 0.4]; % values at each interface
priorpar.background = 0.1; % background value
priorpar.angles = xq';
priorpar.mean = -2; % mean of distance function before taking exp(.)
priorpar.dim = [priorpar.M,priorpar.ninclusions];
priorpar.type = 'starDG';
priorpar.center = [0.37,-0.43;-0.44,0.36]; % center of inclusions
priorpar.std = 0.2; % prior scaling
 
% Sample prior example
% We use M Fourier coefficients and 2 numbers for each inclusion
% Sample random coefficients
xi = priorpar.std*randn(priorpar.dim);
xi = priorpar.lambdahalf.*xi;
xi_center = priorpar.center;

% Make expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;

% Push-forward through star shaped set parametrization
% Here, we evaluate the "inclusion" in meshpar.n points of each element to
% get an averaged value in those elements, where the curve passes through: 
gamma = push_forward_star2D_interp(xi_center,theta,meshpar.xq,meshpar.yq,meshpar.n,meshpar.HN,priorpar);

% Plot high-resolution parameter, and in the DG basis.
figure(1);
plot_from_coef_star(xi,priorpar)

figure(2);
pdesurf(meshpar.p,meshpar.t,gamma')
colormap default
view(2)

%% Make synthetic noisy data and compute approximation error
% To repeat experiment we set the seed:
noiseseed = 1;

% Relative noise level
noiselevel = 0.02;

% Make curves
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.4]; 0.12*cushioncurve+[-0.4;0.4]];
values = [0.2,0.4,0.1];


% Makes data
datapar = make_data(curves,values,noiselevel,fine_hmax,meshpar,noiseseed); 

% Computes approximation error and approximate by N(0,sigma^2 I) 
% based on N samples of the prior
N = 200;
datapar = computeApproxError(datapar,N,priorpar.type);


%% Initialize forward model for sampling

% Computing the number of Dirichlet eigenbasisfunctions to project to
trunc = datapar.N*(2*datapar.N+1);

if strcmpi(priorpar.type,'starDG')
    % Precomputing FEM matrices (independent of diffusion and absorption)
    fmdl = precomputeFEM_DG(meshpar);
    % Precomputing rhs of variational form (indep. of diff. and absorp.)
    fmdl = precomputeRHS_DG(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
    % Fixing diffusion parameter (which is assumed to be known)
    fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');
    % Computing projection matrices to project H = gamma*u to Dirichlet
    % eigenbasis
    fmdl = computeProjectionMatrices_coarse(fmdl,meshpar,datapar.meshpar_fine,priorpar,trunc);
elseif strcmpi(priorpar.type,'level')
    % Precomputing FEM matrices (independent of diffusion and absorption)
    fmdl = precomputeFEM(meshpar);
    % Precomputing rhs of variational form (indep. of diff. and absorp.)
    fmdl = precomputeRHS(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
    % Fixing diffusion parameter (which is assumed to be known)
    fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');
    % Computing projection matrices to project H = gamma*u to Dirichlet
    % eigenbasis
    fmdl = computeProjectionMatrices_coarse(fmdl,meshpar,datapar.meshpar_fine,priorpar,trunc);
end
%% Compute likelihood from a starting guess

if strcmpi(priorpar.type,'starDG')
    % Starting guess from perfect circle
    x0 = zeros(priorpar.dim);
    x0(1,:) = [1 1];
    
    % Plot
    plot_from_coef_star(x0,priorpar);
    
    % Compute likelihood
    ll = compute_log_likelihood_pcn_starDG(x0, datapar, priorpar, fmdl);
elseif strcmpi(priorpar.type,'level')
    % Starting guess from third basis function
    x0 = zeros(priorpar.M,1);
    x0(3) = 2;

    % Plot 
    plot_from_coef_level_smooth(x0,priorpar)

    % Compute likelihood
    ll = compute_log_likelihood_pcn_level(x0, datapar, priorpar, fmdl);
end
%% Markov Chain Monte Carlo (MCMC) with preconditioned Crank-Nicolson (pCN)
% For now we ignore approximation error:
datapar.sigmasq = 0;

% Set sampling seed
rng(1);

% Number of iterations
samplerpar.N_iter = 10000;

% pCN stepsize
if strcmpi(priorpar.type,'starDG')
    % In this case we can take a large step size
    samplerpar.jump_size = 0.03;
elseif strcmpi(priorpar.type,'level')
    % In this case we take a small step size
    samplerpar.jump_size = 0.002;
end

tic;
[LL, N_reject, XR] = pCNsampler(datapar, samplerpar, priorpar, fmdl,x0);
T= toc;
%% Plot
% Last sample for example
xr = XR(end,:);

figure;
if strcmpi(priorpar.type,'starDG')
    
    xr = reshape(xr,priorpar.dim);
    % Make expansion
    theta = priorsample(xr,priorpar);
    theta = priorpar.mean+theta;

    % Pushforward
    gamma = push_forward_star2D_interp(xi_center,theta,meshpar.xq,meshpar.yq,meshpar.n,meshpar.HN,priorpar);

    % Plot
    pdesurf(meshpar.p,meshpar.t,gamma')
    colormap default
    view(2)
elseif strcmpi(priorpar.type,'level')
    % Make expansion
    theta = priorsample(xr',priorpar);
    
    % Push-forward
    gamma = push_forward_levelset2D(theta,priorpar);
    gamma_smooth = push_forward_levelset2D_smooth(theta,priorpar);
    
    % Plot
    figure(1);
    plot_from_gamma(gamma,meshpar)
end
