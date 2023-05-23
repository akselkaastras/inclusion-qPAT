%%
close all;
clc;
clear;
addpath(genpath(pwd))

%% Make data

noiseseed = 1;
noiselevel = 0.02;
priortype = 'level';
x0seed = 1;

iter = 1e6;
stepsize = 0.01;

%% Dont touch
dataSaveid(noiseseed, noiselevel, priortype)

%% Initialize forward mesh

fine_hmax = 0.01;
hmax = 0.0175;
meshpar = mesh_comp(hmax);
close all;

%% Setup prior 'star'
q = 8;
alpha = 1;
tau = 5;
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
priorpar.type = 'id';
priorpar.center = [0.37,-0.43;-0.44,0.36];
priorpar.std = 0.2;

%% Setup prior 'level'

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
xi = xi.*priorpar.lambdahalf;
theta = priorsample(xi,priorpar);

% Push-forward
%gamma = push_forward_levelset2D(theta,priorpar);
gamma = push_forward_levelset2D_smooth(theta,priorpar);

figure(1);
plot_from_gamma(gamma,meshpar);

figure(2);
plot_from_coef_level_smooth(xi,priorpar);


%% Initialize forward model
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');

%% Start guess for sampler
rng(x0seed);

%% computing projection matrix
N = 15;
trunc = (2*N+1)*N;
M = computeMass(meshpar);
E = eigenbasisFEM(meshpar,trunc);
U = M*E;
%%
fmdl.U = U;
%%
%x0 = priorpar.std*randn(priorpar.dim);
x0 = zeros(priorpar.dim);
x0(1,:) = [1 1];
xi_center = priorpar.center;

% Make KL expansion
theta = priorsample(x0,priorpar);
theta = priorpar.mean+theta;

% push-forward
gamma = push_forward_star2D(xi_center,theta,datapar.meshpar.p(1,:)',datapar.meshpar.p(2,:)',priorpar);
plot_from_gamma(gamma,meshpar);

%% Sampling parameters
samplerpar.N_iter = 100000;
samplerpar.jump_size = 0.1;
samplerpar.x0 = x0;

%% Likelihood
ll = compute_log_likelihood_pcn_id(x0, datapar, priorpar, fmdl);
%% Filename
result_filename = strcat('id_',num2str(noiselevel),'.mat');

%% Sample
sampleseed = noiseseed;
rng(sampleseed)
tic;
[LL, N_reject, XR] = pCNsampler(datapar, samplerpar, priorpar, fmdl, x0);
T = toc;
results.LL = LL;
results.N_reject = N_reject;
results.XR = XR;
results.T = T;

%% Sample adaptively
sampleseed = noiseseed;
rng(sampleseed)
tic;
[LL, N_reject, XR] = ApCNsampler(datapar, samplerpar, priorpar, fmdl, x0);
T = toc;
results.LL = LL;
results.N_reject = N_reject;
results.XR = XR;
results.T = T;
%% Plot

figure;
plot_from_coef_star(reshape(XR(end-1000,:),priorpar.M,priorpar.ninclusions),priorpar);
figure;
plot_from_coef_star(reshape(XR(end-3000,:),priorpar.M,priorpar.ninclusions),priorpar);

%% Make folder
if not(isfolder('Results'))
    mkdir('Results');
end
if not(isfolder('Results/Star'))
    mkdir('Results/Star');
end

%% Save
save(strcat('Results/Star/xr_',result_filename),'results')

%% Load eigenfunctions from FEM matrix

filename = 'Data/noise_model/eigenv/eigen_0.0175';
m = 0;
for i = 1:55
    str = strcat(filename,'_',num2str(i),'.mat');
    s = load(str);
    m = m + size(s.x,2);
    i
end

%%
%%
results = burnthin(results,1e4,50);
%% Plot chains
plot_chains_star(results,5)
%% Autocorrelation
plot_autocorr_star(results,5);
%%
plot_from_coef_star_mean(results,priorpar);