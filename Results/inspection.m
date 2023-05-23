%% mesh
clc;
clear;
addpath(genpath(pwd))
%% Make mesh
hmax = 0.0175;
meshpar = mesh_comp(hmax);

%% Load datapar
load('/work3/akara/qPAT-level/Data/data/data_noiseseed_1_0.01_0.0175_star_0.02.mat')

%% Load results
load('/work3/akara/qPAT-level/Results/Star/xr_x0seed_1_noiseseed_1_0.01_0.0175_star_0.02.mat')
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

%%

xr = results.XR(end-2000,:);
xr = reshape(xr,priorpar.dim);
figure;
plot_from_coef_star(xr,priorpar)
xr = results.XR(end-5000,:);
xr = reshape(xr,priorpar.dim);
figure;
plot_from_coef_star(xr,priorpar)

%%
results = burnthin(results,1e5,1);
%% Plot chains
plot_chains_star(results,5)
%% Autocorrelation
plot_autocorr_star(results,5);