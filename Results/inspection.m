%% mesh
clc;
clear;
addpath(genpath(pwd))
%% Make mesh
hmax = 0.015;
meshpar = mesh_comp(hmax);

%% Load datapar
load('/work3/akara/qPAT-level/Data/data/data_noiseseed_1_0.01_0.015_star_0.04.mat')

%% Load results
load('/work3/akara/qPAT-level/Results/Star/xr_x0seed_1_noiseseed_1_0.01_0.015_star_0.12.mat')
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

xr = results.XR(end-100,:);
xr = reshape(xr,priorpar.dim);
plot_from_coef_star(xr,priorpar)
