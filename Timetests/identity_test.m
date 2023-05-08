%%
close all;
clc;
clear;
addpath(genpath(pwd))

%% Make data

noiseseed = 1;
noiselevel = 0.2;
priortype = 'star';
x0seed = 1;

iter = 1e6;
stepsize = 0.01;

%% Dont touch
dataSaveid(noiseseed, noiselevel, priortype)

%% Initialize forward mesh

fine_hmax = 0.01;
hmax = 0.015;
meshpar = mesh_comp(hmax);
close all;
%% Load data
load('/work3/akara/qPAT-level/Data/data/data_id_0.2.mat')

%% Setup prior
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

%% Initialize forward model
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');

%% Start guess for sampler
rng(x0seed);
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
samplerpar.N_iter = iter;
samplerpar.jump_size = 0.05;
samplerpar.x0 = x0;

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

%% Make folder
if not(isfolder('Results'))
    mkdir('Results');
end
if not(isfolder('Results/Star'))
    mkdir('Results/Star');
end

%% Save
save(strcat('Results/Star/xr_',result_filename),'results')

%% Find eigenfunctions from FEM matrix
[K,C] = computeLaplacian(meshpar);
B = speye(size(K));
%%
meshpar.NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));

%%
MM = 10;
%[VV,D] = eigs(K,C',MM,'largestabs');
[VV,D] = pdeeig(K,B,C,[0,300]);
%% Add zeros
V = zeros(pN,size(VV,2));
for i = 1:size(VV,2)
    V(meshpar.NZ,i) = VV(:,i);
end
%%
plot_from_gamma(V(:,2),meshpar);