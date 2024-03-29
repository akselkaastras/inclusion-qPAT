function starPCN(iter,stepsize,noiselevel,x0seed,noiseseed)
%% Initialize forward mesh
fine_hmax = 0.01;
hmax = 0.0175;
meshpar = mesh_comp(hmax);
close all;
%% Load data
filename = strcat('noiseseed_',num2str(noiseseed),'_',num2str(fine_hmax),'_',num2str(hmax),'_star_',num2str(noiselevel), '.mat');
filename_data = strcat('Data/data/data_',filename);
if isfile(filename_data)
    s = load(filename_data);
else
    s1 = 'Error: Please specify noiseseed and noiselevel which exists in data folder.';
    s2 = ' If it does not exist, then compute the data using the function dataSave.';
    s3 = ' No such data.';
    error(strcat(s1,s2,s3));
end
datapar = s.datapar;

%% Setup prior
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

%% Initialize forward model
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');
s = load('Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat');
E = s.E_coarse;
U_proj = fmdl.Carea*E;
fmdl.U_proj = U_proj;

%% Start guess for sampler
%rng(x0seed);
%x0 = priorpar.std*randn(priorpar.dim);
x0 = zeros(priorpar.dim);
x0(1,:) = [1 1];
%% Sampling parameters
samplerpar.N_iter = iter;
samplerpar.jump_size = stepsize;
samplerpar.x0 = x0;

%% Filename
result_filename = strcat('x0seed_',num2str(x0seed),'_','noiseseed_',num2str(noiseseed),'_',num2str(fine_hmax),'_',num2str(hmax),'_star_',num2str(noiselevel), '.mat');

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
