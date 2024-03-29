function starDGPCN(iter,stepsize,noiselevel,x0seed,noiseseed,sampleseed,std,N)
%% Initialize forward mesh
fine_hmax = 0.01;
hmax = 0.0175;
meshpar = mesh_comp(hmax);
close all;
%% Load data
filename = strcat('noiseseed_',num2str(noiseseed),'_',num2str(fine_hmax),'_',num2str(hmax),'_starDG_',num2str(noiselevel),'_',num2str(N), '.mat');
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
%q = 30;
%alpha = 1.75;
%tau = 1;
%maxfreq = 12;
% q = 15;
% alpha = 1.75;
% tau = 0.5;
% maxfreq = 12;
% q = 30;
% alpha = 2;
% tau = 2;
% maxfreq = 12;
%q = 200;
%alpha = 2;
%tau = 2;
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
priorpar.std = std;

%% Initialize forward model
fmdl = precomputeFEM_DG(meshpar);
fmdl = precomputeRHS_DG(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');
trunc = (2*datapar.N+1)*datapar.N;
fmdl = computeProjectionMatrices_coarse(fmdl,meshpar,datapar.meshpar_fine,priorpar,trunc);

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
% todays date
datenow = datestr(now,'ddmm');
result_filename = strcat('x0seed_',num2str(x0seed),'_','noiseseed_',num2str(noiseseed),'_',num2str(fine_hmax),'_',num2str(hmax),'_starDG_',num2str(noiselevel),'_',datenow, '.mat');

%% Sample
rng(sampleseed)
tic;
[LL, N_reject, XR] = pCNsampler(datapar, samplerpar, priorpar, fmdl, x0);
T = toc;
results.LL = LL;
results.N_reject = N_reject;
results.XR = XR;
results.T = T;
results.priorpar = priorpar;
results.datapar = datapar;
results.samplerpar = samplerpar;

%% Make folder
if not(isfolder('Results'))
    mkdir('Results');
end
if not(isfolder('Results/StarDG'))
    mkdir('Results/StarDG');
end

%% Save

save(strcat('Results/StarDG/xr_',result_filename),'results')
