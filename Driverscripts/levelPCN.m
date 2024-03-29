function levelPCN(iter,stepsize,noiselevel,x0seed,noiseseed,sampleseed,std,delta,N)
%% Initialize forward mesh
fine_hmax = 0.01;
hmax = 0.0175;
meshpar = mesh_comp(hmax);
close all;
%% Load data
filename = strcat('noiseseed_',num2str(noiseseed),'_',num2str(fine_hmax),'_',num2str(hmax),'_level_',num2str(noiselevel),'_',num2str(N), '.mat');
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
q = 5^2;
alpha = 1.2;
tau = 10;
% should be 10
maxfreq = 4;


priorpar = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% For levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -1 1 inf]; % contour levels (ninterface+1)
priorpar.v = [0.3 0.1 0.5]; % values at each interface
priorpar.delta = delta; % smoothing factor
priorpar.dim = [priorpar.M,1];
priorpar.type = 'level';
priorpar.std = std;

%% Initialize forward model
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');
trunc = (2*datapar.N+1)*datapar.N;
fmdl = computeProjectionMatrices_coarse(fmdl,meshpar,datapar.meshpar_fine,priorpar,trunc);

%s = load('Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat');
%E_coarse = s.E_coarse;
%E_coarse = E_coarse(:,1:trunc);
%fmdl.U_proj_coarse = fmdl.Carea*E_coarse;
%s = load('Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat');
%E = s.E_coarse;
%U_proj = fmdl.Carea*E;
%fmdl.U_proj = U_proj;

%% Start guess for sampler
x0 = zeros(priorpar.M,1);
x0(3) = 2;
%% Sampling parameters
samplerpar.N_iter = iter;
samplerpar.jump_size = stepsize;
samplerpar.x0 = x0;

%% Filename
% todays date
datenow = datestr(now,'ddmm');
result_filename = strcat('x0seed_',num2str(x0seed),'_','noiseseed_',num2str(noiseseed),'_',num2str(fine_hmax),'_',num2str(hmax),'_level_',num2str(noiselevel),'_',datenow, '.mat');

%% Sample
%sampleseed = noiseseed;
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
if not(isfolder('Results/Level'))
    mkdir('Results/Level');
end

%% Save
save(strcat('Results/Level/xr_',result_filename),'results')
