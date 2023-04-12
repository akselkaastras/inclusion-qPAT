%% Plot
%% Make prior for levelset with matern covariance

meshrefine = 2;
meshpar = mesh_comp(meshrefine);

%% Smooth level sets
q = 1;
alpha = 1;
tau = 1;
maxfreq = 4;
delta = 0.2;

priorpar = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% For levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
priorpar.v = [0.03 0.01 0.05]; % values at each interface
priorpar.delta = delta; % smoothing factor
priorpar.dim = [priorpar.M,1];
priorpar.type = 'level';

% Sample prior
xi = randn(priorpar.M,1);
theta = priorsample(xi,priorpar);

figure(2);
plot_from_coef_level_smooth(xi,priorpar);


