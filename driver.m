clc; clear;
addpath(genpath(pwd))
%% Make mesh
finemeshrefine = 3;
meshrefine = 2;
meshpar = mesh_comp(meshrefine);
%% Make prior for levelset with matern covariance
q = 3;
alpha = 1;
tau = 1;
maxfreq = 3;
delta = 0.2;

priorpar = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% For levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
priorpar.v = [1 0 2]; % values at each interface
priorpar.delta = delta; % smoothing factor

% Sample prior
xi = randn(priorpar.M,1);
theta = priorsample(xi,priorpar);

% Push-forward
%gamma = push_forward_levelset2D(theta,priorpar);
gamma = push_forward_levelset2D_smooth(theta,priorpar);

figure(1);
plot_from_gamma(gamma,meshpar);

figure(2);
[X,Y] = meshgrid(linspace(-1,1,200));
xq = X(:);
yq = Y(:);
plot_from_coef_level_smooth(xi,xq,yq,priorpar);

%% Make prior for star-shaped set with matern covariance
q = 10;
alpha = 1;
tau = 0.5;
maxfreq = 7;
xq = linspace(0,2*pi,256);

priorpar = prior_init(xq,alpha,tau,q,maxfreq);

% For star-shaped inclusion parametrization
priorpar.ninclusions = 2; % number of interfaces
priorpar.v = [1 2]; % values at each interface
priorpar.background = 0;
priorpar.angles = linspace(0,2*pi,256)';
priorpar.mean = -2;


% We use M Fourier coefficients and 2 numbers for each inclusion
xi_total = randn(priorpar.M+2,priorpar.ninclusions);
xi = 0.2*xi_total(1:priorpar.M,:);
xi_center = 0.5*xi_total(priorpar.M+1:end,:);
%xi_center = [0,0;0,0];

% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;

% Push-forward
gamma = push_forward_star2D(xi_center,theta,meshpar.p(1,:)',meshpar.p(2,:)',priorpar);

figure(1);
plot_from_gamma(gamma,meshpar);

figure(2);
[X,Y] = meshgrid(linspace(-1,1,200));
xq = X(:);
yq = Y(:);
plot_from_coef_star(xi_center,xi,xq,yq,priorpar)

%% Simulate data

% kite and cushion
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.3]; 0.12*cushioncurve+[-0.3;0.3]];
values = [2,4];
noiselevel = 0.2;
seed = 0;

datapar = make_data(curves,values,noiselevel,finemeshrefine,meshpar,seed); 

%% Compute approximation error as diagonal normal distribution
N = 10;
datapar = computeApproxError(datapar,N);

%% Initialize forward model
% define source from wfun
%sigma = 0.5;
%m = 0.5*[sqrt(2),sqrt(2)];
%wfun = @(x1,x2) 2*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
%wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 

fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,datapar.wfun,datapar.wfungrad);
fmdl = fixingD(meshpar,fmdl,datapar.D_coarse');

%% Test log-likelihood

for i = 1:100
    i
    xi = randn(priorpar.M,1);
    compute_log_likelihood_pcn_level(xi, datapar, priorpar, fmdl)
end