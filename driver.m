clc; clear;
addpath(genpath(pwd))
%% Make mesh
meshpar = mesh_comp(2);
%% Make forward model
% define source from wfun
wfun = @(x1,x2) 1+x1^2;
wfungrad = @(x1,x2) [2*x1 0*x2]; 
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
%% Make prior for levelset with matern covariance
q = 1e5;
alpha = 3;
tau = 10;
maxfreq = 4;

priorpar = prior_init(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% for levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
priorpar.v = [1 0 2]; % values at each interface

% sample prior
xi = randn(priorpar.M,1);
theta = priorsample(xi,priorpar);

% push-forward
gamma = push_forward_levelset2D(theta,priorpar);

figure(1);
plot_from_gamma(gamma,meshpar);

figure(2);
[X,Y] = meshgrid(linspace(-1,1,200));
xq = X(:);
yq = Y(:);
plot_from_coef(xi,xq,yq,priorpar);

%% Simulate data

% kite and cushion
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 

datapar = make_data([kitecurve; cushioncurve],[2,3],noiselevel,meshrefine); 
