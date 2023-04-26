%% Plot
%% Make prior for levelset with matern covariance
clc;
clear;
addpath(genpath(pwd));
hmax = 0.01;
meshpar = mesh_comp(hmax);
close all;
%% Make tight figure
figure(1);
ha = tight_subplot(1,4,[.01 .03],[.12 .07],[.01 .01]);
%% Smooth level sets
rng(1)
q = 3;
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

axes(ha(2))
plot_from_coef_level_smooth(xi,priorpar);

axes(ha(1))
plot_from_coef_2D(xi,priorpar)

% star-shape
rng(1)
q = 7;
alpha = 1;
tau = 1;
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
% We use M Fourier coefficients and 2 numbers for each inclusion
xi = priorpar.std*randn(priorpar.dim);
%xi_center = 0.5*xi_total(priorpar.M+1:end,:);
xi_center = priorpar.center;
% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;


axes(ha(4))
plot_from_coef_star(xi,priorpar)
axes(ha(3))
c = parula(256);
plot(xq,theta(:,1)-priorpar.mean,'color',c(128,:),'linewidth',4)
hold on
plot(xq,theta(:,2)-priorpar.mean,'color',c(end,:),'linewidth',4)
xticks([0,pi,2*pi])
yticks([0,0.4])
xlim([0,2*pi])
set(gca,'XTickLabel',{'0','\pi','2\pi'})
ax = gca;
ax.FontSize = 16; 
% parula colormap

%% Adjust size of figure
set(gcf, 'Position',  [100, 100, 1320, 360])

%% Export figure
export_fig 'test.eps' -eps -transparent