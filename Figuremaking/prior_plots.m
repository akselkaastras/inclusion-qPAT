%% Plot of prior samples
%% Make prior for levelset with matern covariance
clc;
clear;
addpath(genpath(pwd));
hmax = 0.01;
meshpar = mesh_comp(hmax);
close all;
%% Make tight figure
figure(1);
ha = tight_subplot(1,4,[.01 .01],[.23 .07],[.05 .01]);
%% Smooth level sets
rng(5)

q = 5^2;
alpha = 1.2;
tau = 10;
maxfreq = 4;

std = 1;
delta = 0.1;

priorpar = prior_init_2d(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% For levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -1 1 inf]; % contour levels (ninterface+1)
priorpar.v = [0.3 0.1 0.5]; % values at each interface
priorpar.delta = delta; % smoothing factor
priorpar.dim = [priorpar.M,1];
priorpar.type = 'level';
priorpar.std = std;


xi = randn(priorpar.M,1);
xi = priorpar.lambdahalf.*xi;
theta = priorsample(xi,priorpar);

axes(ha(4))
plot_from_coef_level_smooth(xi,priorpar);

axes(ha(3))
plot_from_coef_2D(xi,priorpar)

%% star-shape
rng(5)
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
priorpar.std = 0.2;


% We use M Fourier coefficients and 2 numbers for each inclusion
xi = priorpar.std*randn(priorpar.dim);
xi = priorpar.lambdahalf.*xi;
%xi_center = 0.5*xi_total(priorpar.M+1:end,:);
xi_center = priorpar.center;
% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;


axes(ha(2))
plot_from_coef_star(xi,priorpar)
axes(ha(1))
c = parula(256);
%plot(xq,theta(:,1)-priorpar.mean,'color',c(128,:),'linewidth',5)
plot(xq,theta(:,1),'color',c(128,:),'linewidth',5)
hold on
%plot(xq,theta(:,2)-priorpar.mean,'color',c(end,:),'linewidth',5)
plot(xq,theta(:,2),'color',c(end,:),'linewidth',5)
xticks([0,pi,2*pi])
yticks(priorpar.mean + [-0.4,0,0.4])
xlim([0,2*pi])
set(gca,'XTickLabel',{'0','\pi','2\pi'})
xlabel('$\vartheta$','Interpreter','latex','FontSize',15)
axis square
ax = gca;
ax.FontSize = 20; 
set(gca, 'Color', [0.8,0.8,0.8] )
% parula colormap

%% Adjust size of figure
set(gcf, 'Position',  [100, 100, 1250, 360])
set(gcf,'color','w');

%% Export figure
%export_fig 'Figures/data_plot/prior_4by1_1.eps' -eps
print('Figures/data_plot/prior_4by1_1.eps','-depsc2')