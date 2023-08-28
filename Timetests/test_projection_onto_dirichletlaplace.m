close all;
clc;
clear;
addpath(genpath(pwd))
%% Extract eigenbasis for Dirichlet Laplacian
hmax = 0.0175;
hmax_fine = 0.01;
meshpar = mesh_comp(hmax);
meshpar_fine = mesh_comp(hmax_fine);
close all;
N = 15;
trunc = (2*N+1)*N;
E = eigenbasisFEM(meshpar,trunc);
E_fine = eigenbasisFEM(meshpar_fine,trunc);

%% Compute mass matrix
M = computeMass(meshpar);
M_fine = computeMass(meshpar_fine);

%% Big Matrix product;

U = M*E;
U_fine = M_fine*E_fine;
%% How fast?
ninclusions = 2;
npoints = length(meshpar.p);
npoints_fine = length(meshpar_fine.p);
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.4]; 0.12*cushioncurve+[-0.4;0.4]];
values = [0.2,0.4,0.1];
Gamma = zeros(ninclusions,npoints);
Gamma_fine = zeros(ninclusions,npoints_fine);
for i = 1:ninclusions
    index = 2*(i-1);
    Gamma(i,:) = inpolygon(meshpar.p(1,:),meshpar.p(2,:),curves(index+1,:),curves(index+2,:));
end
gamma = values(ninclusions+1)+zeros(1,npoints);
for i = 1:ninclusions
    gamma = gamma + values(i)*Gamma(i,:);
end
for i = 1:ninclusions
    index = 2*(i-1);
    Gamma_fine(i,:) = inpolygon(meshpar_fine.p(1,:),meshpar_fine.p(2,:),curves(index+1,:),curves(index+2,:));
end
gamma_fine = values(ninclusions+1)+zeros(1,npoints_fine);
for i = 1:ninclusions
    gamma_fine = gamma_fine + values(i)*Gamma_fine(i,:);
end
figure(1);
plot_from_gamma(gamma,meshpar)
figure(2);
plot_from_gamma(gamma_fine,meshpar_fine)
%%
tic;
for i = 1:1000
    c = gamma*U;
    i
end
T = toc;

%%c_fine = gamma_fine*U_fine;

%% Hvordan ser data ud for identitets forward map??
figure(1);
plot_from_gamma(c*E',meshpar)
figure(2);
plot_from_gamma(c_fine*E',meshpar)
%% Med støj
xi = randn(trunc,1);
rel_noise_level1 = 0.05;
rel_noise_level2 = 0.12;
eps = rel_noise_level1 * max(abs(c))/max(abs(xi));
eps2 = rel_noise_level2 * norm(c,2)/norm(xi,2);
% Y- norm
plot_from_gamma(c*E'+eps*xi'*E',meshpar);
eps^2