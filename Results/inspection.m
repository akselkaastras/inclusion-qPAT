%% mesh
clc;
clear;
addpath(genpath(pwd))
%% Make mesh
hmax = 0.0175;
meshpar = mesh_comp(hmax);

%% What does the data look like?

load('/work3/akara/qPAT-level/Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat')
N = 13;
trunc = (2*N+1)*N;
plot_from_gamma(datapar.bq'*E_coarse(:,1:trunc)',meshpar)

%% Prior starDG


%%

xr = results.XR(end,:);
xr = reshape(xr,results.priorpar.dim);
figure;
plot_from_coef_star(xr,results.priorpar)
xr = results.XR(end-15000,:);
xr = reshape(xr,results.priorpar.dim);
figure;
plot_from_coef_star(xr,results.priorpar)

%%
results = burnthin(results,5e5,1);
%% Plot chains
plot_chains_star(results,5)
%% Autocorrelation
plot_autocorr_star(results,5);
%% Plot of mean
plot_from_coef_star_mean(results,results.priorpar)


%%


xr = results.XR(end,:);
figure;
plot_from_coef_level_smooth(xr',results.priorpar)
%figure;
%plot_from_coef_2D(xr',results.priorpar)
%colorbar;
%set(gca, 'clim', [-2 2]);
xr = results.XR(1000000,:);
figure;
plot_from_coef_level_smooth(xr',results.priorpar)

%%
