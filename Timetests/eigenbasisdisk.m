%% Test noise model with eigenbasis in L^2(D)
clc;
clear;
% Mesh to project basis functions onto
hmax = 0.03;
meshpar = mesh_comp(hmax);
fmdl = precomputeFEM(meshpar);

% Define number of basis functions to include
MBessel = 20;
Nzeros = 20;

[E,Lambda] = eigenbasisLaplacianDisk2D(MBessel,Nzeros,meshpar);
[xi, norm_noise, norm_noise_fem] = make_noise(meshpar,fmdl,E,Lambda);
