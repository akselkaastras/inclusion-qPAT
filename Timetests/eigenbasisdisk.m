%% Test noise model with eigenbasis in L^2(D)
clc;
clear;
% Mesh to project basis functions onto
hmax = 0.05;
meshpar = mesh_comp(hmax);

% Define number of basis functions to include
MBessel = 50;
Nzeros = 50;

[E,Lambda] = eigenbasisLaplacianDisk2D(MBessel,Nzeros,meshpar);

%%
for j = 1:1000
    plot_from_gamma(E(j,:),meshpar)
    pause(0.1);
end
%% Noise
e = randn(1,size(E,1));
v = e * E;
plot_from_gamma(v,meshpar)

%% Orthogonal?
fmdl = precomputeFEM(meshpar);
%%
for j = 1:100
    E(j,:)*fmdl.Carea*E(j,:)'
    pause(0.2);
end
%% Orthonormal?
for j = 1:size(E,1)
    if abs(E(j,:)*fmdl.Carea*E(380,:)') > 0.01 && j~=380
        j
        error('hey');
    end
end
E(2972,:)*fmdl.Carea*E(380,:)'