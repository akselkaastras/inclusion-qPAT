% Does K match K from precompute?
% Does C match C from precompute?
clc; clear;
addpath(genpath(pwd))
%%
meshpar = mesh_comp(2);
%% Define parameter v
v = ones(1,length(meshpar.p));
%% Load mesh
p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
NN = size(Nbound,1);
HN = size(H,1);

pNN = pN-NN;
%%
C1 = sparse(pN, pN);
K = sparse(pN, pN);
for ii = 1:pN
    ind = H(ii, :);
    gg = p(ind, :);
    %   Wind = W(ii);
    int3 = triangint3(gg, v(ind));
    int3grad = triangint3grad(gg, v(ind));
    C1(ind, ind) = C1(ind, ind) + int3;
    K(ind, ind) = K(ind, ind) + int3grad;
end
%%
ones(1,pN)*C1*ones(pN,1)
%%
fmdl = precomputeFEM(meshpar);
C = reshape(sum(fmdl.Aarea,2),pN,pN);
C = 1/2*(C'+C);
%% 
ones(1,pN)*C*ones(pN,1)