clc; clear;
addpath(genpath(pwd))

%% prior samples

meshpar = mesh_comp(2);
p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
HN = size(H,1);
kappa = ones(pN,1);
K = sparse(pN,pN);
C = sparse(pN+1,pN+1);
A = zeros(pN,1);
for ii = 1:HN
    ind = H(ii, :);
    gg = p(ind, :);
    int3grad = triangint3grad(gg, kappa(ind));
    int3 = triangint3area(gg, kappa(ind));
    K(ind, ind) = K(ind, ind) + int3grad;
    A(ind) = A(ind) + int3;
end
C(1:pN,1:pN) = K;
C(pN+1,1:pN) = A';
C(1,pN+1) = 1;
%% Does a solution to the Neumann problem really integrate to 1?
b = abs(randn(pN+1,1));
u = C\b;
%%
[V,d] = eigs(K,100,'smallestabs');

trisurf(meshpar.t(1:3,:)', meshpar.p(1, :), meshpar.p(2, :), V(:,4)','EdgeColor','none','FaceColor','interp'), view(2)
