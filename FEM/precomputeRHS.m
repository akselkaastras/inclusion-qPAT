function fmdl = precomputeRHS(meshpar,fmdl)
%   Computes and saves integrals for rhs of variational form of
%       - gamma nabla u . nabla v + quv    =  0
%                                       u  =  f
%   for f = x_1 and f = x_2 in two matrices such that U (projection of u) 
%   solves
%       gamma' * Agrad * U + q' * Aq * U   =  LHS(f,gamma,q)  
%   with
%       LHS(f,gamma,q) = - gamma' * L_1 - q * L_2
%   where
%       L_1 = e_i . nabla v
%       L_2 = x_i . nabla v
%   Niko Hänninen 2019
%   Aksel Rasmussen 2021

% $$$ keyboard

%% Load mesh
p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
HN = size(H,1);

%% Build right-hand side  
disp('- Building rhs')
for ii = 1:HN
    ind = mesh.H(ii, :);
    gg = mesh.p(ind, :);
    int3grad = triangint3gradpar(gg, kappa(ind));
    C(ind, ind) = C(ind, ind) + int3;
    K(ind, ind) = K(ind, ind) + int3grad;
end
