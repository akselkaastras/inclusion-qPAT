function fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad)
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
L1 = sparse(pN,pN);
L2 = sparse(pN,pN);
disp('- Building rhs')
for kk = 1:pN
    if rem(kk,round(pN/20)) == 0
        disp(['- - ', num2str(round((kk/pN)*100)),' %'])
    end
    K = sparse(pN,1);
    C = sparse(pN,1);
    nodeInd = kk;
    for jj = 1:3
        elInd = find(H(:,jj) == nodeInd);
        for ii = elInd'
            ind = H(ii,:);
            gg = p(ind,:);
            elIntGrad = triangint3gradpar(gg, jj, wfungrad);
            elInt = triangint3par(gg, 7, jj, wfun);
            C(ind) = C(ind) + elInt';
            K(ind) = K(ind) + elIntGrad';
        end
    end
    L1(:,kk) = K;
    L2(:,kk) = C;
end

fmdl.L1 = L1(fmdl.phi,:);
fmdl.L2 = L2(fmdl.phi,:);