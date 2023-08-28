function fmdl = computeRHStest(meshpar,fmdl,wfun)
%   Computes and saves integrals for rhs of variational form of
%       - gamma nabla u . nabla v + quv  =  0
%                                     u  =  f
%   for f = x_1 and f = x_2 in two matrices such that U (projection of u) 
%   solves
%       gamma' * Agrad * U + q' * Aq * U   =  LHS(f,gamma,q)  
%   with
%       LHS(f,gamma,q) = - gamma' * L_1 - q * L_2
%   Aksel Rasmussen 2021

% $$$ keyboard

%% Load mesh
p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
NN = size(Nbound,1);
HN = size(H,1);

pNN = pN-NN;
%pNN = pN;
%% Build right-hand side 
L1 = sparse(pNN,pN);
L2 = sparse(pNN,pN);
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
            elInt = triangint3par(gg, 7, jj, wfun);
            C(ind) = C(ind) + elInt';
        end
    end
    % Add dirichlet boundary conditions to C and K before adding using
    % Nbound
    C(Nbound) = [];
    L2(:,kk) = C;
end
%% save in struct and use reordering from LHS
%fmdl.L1 = L1;
%fmdl.L2 = L2;
fmdl.L2 = L2(fmdl.phi,:);