function fmdl = precomputeFEM(meshpar)
%   Computes and saves integrals for lhs of variational form of
%         gamma nabla u . nabla v + quv    =  0
%                                       u  =  f
%   in matrices Agrad and Aint and for which U solves
%       gamma' * Agrad * U + q' * Aint * U   =  LHS(f,gamma,q)           
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

%% Build stiffness and mass matrix
Agrad = sparse(pN*pN,pN);
Aint = sparse(pN*pN,pN);

for kk = 1:pN
    if rem(kk,round(pN/20)) == 0
        disp(['- - ', num2str(round((kk/pN)*100)),' %'])
    end
    K = sparse(pN,pN);
    C = sparse(pN,pN);
    nodeInd = kk;
    for jj = 1:3
        elInd = find(H(:,jj) == nodeInd);
        for ii = elInd'
            ind = H(ii,:);
            gg = p(ind,:);
            elInt = triangint3Pre(gg, 7, jj);
            elIntGrad = triangint3gradPre(gg, 3, jj);
            C(ind, ind) = C(ind, ind) + elInt;
            K(ind, ind) = K(ind, ind) + elIntGrad;
        end
    end
    Aint(:,kk) = C(:);
    Agrad(:,kk) = K(:);
end

%% Build sparsity pattern

% Insert inclusion in (0,0)
fcn = @(x,y) abs(x+1i*y)<0.3;
gamma = 1+feval(fcn,meshpar.p(1,:),meshpar.p(2,:))';
N = length(gamma);
gammatrix = spdiags(gamma, 0, N, N);

% Build stiffness matrix from inclusion  
K = Agrad*gammatrix;
K = reshape(sum(K,2),N,N);
K = 1/2*(K'+K);
%K = K+1e-12*speye(size(K));

% Build mass matrix from inclusion (this will be changed to prior sample)
C = Aint*gammatrix;
C = reshape(sum(C,2),N,N);
C = 1/2*(C'+C);

% Optimal reordering
[phi,stats] = amd(K+C); 
r(phi) = 1:max(size(phi));
%Qperm = Q(phi,:);

%% Save matrices in struct for output
fmdl.Agrad = Agrad;
fmdl.Aint = Aint;
%fmdl.Q = Q;
%fmdl.Qperm = Qperm;
fmdl.phi = phi;
fmdl.r = r;

end

