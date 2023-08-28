function [K, C, M] = computeLaplacian(meshpar)


%% Load mesh
p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
NN = size(Nbound,1);
HN = size(H,1);

pNN = pN-NN;
%% Build stiffness and mass matrix
Agrad = sparse(pNN*pNN,pN);
Aarea = sparse(pN*pN,pN);
Aint = sparse(pNN*pNN,pN);


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
    Aarea(:,kk) = C(:);
    C(Nbound,:) = [];
    C(:,Nbound) = [];
    K(Nbound,:) = [];
    K(:,Nbound) = [];

    Aint(:,kk) = C(:);
    Agrad(:,kk) = K(:);
end

%% Build sparsity pattern


eyematrix = speye(pN, pN);

% Build stiffness matrix
K = Agrad*eyematrix;
K = reshape(sum(K,2),pNN,pNN);
K = 1/2*(K'+K);

% Build mass matrix
C = Aint*eyematrix;
C = reshape(sum(C,2),pNN,pNN);
C = 1/2*(C'+C);

% Build area matrix
M = Aarea*eyematrix;
M = reshape(sum(M,2),pN,pN);
M = 1/2*(M'+M);


end

