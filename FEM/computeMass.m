function M = computeMass(meshpar)


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

Aarea = sparse(pN*pN,pN);



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
            C(ind, ind) = C(ind, ind) + elInt;
        end
    end
    Aarea(:,kk) = C(:);

end

%% Build sparsity pattern


eyematrix = speye(pN, pN);

% Build area matrix
M = Aarea*eyematrix;
M = reshape(sum(M,2),pN,pN);
M = 1/2*(M'+M);


end

