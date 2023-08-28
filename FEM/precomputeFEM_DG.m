function fmdl = precomputeFEM_DG(meshpar)
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
NN = size(Nbound,1);
HN = size(H,1);

pNN = pN-NN;
%% Build stiffness and mass matrix
Agrad = sparse(pNN*pNN,pN);
Aint = sparse(pNN*pNN,HN);
Aarea = sparse(pN*pN,pN);

disp(' - - Building mass matrix');
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
    K(Nbound,:) = [];
    K(:,Nbound) = [];
    Agrad(:,kk) = K(:);
end

disp('Building FEM matrices for DG basis functions')
for kk = 1:HN
    if rem(kk,round(pN/20)) == 0
        disp(['- - ', num2str(round((kk/HN)*100)),' %'])
    end
    C = sparse(pN,pN);
    elInd = kk;

    ind = H(elInd,:);
    gg = p(ind,:);
    elInt = triangint3(gg,[1,1,1]);

    C(ind, ind) = C(ind, ind) + elInt;

    % Add dirichlet boundary conditions to C and K before adding using
    % Nbound
    C(Nbound,:) = [];
    C(:,Nbound) = [];
    Aint(:,kk) = C(:);

end

%% Build sparsity pattern

% Insert inclusion in (0,0)
inc1 = @(x,y) abs((x-0.4)+1i*(y+0.4))<0.25;
inc2 = @(x,y) abs((x+0.4)+1i*(y-0.4))<0.25;
gammaDG = 1+feval(inc1,meshpar.interp_x,meshpar.interp_y)+feval(inc2,meshpar.interp_x,meshpar.interp_y);
gamma = 1+feval(inc1,meshpar.p(1,:),meshpar.p(2,:))'+feval(inc2,meshpar.p(1,:),meshpar.p(2,:))';
gammaDG = mean(gammaDG,2);
N = length(gammaDG);
NL1 = length(gamma);
gamDGmatrix = spdiags(gammaDG, 0, N, N);
gammatrix = spdiags(gamma, 0, NL1, NL1);

% Build stiffness matrix from inclusion  
K = Agrad*gammatrix;
K = reshape(sum(K,2),pNN,pNN);
K = 1/2*(K'+K);
%K = K+1e-12*speye(size(K));

% Build mass matrix from inclusion (this will be changed to prior sample)
C = Aint*gamDGmatrix;
C = reshape(sum(C,2),pNN,pNN);
C = 1/2*(C'+C);

% Build "area" matrix. Simply mass matrix corresponding to gamma = 1
Carea = reshape(sum(Aarea,2),pN,pN);
Carea = 1/2*(Carea'+Carea);

%% Optimal reordering
[phi,stats] = amd(K+C); 
r(phi) = 1:max(size(phi));
%Agrad = Agrad(phi,:);
%Aint = Aint(phi,:);

%% Save matrices in struct for output
fmdl.Agrad = Agrad;
fmdl.Aint = Aint;
fmdl.phi = phi;
fmdl.r = r;
fmdl.Carea = Carea;

end

