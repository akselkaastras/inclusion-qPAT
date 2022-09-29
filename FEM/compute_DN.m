function DN = compute_DN(cond,meshpar,fmdl)

N = length(cond);
condmatrix = spdiags(cond, 0, N, N);
K = fmdl.Agrad*condmatrix;
K = reshape(sum(K,2),N,N);
K = 1/2*(K+K');
K = K + fmdl.I;

R = chol(K(fmdl.phi,fmdl.phi));
u = R\(R'\fmdl.Qperm);
u = u(fmdl.r,:);

%% Compute trace of solution

u_tr = u(meshpar.e(1,:),:);

% Nodes on the boundary
NtoD = 1/sqrt(2*pi)*meshpar.Dfii*(exp(-1i*fmdl.Nvec'*meshpar.theta)*u_tr);

% Compute DN map and project to symmetric matrix
DN = inv(NtoD);
DN = 1/2*(DN'+DN);

