%% Test fourier prior against FEM prior on square domain
k = 4;
M = 128;
xi = randn(1,N);
%% FEM

meshpar = meshrect(-1,-1,1,1,4);

% Make FEM discretization of operator (-Delta)
p = meshpar.p;
H = meshpar.H;

pN = size(p,1);
HN = size(H,1);
kappa = ones(pN,1);
K = sparse(pN,pN);
A = zeros(pN,1);
for ii = 1:HN
    ind = H(ii, :);
    gg = p(ind, :);
    int3grad = triangint3grad(gg, kappa(ind));
    int3 = triangint3area(gg, kappa(ind));
    K(ind, ind) = K(ind, ind) + int3grad;
    A(ind) = A(ind) + int3;
end

% Extract N eigenvectors corresponding to the N smallest eigenvalues 
[Psi,lambda] = eigs(K,N,'smallestabs');
Psi = Psi';
lambda = diag(lambda).';

% Eigenvalues of q*(tau^2 - Delta)^{-alpha} is q*(tau^2 + lambda)^{-alpha}
lambda = q*(tau^2+lambda).^(-alpha);

% Make KL expansion to plot example
u = xi.*(lambda).^(1/2) * Psi;

figure(1);
trisurf(meshpar.t(1:3,:)', meshpar.p(1, :), meshpar.p(2, :), u,'EdgeColor','none','FaceColor','interp')
view(2)

%% Fourier

A = zeros(N,N);
A(1,1) = 1;
G = N^2*ifft2(A);
