function priorpar = prior_init_from_mesh(meshpar,alpha,tau,q,N)

% Finds eigenbasis (lambda_n, \psi_n) of (-Delta) on mesh provided by
% mesh_par.
% Input:    
%       mesh_par: struct of mesh parameters, points, connectivity etc.
%       alpha   : regularity of prior samples (decay of eigenvalues) > 0
%       tau     : inverse length scale > 0
%       q       : amplitude
%       N       : number of basisfunctions to include

% Matern covariance parameters
priorpar.alpha = alpha;
priorpar.tau = tau;
priorpar.q = q;

% Make FEM discretization of operator (-Delta)
p = meshpar.p';
H = meshpar.t(1:3,:)';

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
xi = randn(1,N);
u = xi.*(lambda).^(1/2) * Psi;

figure(1);
trisurf(meshpar.t(1:3,:)', meshpar.p(1, :), meshpar.p(2, :), u,'EdgeColor','none','FaceColor','interp')
view(2)

% Make also a plot of the restriction to a line y = 0;
F = scatteredInterpolant(meshpar.p(1, :)', meshpar.p(2, :)', u');
t = linspace(-1,1,1000);
theta = pi/4;
a = [cos(theta),sin(theta)];
figure(2);
plot(t,F(a(1)*t,a(2)*t));
figure(1);
hold on
plot3(a(1)*t,a(2)*t,100+0.*t,'r-','linewidth',2)
hold off

% Save eigenvalues and eigenvectors
priorpar.Psi = Psi;
priorpar.lambda = lambda;
priorpar.u = u;
