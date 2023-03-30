function priorpar = prior_init(xq,alpha,tau,q,maxfreq)

% Finds eigenbasis (lambda_n, \psi_n) of q*(tau^2-Delta)^(-alpha) on
% [0,2*pi] from fourier basis
% mesh_par.
% Input:    
%       [xq]    : points in which Fourier basis is evaluated
%       alpha   : regularity of prior samples (decay of eigenvalues) > 0
%       tau     : inverse length scale > 0
%       q       : amplitude
%       maxfreq : maximum freq. to include in basis

% Matern covariance parameters
priorpar.alpha = alpha;
priorpar.tau = tau;
priorpar.q = q;

evalpar = makeEvalMatrixFourier(xq,maxfreq);
priorpar.M = evalpar.M;
priorpar.maxfreq = evalpar.maxfreq;

% Eigenvalues of q*(tau^2 - Delta)^{-alpha} is q*(tau^2 + lambda)^{-alpha}
lambda = q*(tau^2+evalpar.Lambda).^(-alpha);


% Make KL expansion to plot example
xi = randn(priorpar.M,1);
u = evalpar.B * (xi.*(lambda).^(1/2));


% Save eigenvalues and eigenvectors
priorpar.B = evalpar.B;
priorpar.lambda = lambda;
priorpar.lambdahalf = lambda.^(1/2);
priorpar.u = u;
priorpar.res = evalpar.res;
