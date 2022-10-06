function theta = priorsample(xi,priorpar)


% Make KL expansion
theta = xi.*(priorpar.lambda).^(1/2) * priorpar.Psi;

