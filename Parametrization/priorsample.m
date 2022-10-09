function theta = priorsample(xi,priorpar)


% Make KL expansion
theta = priorpar.B * (xi.*priorpar.lambdahalf);

