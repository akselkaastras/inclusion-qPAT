function theta = priorsample(xi,priorpar)


%theta = priorpar.B * (xi.*priorpar.lambdahalf);

theta = priorpar.B * xi;

