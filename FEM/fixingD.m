function fmdl = fixingD(meshpar,fmdl,D)

N = length(meshpar.p);
Dmatrix = spdiags(D', 0, N, N);


K = fmdl.Agrad*Dmatrix;
K = reshape(sum(K,2),N,N);
K = 1/2*(K'+K);

Q1 = fmdl.L1*Dmatrix;
Q1 = sum(Q1,2);

fmdl.K = K(fmdl.phi,fmdl.phi);
fmdl.Q1 = Q1;