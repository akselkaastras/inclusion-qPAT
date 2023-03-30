function fmdl = fixingD(meshpar,fmdl,D)

pN = length(meshpar.p);
pNN = pN-size(meshpar.e(1,:),2);
Dmatrix = spdiags(D', 0, pN, pN);


K = fmdl.Agrad*Dmatrix;
K = reshape(sum(K,2),pNN,pNN);
K = 1/2*(K'+K);

Q1 = fmdl.L1*Dmatrix;
Q1 = sum(Q1,2);

fmdl.K = K(fmdl.phi,fmdl.phi);
fmdl.Q1 = Q1;