function u = evalFowardModel(fmdl,meshpar,q)

pN = length(meshpar.p);
pNN = pN-size(meshpar.e(1,:),2);
qmatrix = spdiags(q, 0, pN, pN);

% Build mass matrix
C = fmdl.Aint*qmatrix;
C = reshape(sum(C,2),pNN,pNN);
C = 1/2*(C'+C);

% System matrix A and reordering of C
A = fmdl.K+C(fmdl.phi,fmdl.phi);

% Build rhs
Q2 = fmdl.L2*qmatrix;
Q2 = sum(Q2,2);

Q = - fmdl.Q1 - Q2;

% Solve
R = chol(A);
u = R\(R'\Q);
u = u(fmdl.r);
