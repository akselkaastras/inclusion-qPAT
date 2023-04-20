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
%A = fmdl.K+C;

% Build rhs
Q2 = fmdl.L2*qmatrix;
Q2 = sum(Q2,2);
%fmdl.Q1 = fmdl.Q1(fmdl.phi);
%Q2 = Q2(fmdl.phi);

Q =  - fmdl.Q1 - Q2;

% Solve
R = chol(A);
u = R\(R'\Q);
u = u(fmdl.r);
