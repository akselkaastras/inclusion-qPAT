meshpar = mesh_comp(3);

wfun = @(x1,x2) sin(4*x1);
wfungrad = @(x1,x2) [4*cos(4*x1) 0*x2]; 

fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);

N = length(meshpar.p);
%condmatrix = spdiags(cond, 0, N, N);
condmatrix = speye(N);
qmatrix = 16*speye(N);

K = fmdl.Agrad*condmatrix;
K = reshape(sum(K,2),N,N);
%K = K + 1e-10*speye(N);
C = fmdl.Aint*qmatrix;
C = reshape(sum(C,2),N,N);

A = K+C;

Q1 = fmdl.L1*condmatrix;
Q1 = sum(Q1,2);
Q2 = fmdl.L2*qmatrix;
Q2 = sum(Q2,2);

Q = Q1 + Q2;

R = chol(A(fmdl.phi,fmdl.phi));
u = R\(R'\Q);
u = u(fmdl.r);

figure(1);
trisurf(meshpar.t(1:3,:)',meshpar.p(1,:)',meshpar.p(2,:)',full(u)-sin(4*meshpar.p(1,:)'),'facecolor','interp')
view(2)
figure(2);
trisurf(meshpar.t(1:3,:)',meshpar.p(1,:)',meshpar.p(2,:)',full(u),'facecolor','interp')
view(2)



