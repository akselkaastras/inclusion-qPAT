clc; clear;
addpath(genpath(pwd))
%%
meshpar = mesh_comp(2);
pN = length(meshpar.p);
meshpar.NZ = setdiff(1:pN,meshpar.e(1,:));
wfun = @(x1,x2) exp(2*x1);
wfungrad = @(x1,x2) [2*exp(2*x1) 0*x2]; 

fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
%%

pNN = pN-size(meshpar.e(1,:),2);

%condmatrix = spdiags(cond, 0, N, N);
condmatrix = speye(pN);
qmatrix = 4*speye(pN);

K = fmdl.Agrad*condmatrix;
K = reshape(sum(K,2),pNN,pNN);
K = 1/2*(K'+K);
%K = K + 1e-5*speye(N);

C = fmdl.Aint*qmatrix;
C = reshape(sum(C,2),pNN,pNN);
C = 1/2*(C'+C);
%%
A = K+C;

Q1 = fmdl.L1*condmatrix;
Q1 = sum(Q1,2);
Q2 = fmdl.L2*qmatrix;
Q2 = sum(Q2,2);

Q = Q1 + Q2;

R = chol(A(fmdl.phi,fmdl.phi));
u = R\(R'\Q);
u = u(fmdl.r);

%% Add zero nodes again
U1 = zeros(pN,1);
U1(meshpar.NZ) = u;
U2 = U1 + wfun(meshpar.p(1,:)',meshpar.p(2,:)');
%%
figure(1);
trisurf(meshpar.t(1:3,:)',meshpar.p(1,:)',meshpar.p(2,:)',full(U1),'EdgeColor','none','facecolor','interp')
view(2)
%% Integral of solution squared
fun = @(r,theta) r.*sin(4*r.*cos(theta));
normsq = integral2(fun,0,1,0,2*pi);
% u' C u
%%
figure(3)
T = linspace(-1,1,100);
[X,Y] = meshgrid(T);
vq = interpolateMesh(u,X(:),Y(:),meshpar);
vq = reshape(vq,size(X));
surf(X,Y,vq,'EdgeColor','none')
view(2)
%% test prior
q = 10;
alpha = 1;
tau = 15;
maxfreq = 4;
usum = zeros(1,length(meshpar));

priorpar = prior_init(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,4);
trisurf(meshpar.t(1:3,:)', meshpar.p(1, :), meshpar.p(2, :), priorpar.u,'EdgeColor','none','FaceColor','interp')
view(2)