clc; clear;
addpath(genpath(pwd))
%%
meshpar = mesh_comp(2);
pN = length(meshpar.p);
meshpar.NZ = setdiff(1:pN,meshpar.e(1,:));

% define source from wfun
sigma = 0.5;
m = 0.5*[sqrt(2),sqrt(2)];
%m = [0,0];
%figure(2);
% Source is Gaussian
wfun = @(x1,x2) 2*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 
%wfun = @(x1,x2) x1*0+1;
%wfungrad = @(x1,x2) [x1*0, x2*0];
%wfungrad = @(x1,x2) [2*exp(2*x1) 0*x2]; 

fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
%%



pNN = pN-size(meshpar.e(1,:),2);

%condmatrix = spdiags(cond, 0, N, N);
condmatrix = 0*speye(pN);
qmatrix = 1*speye(pN);

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

Q = - Q1 - Q2;

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
funxy = @(x,y) 0*x+1;
fun = @(r,theta) r.*(r.*cos(theta));
normsq = integral2(fun,0,1,0,2*pi);
u = funxy(meshpar.p(1,:),meshpar.p(2,:));
gam = spdiags(ones(length(meshpar.p),1),0,length(meshpar.p),length(meshpar.p));
%C = fmdl.Aarea*gam;
%C = reshape(sum(C,2),length(meshpar.p),length(meshpar.p));
%C = 1/2*(C'+C);
u*(C*u')
%% Test rhs against integrals
pp = meshpar.p;
pp1 = meshpar.p(:,meshpar.NZ);
f = @(x1,x2) (sqrt(x1.^2+x2.^2)+1);
%f = @(x1,x2) x1*0+1;
ff = f(pp(1,:),pp(2,:));
gg = wfun(pp(1,:),pp(2,:));

area = gg*(fmdl.L1*ff');