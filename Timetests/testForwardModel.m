clear;
clc;
%%
meshrefine = 2;
meshpar = mesh_comp(meshrefine);
pN = length(meshpar.p(1,:));

wfun = @(x1,x2) 3+2*x1;
wfungrad = @(x1,x2) 2*[0*x1+1 0*x2];

center = [-0.1,0.1];
radius = 0.3;
ind_D = (meshpar.p(1,:)-center(1)).^2 + (meshpar.p(2,:)-center(2)).^2 <= radius.^2;
ind_gamma = (meshpar.p(1,:)+center(1)).^2 + (meshpar.p(2,:)+center(2)).^2 <= radius.^2;
D = ones(pN,1);
D(ind_D) = 2;
gamma = 1*ones(pN,1);
gamma(ind_gamma) = 2;
figure(1);
plot_from_gamma(D,meshpar);
figure(2);
plot_from_gamma(gamma,meshpar);
%%

fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
fmdl = fixingD(meshpar,fmdl,0.1*D');
%%
U = zeros(pN,1);
u = evalFowardModel(fmdl,meshpar,gamma);
meshpar.NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));
U(meshpar.NZ) = u;
u = U+wfun(meshpar.p(1,:)',meshpar.p(2,:)');
figure(3);
plot_from_gamma(u,meshpar);



%% Does Aint work as intended?
Nbound = meshpar.e(1,:)';
pN = size(meshpar.p,2);
NN = size(Nbound,1);

pNN = pN-NN;
gammamatrix = spdiags(gamma, 0, pN, pN);
%qmatrix = 0.001*speye(pN);

% Build mass matrix
C = fmdl.Aint*gammamatrix;
C = reshape(sum(C,2),pNN,pNN);
C = 1/2*(C'+C);
v = 1*ones(pNN,1);

0.01*v'*C*v
gamma'*fmdl.Carea*gamma

%%
gamfun = @(x,y) 3*x.^2-y+2;
gamma = gamfun(meshpar.p(1,:)',meshpar.p(2,:)');
gamma'*fmdl.Carea*gamma
fun = @(r,theta) r.*((3*abs(r.*cos(theta)).^2 - (r.*sin(theta))+2)).^2;
normsq = integral2(fun,0,1,0,2*pi)