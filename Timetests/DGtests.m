
%% Initialize mesh
meshpar = mesh_comp(0.0175);

%% Make phantom (constant = 1)
% Insert inclusion in (0,0)
inc1 = @(x,y) abs((x-0.4)+1i*(y+0.4))<0;
inc2 = @(x,y) abs((x+0.4)+1i*(y-0.4))<0;
gammaDG = 1+feval(inc1,meshpar.interp_x,meshpar.interp_y)+feval(inc2,meshpar.interp_x,meshpar.interp_y);
gamma = 1+feval(inc1,meshpar.p(1,:),meshpar.p(2,:))'+feval(inc2,meshpar.p(1,:),meshpar.p(2,:))';
gammaDG = mean(gammaDG,2);
N = length(gammaDG);
NL1 = length(gamma);
gamDGmatrix = spdiags(gammaDG, 0, N, N);
gammatrix = spdiags(gamma, 0, NL1, NL1);
p = meshpar.p';
H = meshpar.t(1:3,:)';
Ebound = meshpar.btri_ind;
Nbound = meshpar.e(1,:)';

pN = size(p,1);
NN = size(Nbound,1);
HN = size(H,1);

pNN = pN-NN;
%% Precompute for gamma = 1 projected to L1
fmdl = precomputeFEM(meshpar);
% Build mass matrix
C = fmdl.Aint*gammatrix;
C = reshape(sum(C,2),pNN,pNN);
C = 1/2*(C'+C);
%% Precompute for gamma = 1 projected to DG
fmdl_DG = precomputeFEM_DG(meshpar);
C_DG = fmdl_DG.Aint*gamDGmatrix;
C_DG = reshape(sum(C_DG,2),pNN,pNN);
C_DG = 1/2*(C_DG'+C_DG);

%% RHS
% define source from wfun
sigma = 0.5;
m1 = 0.5*[sqrt(2),sqrt(2)];
m2 = -0.5*[sqrt(2),sqrt(2)];
scale = 10;
figure(2);

% Source is Gaussian
wfunfun = @(x1,x2,m) scale*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
wfungradfun = @(x1,x2,m) 1/(sigma^2)*wfunfun(x1,x2,m).*[m(1)-x1 m(2)-x2]; 
wfun = @(x1,x2) wfunfun(x1,x2,m1) + wfunfun(x1,x2,m2);
wfungrad = @(x1,x2) wfungradfun(x1,x2,m1) + wfungradfun(x1,x2,m2);
%% RHS L1
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
Q2 = fmdl.L2*gammatrix;
Q2 = sum(Q2,2);
%% RHS DG
fmdl_DG = precomputeRHS_DG(meshpar,fmdl_DG,wfun,wfungrad);
Q2_DG = fmdl_DG.L2*gamDGmatrix;
Q2_DG = sum(Q2_DG,2);

%% whole routine

% Precomputing stiffness
fmdl = fixingD(meshpar,fmdl,gamma');

% Precomputing stiffness
fmdl_DG = fixingD(meshpar,fmdl_DG,gamma');

%%
q = 10;
alpha = 1;
tau = 0.5;
maxfreq = 10;
xq = linspace(0,2*pi,512);

priorpar = prior_init(xq,alpha,tau,q,maxfreq);

% For star-shaped inclusion parametrization
priorpar.ninclusions = 2; % number of interfaces
priorpar.v = [0.2 0.4]; % values at each interface

priorpar.background = 0.1;
priorpar.angles = linspace(0,2*pi,512)';
priorpar.mean = -2;
priorpar.dim = [priorpar.M,priorpar.ninclusions];
priorpar.type = 'star';
priorpar.center = [0.37,-0.43;-0.44,0.36];
priorpar.std = 0.2;

%%

% We use M Fourier coefficients and 2 numbers for each inclusion
xi = priorpar.std*randn(priorpar.dim);
xi = priorpar.lambdahalf.*xi;
xi_center = priorpar.center;


%% Evaluating forward model from precomputed matrices
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;
n = 10;


gamma = push_forward_star2D(xi_center,theta,meshpar.p(1,:)',meshpar.p(2,:)',priorpar);
gammaDG = push_forward_star2D_interp(xi_center,theta,meshpar.xq,meshpar.yq,meshpar.n,meshpar.HN,priorpar);
plot_from_gamma(gamma,meshpar)
figure;
pdesurf(meshpar.p,meshpar.t,gammaDG')
view(2)
colormap default
%%
tic;
for i = 1:1000
    gammaDG = push_forward_star2D_interp(xi_center,theta,meshpar.xq,meshpar.yq,meshpar.n,meshpar.HN,priorpar);
    uDG = evalFowardModel(fmdl_DG,meshpar,gammaDG);
    i
end
T1 = toc;
tic;
for i = 1:1000
    % Push-forward
    gamma = push_forward_star2D(xi_center,theta,meshpar.p(1,:)',meshpar.p(2,:)',priorpar);
    u = evalFowardModel(fmdl,meshpar,gamma);
    i
end
T2 = toc;