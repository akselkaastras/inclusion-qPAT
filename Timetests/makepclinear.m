meshpar = mesh_comp(0.0175);
%%
% Plane matrix
G = computePlanes(meshpar);

% These are coordinates of each point sorted by elements
points = meshpar.p(:,meshpar.t(1:3,:));

%% Prior samples
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
priorpar.type = 'star_DG';
priorpar.center = [0.37,-0.43;-0.44,0.36];
priorpar.std = 0.2;

% We use M Fourier coefficients and 2 numbers for each inclusion
xi = priorpar.std*randn(priorpar.dim);
xi = priorpar.lambdahalf.*xi;
xi_center = priorpar.center;

theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;
n = 10;

gammaDG = push_forward_star2D_interp(xi_center,theta,meshpar.xq,meshpar.yq,meshpar.n,meshpar.HN,priorpar);
gamma = push_forward_star2D(xi_center,theta,meshpar.p(1,:)',meshpar.p(2,:)',priorpar);
pdesurf(meshpar.p,meshpar.t,gammaDG')
colormap default;

%% Precompute matrices
meshpar.NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));
fmdl_DG = precomputeFEM_DG(meshpar);
fmdl = precomputeFEM(meshpar);

% RHS
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

% RHS DG
fmdl_DG = precomputeRHS_DG(meshpar,fmdl_DG,wfun,wfungrad);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);

% Precomputing stiffness
fmdl_DG = fixingD(meshpar,fmdl_DG,gamma');
fmdl = fixingD(meshpar,fmdl,gamma');
%% Get solution
uDG = evalFowardModel(fmdl_DG,meshpar,gammaDG);
U = zeros(size(meshpar.p,1),1);
U(meshpar.NZ) = uDG;
U = U + wfun(meshpar.p(1,:)',meshpar.p(2,:)');
%% Get planes of solution over each element
U_elem = U(meshpar.t(1:3,:));

%gammaDG_elem = gammaDG(meshpar.t(1:3,:));
v = G*U_elem(:);
v = reshape(v,[3 meshpar.HN]);

%% Build projection matrices
N = 15;
trunc = N*(2*N+1);
[A,B,C] = computeProjMatrix(meshpar);
E = eigenbasisFEM(meshpar,trunc);

%%
A1 = A*E;
B1 = B*E;
C1 = C*E;
tic;
for i = 1:100
    xi = priorpar.std*randn(priorpar.dim);
    xi = priorpar.lambdahalf.*xi;
    theta = priorsample(xi,priorpar);
    theta = priorpar.mean+theta;
    gammaDG = push_forward_star2D_interp(xi_center,theta,meshpar.xq,meshpar.yq,meshpar.n,meshpar.HN,priorpar);
    gamma = push_forward_star2D(xi_center,theta,meshpar.p(1,:)',meshpar.p(2,:)',priorpar);

    u = evalFowardModel(fmdl,meshpar,gamma);
    U2 = zeros(size(meshpar.p,1),1);
    U2(meshpar.NZ) = uDG;
    U2 = U2 + wfun(meshpar.p(1,:)',meshpar.p(2,:)');
    U2 = gamma.*U2;

    uDG = evalFowardModel(fmdl_DG,meshpar,gammaDG);
    U1 = zeros(size(meshpar.p,1),1);
    U1(meshpar.NZ) = uDG;
    U1 = U1 + wfun(meshpar.p(1,:)',meshpar.p(2,:)');

    U_elem = U1(meshpar.t(1:3,:));

    v = G*U_elem(:);
    v = reshape(v,[3 meshpar.HN]);


    v = gammaDG'.*v;
    c = v(1,:)*A1 + v(2,:)*B1 + v(3,:)*C1;
    c2 = U2'*M;
    norm(c-c2,2)
    pause(0.5);
    i
end
T1 = toc;
%%


%% Hvor stor forskel er der normalt?
