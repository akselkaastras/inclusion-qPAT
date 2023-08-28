function datapar = make_data(curves,values,noiselevel,fine_hmax,meshpar,seed)
%% Simulates data for the qpat problem with known diffusion coefficient
% Assumptions:
%       Domain size: radius 5 mm
%       Background absorption: 0.01 mm^-1
%       Background reduced scattering 1mm^-1
%
% Input:
%   curves [2*ninclusions x ncurve] :  (x,y) values of boundary
%   values [ninclusions x 1]: absorption level of inclusion when bg = 1
%                             around 2-5 times of background
%   noiselevel: relative noiselevel of Gaussian noise added
%   meshrefine: fineness of FEM mesh, which data is simulated on
%   seed: to replicate results
%
% Output:
%   datapar [struct]: contains data qu, noisy data, relative noiselevel and
%                     absolute noiselevel, seed
%% Read user specified curves
ninclusions = size(curves,1)/2;
ncurve = size(curves,2);

%% Initialize mesh
meshpar_fine = mesh_comp(fine_hmax);

% Find mesh on boundary
npoints = length(meshpar_fine.p);
npoints_coarse = length(meshpar.p);
pN = length(meshpar_fine.p);
pNN = pN-size(meshpar_fine.e(1,:),2);
meshpar.NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));
meshpar_fine.NZ = setdiff(1:pN,meshpar_fine.e(1,:));

%% Make gamma from curves
Gamma = zeros(ninclusions,meshpar_fine.HN);
for i = 1:ninclusions
    index = 2*(i-1);
    nodes = [curves(index+1,:)',curves(index+2,:)'];
    Gamma(i,:) = Gamma(i,:) + mean(reshape(inpoly2([meshpar_fine.xq, meshpar_fine.yq],nodes),[meshpar_fine.HN meshpar_fine.n]),2)';
end
gamma = values(ninclusions+1)+zeros(1,meshpar_fine.HN);
for i = 1:ninclusions
    gamma = gamma + values(i)*Gamma(i,:);
end


Mua = zeros(ninclusions,length(meshpar_fine.p));
for i = 1:ninclusions
    index = 2*(i-1);
    nodes = [curves(index+1,:)',curves(index+2,:)'];
    Mua(i,:) = Mua(i,:) + inpoly2([meshpar_fine.p(1,:)', meshpar_fine.p(2,:)'],nodes)';
end
mua = values(ninclusions+1)+zeros(1,length(meshpar_fine.p));
for i = 1:ninclusions
    mua = mua + values(i)*Mua(i,:);
end

%% Make smoothened scattering mus with Gaussian smoothing
nimage = 300;
[X,Y] = meshgrid(linspace(-1,1,nimage));
x = X(:);
y = Y(:);
Musimage = zeros(ninclusions, nimage*nimage);

% Make image
for i = 1:ninclusions
    index = 2*(i-1);
    Musimage(i,:) = inpolygon(x,y,curves(index+1,:),curves(index+2,:));
end

% Include factor 100
musimage = 100*values(ninclusions+1)+zeros(1,nimage*nimage);
for i = 1:ninclusions
    musimage = musimage + 100*values(i)*Musimage(i,:);
end
musimage = reshape(musimage,nimage,nimage);

% Gaussian filtering
musimage = imgaussfilt(musimage,15);

% Interpolate onto mesh
mus = interp2(X,Y,musimage,meshpar_fine.p(1,:),meshpar_fine.p(2,:));


%% Make diffusion from q and mus
constantg = 0.8;   % Anisotropy factor of the Heneyey-Greenstein scattering function
constantkind = 2;  % Dimension

murs = mus.*(1-constantg);
D = 1./constantkind * (1./(mua + murs));

%% plot parameters
figure(1), clf
subplot(121)
pdesurf(meshpar_fine.p,meshpar_fine.t,gamma)
colormap default
view(2)
title('Absorption coefficient, $\gamma$','interpreter','latex','fontsize',16)
colorbar
%subplot(132)
%plot_from_gamma(full(mus),meshpar_fine)
%title('Scattering coefficient')
%colorbar
subplot(122)
plot_from_gamma(full(D),meshpar_fine)
title('Diffusion coefficient, $a$','interpreter','latex','fontsize',16)
colorbar
%% evaluate forward model
% define source from wfun
sigma = 0.5;
m1 = 0.5*[sqrt(2),sqrt(2)];
m2 = 0.5*[-sqrt(2),sqrt(2)];
scale1 = 10;
scale2 = 5;
scale3 = 2;
figure(2);

% Source is Gaussian
wfunfun = @(x1,x2,m,scale) scale*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
wfungradfun = @(x1,x2,m,scale) 1/(sigma^2)*wfunfun(x1,x2,m,scale).*[m(1)-x1 m(2)-x2]; 
wfun = @(x1,x2) wfunfun(x1,x2,m1,scale1) + wfunfun(x1,x2,-m1,scale2) + wfunfun(x1,x2,m2,scale3);  %wfunfun(x1,x2,-m2,scale2);
wfungrad = @(x1,x2) wfungradfun(x1,x2,m1,scale1) + wfungradfun(x1,x2,-m1,scale2) + wfungradfun(x1,x2,m2,scale3);  %wfungradfun(x1,x2,-m2,scale2);


% L1
% % Precomputing finite element matrices and rhs
% fmdl = precomputeFEM(meshpar_fine);
% fmdl = precomputeRHS(meshpar_fine,fmdl,wfun,wfungrad);
% 
% % Precomputing stiffness
% fmdl = fixingD(meshpar_fine,fmdl,D);
% 
% % Evaluating forward model from precomputed matrices
% u = evalFowardModel(fmdl,meshpar_fine,gamma');

% DG
% Precomputing finite element matrices and rhs
fmdl = precomputeFEM_DG(meshpar_fine);
fmdl = precomputeRHS_DG(meshpar_fine,fmdl,wfun,wfungrad);

% Precomputing stiffness
fmdl = fixingD(meshpar_fine,fmdl,D);

% Evaluating forward model from precomputed matrices
N = 13;
trunc = N*(2*N+1);
fmdl = computeProjectionMatrices_fine(fmdl,meshpar_fine,trunc);
u = evalFowardModel(fmdl,meshpar_fine,gamma');

%% Add source back again
U = zeros(pN,1);
U(meshpar_fine.NZ) = u;

u = U + wfun(meshpar_fine.p(1,:)',meshpar_fine.p(2,:)');
figure;
plot_from_gamma(full(u),meshpar_fine)

%% Data
% L1
%data = full(u).*gamma';

% Mixed L1 and DG
u = u(meshpar_fine.t(1:3,:));
u = fmdl.G*u(:);
u = reshape(u,[3 meshpar_fine.HN]);
        
        
data = gamma.*u;
dataq = data(1,:)*fmdl.U_proj1 + data(2,:)*fmdl.U_proj2 + data(3,:)*fmdl.U_proj3;
%% plot solution and data
figure(2), clf
subplot(121)
plot_from_gamma(U + wfun(meshpar_fine.p(1,:)',meshpar_fine.p(2,:)'),meshpar_fine)
title('Solution to PDE')
colorbar
subplot(122)
plot_from_gamma(dataq*fmdl.E',meshpar_fine)
title('Data')
colorbar

%% Project data to coarse mesh
% We do this by evaluating data in coarse mesh points
% This corresponds to minimizing approx. L^2(D) functional over 
% coarse hat-basis-functions.
%uq = interpolateMesh(full(u),meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
%dataq = uq.*gamma_coarse';

D_coarse = interpolateMesh(D',meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
gamma_true = gamma';

%% add noise
rng(seed);
%pN = length(meshpar.p);
% Define number of basis functions to include
%NN = ceil(1/4*(sqrt(1+8*pN)-1));
%MBessel = NN;
%Nzeros = NN;

% Compute norm squared of data and noisesample
%fmdl_coarse = precomputeFEM(meshpar);

% Make eigenfunctions of L^2(D)
%[E,Lambda] = eigenbasisLaplacianDisk2D(MBessel,Nzeros,meshpar);
%R = chol(fmdl_coarse.Carea);
%RR = R\speye(pN);

%xi_vec = randn(pN,1);

% Noise in basis consisting of scaled eigenvectors of mass matrix
%xi = RR*xi_vec;
% Truncation of noise
% This corresponds to NN = 40;
%NN = 40;
%trunc = NN*(2*NN+1);

% Make some noise and compute norm
%[xi, norm_noise, norm_noise_fem, ~] = make_noise(meshpar,fmdl_coarse,E,Lambda,trunc);



%N = 15;
%trunc = (2*N+1)*N;
%M_coarse = fmdl_coarse.Carea;
%E = eigenbasisFEM(meshpar_fine,trunc);
%U = fmdl.Carea*E;
%s = load('Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat');
%E_coarse = s.E_coarse;
%U_coarse = M_coarse*E_coarse;


%dataq = data'*U;


xi = randn(trunc,1);
%rel_noise_level = 0.01;
eps = noiselevel * norm(dataq,2) / norm(xi,2);
epssq = eps^2;


%epssq = noiselevel^2 * normdata / norm_noise;
%eps = noiselevel * max(abs(data)) / max(abs(xi));

bq = dataq' + eps*xi;
%figure;
plot_from_gamma(bq'*fmdl.E',meshpar_fine)
title('observation')

%% save in struct
datapar.bq = bq;
%datapar.b = b;
datapar.data = data;
datapar.gamma_true = gamma_true;
datapar.epssq = epssq;
datapar.noiselevel = noiselevel;
%datapar.epssqinf = epssqinf;
datapar.D = D;
datapar.D_coarse = D_coarse;
datapar.meshpar_fine = meshpar_fine;
datapar.meshpar = meshpar;
datapar.W = wfun(meshpar.p(1,:)',meshpar.p(2,:)');
datapar.wfun = wfun;
datapar.wfungrad = wfungrad;
datapar.U = zeros(length(meshpar.p),1);
datapar.N = N;