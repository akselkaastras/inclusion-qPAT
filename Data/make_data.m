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
Gamma = zeros(ninclusions,npoints);
for i = 1:ninclusions
    index = 2*(i-1);
    Gamma(i,:) = inpolygon(meshpar_fine.p(1,:),meshpar_fine.p(2,:),curves(index+1,:),curves(index+2,:));
end
gamma = values(ninclusions+1)+zeros(1,npoints);
for i = 1:ninclusions
    gamma = gamma + values(i)*Gamma(i,:);
end

% Gamma on coarse grid
Gamma_coarse = zeros(ninclusions,npoints_coarse);
for i = 1:ninclusions
    index = 2*(i-1);
    Gamma_coarse(i,:) = inpolygon(meshpar.p(1,:),meshpar.p(2,:),curves(index+1,:),curves(index+2,:));
end
gamma_coarse = values(ninclusions+1)+zeros(1,npoints_coarse);
for i = 1:ninclusions
    gamma_coarse  = gamma_coarse + values(i)*Gamma_coarse(i,:);
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
D = 1./constantkind * (1./(gamma + murs));

%% plot parameters
figure(1), clf
subplot(131)
plot_from_gamma(full(gamma),meshpar_fine)
title('Absorption coefficient')
colorbar
subplot(132)
plot_from_gamma(full(mus),meshpar_fine)
title('Scattering coefficient')
colorbar
subplot(133)
plot_from_gamma(full(D),meshpar_fine)
title('Diffusion coefficient')
colorbar
%% evaluate forward model
% define source from wfun
sigma = 0.5;
m = 0.5*[sqrt(2),sqrt(2)];
scale = 10;
figure(2);

% Source is Gaussian
wfun = @(x1,x2) scale*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 

% Precomputing finite element matrices and rhs
fmdl = precomputeFEM(meshpar_fine);
fmdl = precomputeRHS(meshpar_fine,fmdl,wfun,wfungrad);

% Precomputing stiffness
fmdl = fixingD(meshpar_fine,fmdl,D);

% Evaluating forward model from precomputed matrices
u = evalFowardModel(fmdl,meshpar_fine,gamma');

%% Add source back again
U = zeros(pN,1);
U(meshpar_fine.NZ) = u;

u = U + wfun(meshpar_fine.p(1,:)',meshpar_fine.p(2,:)');
figure;
plot_from_gamma(full(u),meshpar_fine)

%% Data
data = full(u).*gamma';

%% plot solution and data
figure(2), clf
subplot(121)
plot_from_gamma(full(u),meshpar_fine)
title('Solution to PDE')
colorbar
subplot(122)
plot_from_gamma(data,meshpar_fine)
title('Data')
colorbar

%% Project data to coarse mesh
% We do this by evaluating data in coarse mesh points
% This corresponds to minimizing approx. L^2(D) functional over 
% coarse hat-basis-functions.
uq = interpolateMesh(full(u),meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
dataq = uq.*gamma_coarse';

D_coarse = interpolateMesh(D',meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
gamma_true = gamma_coarse';

%% add noise
rng(seed);
pN = length(meshpar.p);
% Define number of basis functions to include
NN = ceil(1/4*(sqrt(1+8*pN)-1));
MBessel = NN;
Nzeros = NN;

% Compute norm squared of data and noisesample
fmdl_coarse = precomputeFEM(meshpar);

% Make eigenfunctions of L^2(D)
[E,Lambda] = eigenbasisLaplacianDisk2D(MBessel,Nzeros,meshpar);

% Make some noise and compute norm
[xi, norm_noise, norm_noise_fem, ~] = make_noise(meshpar,fmdl_coarse,E,Lambda);


normdata = normFEM(data,fmdl);

epssq = noiselevel^2 * normdata / norm_noise;

bq = dataq + sqrt(epssq)*xi';
figure;
plot_from_gamma(bq,meshpar)
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