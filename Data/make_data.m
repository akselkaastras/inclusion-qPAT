function datapar = make_data(curves,values,noiselevel,meshrefine,meshpar,seed)
%% Simulates data for the qpat problem with known diffusion coefficient
% Assumptions:
%       Domain size: radius 5 mm
%       Background absorption: 0.01mm^-1
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
%% read curves
ninclusions = size(curves,1)/2;
ncurve = size(curves,2);
%% initialize mesh
meshpar_fine = mesh_comp(meshrefine);

npoints = length(meshpar_fine.p);
pN = length(meshpar_fine.p);
pNN = pN-size(meshpar_fine.e(1,:),2);
meshpar.NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));
meshpar_fine.NZ = setdiff(1:pN,meshpar_fine.e(1,:));
%% make gamma from curves
Gamma = zeros(ninclusions,npoints);
for i = 1:ninclusions
    index = 2*(i-1);
    Gamma(i,:) = inpolygon(meshpar_fine.p(1,:),meshpar_fine.p(2,:),curves(index+1,:),curves(index+2,:));
end
gamma = 1+zeros(1,npoints);
for i = 1:ninclusions
    gamma = gamma + values(i)*Gamma(i,:);
end
% scaling
gamma = gamma*0.01;
%% make smoothened scattering mus
nimage = 300;
[X,Y] = meshgrid(linspace(-1,1,nimage));
x = X(:);
y = Y(:);
Musimage = zeros(ninclusions, nimage*nimage);
for i = 1:ninclusions
    index = 2*(i-1);
    Musimage(i,:) = inpolygon(x,y,curves(index+1,:),curves(index+2,:));
end
musimage = 1+zeros(1,nimage*nimage);
for i = 1:ninclusions
    musimage = musimage + values(i)*Musimage(i,:);
end
musimage = reshape(musimage,nimage,nimage);
musimage = imgaussfilt(musimage,15);
imagesc(musimage);
% interpolate onto mesh
mus = interp2(X,Y,musimage,meshpar_fine.p(1,:),meshpar_fine.p(2,:));


%% make diffusion from q and mus
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
figure(2);

% Source is Gaussian
wfun = @(x1,x2) 2*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
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

%% plot solution and data
figure(2), clf
subplot(121)
plot_from_gamma(full(u),meshpar_fine)
title('Solution to PDE')
colorbar
subplot(122)
plot_from_gamma(full(u),meshpar_fine)
title('Data')
colorbar
%% add noise
rng(seed);
data = full(u).*gamma';
xi = randn(pN,1);
% compute norm squared of data and noisesample
normdata = normFEM(data,fmdl);
normxi = normFEM(xi,fmdl);
epssq = noiselevel^2 * normdata / normxi;

b = data + sqrt(epssq)*xi;

plot_from_gamma(b,meshpar_fine)

%% Project data to coarse mesh
% We do this by evaluating data in coarse mesh points
% This corresponds to minimizing approx. L^2(D) functional over 
% coarse hat-basis-functions.
bq = interpolateMesh(b,meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
D_coarse = interpolateMesh(D',meshpar.p(1,:)',meshpar.p(2,:)',meshpar_fine);
figure;
plot_from_gamma(bq,meshpar);
%% save in struct
datapar.bq = bq;
datapar.b = b;
datapar.data = data;
datapar.epssq = epssq;
datapar.D = D;
datapar.D_coarse = D_coarse;
datapar.meshpar_fine = meshpar_fine;
datapar.meshpar = meshpar;
datapar.W = wfun(meshpar.p(1,:)',meshpar.p(2,:)');
datapar.wfun = wfun;
datapar.wfungrad = wfungrad;
datapar.U = zeros(length(meshpar.p),1);