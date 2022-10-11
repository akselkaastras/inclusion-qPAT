function datapar = make_data(curves,values,noiselevel,meshrefine,seed); 
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
meshpar = mesh_comp(meshrefine);
npoints = length(meshpar.p);
%% make q from curves
Q = zeros(ninclusions,npoints);
for i = 1:ninclusions
    index = 2*(i-1);
    Q(i,:) = inpolygon(meshpar.p(1,:),meshpar.p(2,:),curves(index+1,:),curves(index+2,:));
end
q = 1+zeros(1,npoints);
for i = 1:ninclusions
    q = q + values(i)*Q(i,:);
end
% scaling
q = q*0.01;
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
mus = interp2(X,Y,musimage,meshpar.p(1,:),meshpar.p(2,:));


%% make diffusion from q and mus
constantg = 0.8;   % Anisotropy factor of the Heneyey-Greenstein scattering function
constantkind = 2;  % Dimension

murs = mus.*(1-constantg);
D = 1./constantkind * (1./(q + murs));

%% plot parameters
figure(1), clf
subplot(131)
plot_from_gamma(full(q),meshpar)
title('Absorption coefficient')
colorbar
subplot(132)
plot_from_gamma(full(mus),meshpar)
title('Scattering coefficient')
colorbar
subplot(133)
plot_from_gamma(full(D),meshpar)
title('Diffusion coefficient')
colorbar
%% evaluate forward model
% define source from wfun
sigma = 0.5;
m = 0.5*[sqrt(2),sqrt(2)];
figure(2);
%wfun = @(x1,x2) 2*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
%wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 
wfun = @(x1,x2) 0*x1 + 1;
wfungrad = @(x1,x2) [0*x1 0*x2];
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
fmdl = fixingD(meshpar,fmdl,D);
u = evalFowardModel(fmdl,meshpar,q');
%% plot solution and data
figure(2), clf
subplot(121)
plot_from_gamma(full(u),meshpar)
title('Solution to PDE')
colorbar
subplot(122)
plot_from_gamma(full(u).*q',meshpar)
title('Data')
colorbar
%% add noise
%sigma = noiselevel * norm(Ax) / sqrt(N);
%eta =  sigma * randn(N,1);
%b = Ax + eta;
%datapar = 0;
%% save in struct
