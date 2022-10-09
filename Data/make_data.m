datapar = make_data(curves,values,noiselevel,meshrefine,seed); 
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

%% initialize mesh

%% make q from curves

%% make smoothened scattering mus

%% make diffusion from q and mus
% Background
mua = 0.01 * ones(size(mesh.r, 1), 1);
mus = 2 * ones(size(mesh.r, 1), 1);

% Create inclusions 
mua(vecnorm(mesh.r' - [5, 6]') < 2) = 0.1;
mus(vecnorm(mesh.r' - [2,3]') <1) = 6;

murs = mus.*(1-constant.g);
kappa = 1./constant.kind * (1./(mua + murs));
%% evaluate forward model

%% add noise, see bardsley

