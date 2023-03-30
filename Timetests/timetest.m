meshpar = mesh_comp(2);

%% Make prior for levelset with matern covariance
q = 1e5;
alpha = 3;
tau = 10;
maxfreq = 4;

priorpar = prior_init(meshpar.p(1, :), meshpar.p(2, :),alpha,tau,q,maxfreq);

% for levelset
priorpar.ninterface = 3; % number of interfaces
priorpar.c = [-inf -2 2 inf]; % contour levels (ninterface+1)
priorpar.v = [1 0 2]; % values at each interface

% sample prior
xi = randn(priorpar.M,1);
theta = priorsample(xi,priorpar);

% push-forward
gamma = push_forward_levelset2D(theta,priorpar);
%%
%% read curves
ninclusions = size(curves,1)/2;
ncurve = size(curves,2);
%% initialize mesh
meshpar = mesh_comp(meshrefine);
cmeshpar = mesh_comp(coarsemeshrefine);
npoints = length(meshpar.p);
pN = length(meshpar.p);
pNN = pN-size(meshpar.e(1,:),2);
meshpar.NZ = setdiff(1:pN,meshpar.e(1,:));
%% make gamma from curves
Gamma = zeros(ninclusions,npoints);
for i = 1:ninclusions
    index = 2*(i-1);
    Gamma(i,:) = inpolygon(meshpar.p(1,:),meshpar.p(2,:),curves(index+1,:),curves(index+2,:));
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
mus = interp2(X,Y,musimage,meshpar.p(1,:),meshpar.p(2,:));


%% make diffusion from q and mus
constantg = 0.8;   % Anisotropy factor of the Heneyey-Greenstein scattering function
constantkind = 2;  % Dimension

murs = mus.*(1-constantg);
D = 1./constantkind * (1./(gamma + murs));

%%
%% Make forward model
% define source from wfun
sigma = 0.5;
m = 0.5*[sqrt(2),sqrt(2)];
wfun = @(x1,x2) 2*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
wfungrad = @(x1,x2) 1/(sigma^2)*wfun(x1,x2).*[m(1)-x1 m(2)-x2]; 
%wfun = @(x1,x2) 0.1+0.1*x1.^2;
%wfungrad = @(x1,x2) [0.2*x1 0*x2]; 
fmdl = precomputeFEM(meshpar);
fmdl = precomputeRHS(meshpar,fmdl,wfun,wfungrad);
fmdl = fixingD(meshpar,fmdl,D);
u = evalFowardModel(fmdl,meshpar,gamma');

%%
U = zeros(pN,1);
w = wfun(meshpar.p(1,:)',meshpar.p(2,:)');


for i = 1:1000
    U(meshpar.NZ) = evalFowardModel(fmdl,meshpar,gamma');
    u = U + w;
    normFEM(u,fmdl);
    i
end