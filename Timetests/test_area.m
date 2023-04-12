%% Is Carea correct?
meshpar = mesh_comp(3);
fmdl = precomputeFEM(meshpar);
%% Generate fine uniform grid on [-1,1]^2, define function on it
%rng(0);
h = 0.005;
[R,Theta] = meshgrid(0:h:1,0:h:2*pi);
N1 = size(R,1);
N2 = size(R,2);
f = @(r,theta) 2*rand*r.^2+2*cos(10*r);
r = reshape(R,N1*N2,1);
theta = reshape(Theta,N1*N2,1);
z = f(r,theta);


%% Interpolate onto mesh and area
Xq = meshpar.p(1,:)';
Yq = meshpar.p(2,:)';
[thetaq,rq] = cart2pol(Xq,Yq);
Vq = f(rq,thetaq);
area1 = normFEM(Vq,fmdl);

%% Area from fine mesh
area2 = h^2*sum(z.^2.*r);


%% Plots

plot_from_gamma(Vq,meshpar)

%%
Area_error = area1-area2