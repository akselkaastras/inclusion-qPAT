meshpar = mesh_comp(0.0175);

%%
% Make interpolation points of elements
% Based on EIDORS elem_select
n_interp = 3;
el_dim = 2;
% Get element nodes, and reshape
% need to be of size n_dims_1 x (n_elems*n_dims) for reshape
t = meshpar.t(1:3,:)';
p = meshpar.p';
el_nodes= p(t',:);
el_nodes= reshape(el_nodes, el_dim+1, []);

% Get interpolation matrix
interp= zeros(0,el_dim+1);
for i=0:n_interp
    for j=0:n_interp-i
        interp= [interp;i,j,n_interp-i-j];
    end
end
interp= (interp + 1/(el_dim+1) )/(n_interp+1);

l_interp = size(interp,1);
mdl_pts = interp*el_nodes;
mdl_pts = reshape(mdl_pts, l_interp, length(meshpar.t(1:3,:)), el_dim);
mdl_pts = permute(mdl_pts, [2,3,1]);

interp_x = squeeze(mdl_pts(:,1,:));
interp_y = squeeze(mdl_pts(:,2,:));

%% inpolygon approach
n = 10;
HN = size(meshpar.t,2);

xq = reshape(interp_x,[n*HN 1]);
yq = reshape(interp_y,[n*HN 1]);

% We use M Fourier coefficients and 2 numbers for each inclusion
%xi = priorpar.std*randn(priorpar.dim);
%xi = priorpar.lambdahalf.*xi;
%xi_center = priorpar.center;


% Make KL expansion


tic;
for i = 1:100
    theta = priorsample(xi,priorpar);
    theta = priorpar.mean+theta;
    gamma = push_forward_star2D_interp(xi_center,theta,xq,yq,n,HN,priorpar);
    i
end
T = toc;
%% plot
figure;
pdesurf(meshpar.p,meshpar.t,gamma')
%figure;
%plot_from_coef_star(xi,priorpar)
%% feval approach

tic;
for i = 1:100
    theta = priorsample(xi,priorpar);
    theta = priorpar.mean+theta;
    gamma = xqnew.^2 + yqnew.^2 <= exp(theta(:,1));
    gamma = mean(reshape(gamma,HN,10),2);
    i
end
T = toc;
