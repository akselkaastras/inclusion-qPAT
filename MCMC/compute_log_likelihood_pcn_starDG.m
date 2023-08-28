function ll = compute_log_likelihood_pcn_starDG(xi, datapar, priorpar, fmdl)

xi_center = priorpar.center;

% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;

% push-forward
gamma = push_forward_star2D_interp(xi_center,theta,datapar.meshpar.xq,datapar.meshpar.yq,datapar.meshpar.n,datapar.meshpar.HN,priorpar);

% evaluate forward model
u = evalFowardModel(fmdl,datapar.meshpar,gamma);

% add source
datapar.U(datapar.meshpar.NZ) = u;
u = datapar.U + datapar.W;

% rearrange and compute for each mesh the planes equation
u = u(datapar.meshpar.t(1:3,:));
u = fmdl.G*u(:);
u = reshape(u,[3 datapar.meshpar.HN]);

% find data
data = gamma'.*u;
          
% find projection
data_proj = data(1,:)*fmdl.U_proj_coarse1 + data(2,:)*fmdl.U_proj_coarse2 + data(3,:)*fmdl.U_proj_coarse3;
       
% compute log-likelihood
%v = data_proj'+datapar.m-datapar.bq;
v = data_proj'-datapar.bq;

%v = gamma - datapar.bq;
ll = -1/(2*(datapar.epssq + datapar.sigmasq))*sum(v.^2);
