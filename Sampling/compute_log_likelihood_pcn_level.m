function ll = compute_log_likelihood_pcn_level(xi, datapar, priorpar, fmdl)

% sample prior
theta = priorsample(xi,priorpar);


% push-forward
gamma = push_forward_levelset2D_smooth(theta,priorpar);
%plot_from_gamma(gamma,datapar.meshpar)

% evaluate forward model
u = evalFowardModel(fmdl,datapar.meshpar,gamma);

% add source
%U = zeros(length(datapar.meshpar.p),1);
datapar.U(datapar.meshpar.NZ) = u;
u = datapar.U + datapar.W;

% find data
data = u.*gamma;

% Find projection
c = fmdl.U_proj_coarse'*data;

% compute log-likelihood
v = c-datapar.bq;
ll = -1/(2*(datapar.epssq + datapar.sigmasq))*sum(v.^2);

