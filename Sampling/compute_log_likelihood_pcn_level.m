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

% compute log-likelihood
v = u.*gamma-datapar.bq;
%v = gamma - datapar.bq;
ll = -1/(2*datapar.epssq_approx)*normFEM(v,fmdl);
