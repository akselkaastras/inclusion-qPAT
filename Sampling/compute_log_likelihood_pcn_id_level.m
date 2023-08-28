function ll = compute_log_likelihood_pcn_id_level(xi, datapar, priorpar, fmdl)


% Make KL expansion
theta = priorsample(xi,priorpar);

% push-forward
gamma = push_forward_levelset2D_smooth(theta,priorpar);

% Find projection
c = fmdl.U'*gamma;

% compute log-likelihood
v = c-datapar.bq;
%plot_from_gamma(v,datapar.meshpar)



%ll = -1/(2*(datapar.epssq+5*1e-8))*normFEM(v,fmdl);
ll = -1/(2*(datapar.epssq+2*1e-7))*sum(v.^2);
%ll = 1/(2*(datapar.epssq+5*1e-8))*(datapar.bq'*fmdl.Carea*gamma-normFEM(gamma,fmdl));
