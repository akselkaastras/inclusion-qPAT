function ll = compute_log_likelihood_pcn_id(xi, datapar, priorpar, fmdl)

%xi = 0.2*xi_total(1:priorpar.M,:);
xi_center = priorpar.center;

% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;

% push-forward
gamma = push_forward_star2D(xi_center,theta,datapar.meshpar.p(1,:)',datapar.meshpar.p(2,:)',priorpar);



% compute log-likelihood
v = gamma-datapar.bq;
%plot_from_gamma(v,datapar.meshpar)

ll = -1/(2*(1e-5))*normFEM(v,fmdl);
