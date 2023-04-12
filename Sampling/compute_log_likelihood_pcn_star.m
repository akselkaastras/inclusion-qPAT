function ll = compute_log_likelihood_pcn_star(xi, datapar, priorpar, fmdl)

%xi = 0.2*xi_total(1:priorpar.M,:);
xi_center = priorpar.center;

% Make KL expansion
theta = priorsample(xi,priorpar);
theta = priorpar.mean+theta;

% push-forward
gamma = push_forward_star2D(xi_center,theta,datapar.meshpar.p(1,:)',datapar.meshpar.p(2,:)',priorpar);
%plot_from_gamma(gamma,datapar.meshpar)

% evaluate forward model
u = evalFowardModel(fmdl,datapar.meshpar,gamma);

% add source
%U = zeros(length(datapar.meshpar.p),1);
datapar.U(datapar.meshpar.NZ) = u;
u = datapar.U + datapar.W;

% compute log-likelihood
v = u.*gamma-datapar.bq;
figure(1)
plot_from_gamma(u.*gamma,datapar.meshpar)
figure(2)
plot_from_gamma(datapar.bq,datapar.meshpar)
%v = gamma - datapar.bq;
ll = -1/(2*datapar.epssq_approx)*normFEM(v,fmdl);
