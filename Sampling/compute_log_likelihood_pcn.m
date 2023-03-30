function ll = compute_log_likelihood_pcn_level(xi, datapar, priorpar, fmdl)

% sample prior
theta = priorsample(xi,priorpar);

% push-forward
gamma = push_forward_levelset2D_smooth(theta,priorpar);

u = forwardFEM(gamma',fmdl);

ll = -1/(2*datapar.delta^2)*sum(abs(u-datapar.u_eps).^2);

