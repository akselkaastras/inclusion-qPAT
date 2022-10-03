function ll = compute_log_likelihood_pcn(cn, datapar, priorpar, fmdl)

[x,w,gamma] = push_forward(cn,priorpar);


u = forwardFEM(gamma',fmdl);

ll = -1/(2*datapar.delta^2)*sum(abs(u-datapar.u_eps).^2);

