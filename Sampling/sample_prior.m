function [x,w,gamma] = sample_prior(prior_par)
% Samples from the prior with parameters alpha, tau and q.


% sample pre-prior
cn = sample_preprior(prior_par);
% push forward
[x,w,gamma] = push_forward(cn,prior_par);

