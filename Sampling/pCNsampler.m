function [LL, N_reject, XR] = pCNsampler(datapar, samplerpar, priorpar, fmdl)

% Prior parameters
N = priorpar.N;

% Sampler parameters
N_iter = samplerpar.N_iter;
jump_size = samplerpar.jump_size;

% Initialize and compute initial log likelihood
xr = samplerpar.x0;
ll = compute_log_likelihood_pcn(xr, datapar, priorpar, fmdl);

% Initialize records
XR = zeros(N_iter,N);
N_reject = 0;          % number of rejections
LL = zeros(N_iter,1);  % likelihood
LL(1) = ll;

for k = 2:N_iter
    
    xi = sample_preprior(priorpar)';
    
    % Update with preconditioned Crank-Nicolson
    xr_new = sqrt(1-jump_size^2)*xr + jump_size*xi;
    
    ll_new = compute_log_likelihood_pcn(xr_new, datapar, priorpar, fmdl);
    
    % Compute acceptance probability
    a = min(0,ll_new-ll);

    % Accept or reject
    if log(rand) < a
        xr = xr_new;
        ll = ll_new;
    else
        N_reject = N_reject + 1; % do nothing
    end

    % Current posterior probability
    LL(k) = ll;
    XR(k,:) = xr;
    
    % Print progression every 1/100'th time increment
    if ~mod(k,ceil(N_iter/1000))
        disp([num2str(k), '/', num2str(N_iter)]);
        disp(['Rejections: ', num2str(N_reject)])
    end  
    
end
