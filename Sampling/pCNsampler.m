function [LL, N_reject, XR] = pCNsampler(datapar, samplerpar, priorpar, fmdl)

% Prior parameters
dim = priorpar.dim;

% Sampler parameters
N_iter = samplerpar.N_iter;
jump_size = samplerpar.jump_size;

% Initialize and compute initial log likelihood
xr = samplerpar.x0;
if strcmp(priorpar.type,'level')
    ll = compute_log_likelihood_pcn_level(xr, datapar, priorpar, fmdl);
elseif strcmp(priorpar.type,'star')
    ll = compute_log_likelihood_pcn_star(xr, datapar, priorpar, fmdl); 
end

% Initialize records
NN = prod(dim);
XR = zeros(N_iter,NN);
XR(1,:) = reshape(xr,NN,1);
N_reject = 0;          % number of rejections
LL = zeros(N_iter,1);  % likelihood
LL(1) = ll;
acc = zeros(N_iter,1); % acceptance history
acc(1) = 1;        % initial guess is accepted

% Adaptation initialize
Na = floor(0.05*N_iter);
hat_acc = zeros(floor(N_iter/Na),1);
star_acc = 0.3; % target acceptance rate
i = 1;
idx = 1;


for k = 2:N_iter
    
    % Sample preprior
    xi = priorpar.std*randn(dim);

    % Update with preconditioned Crank-Nicolson
    xr_new = sqrt(1-jump_size^2)*xr + jump_size*xi;
    
    % Compute log-likelihood
    if strcmp(priorpar.type,'level')
        ll_new = compute_log_likelihood_pcn_level(xr_new, datapar, priorpar, fmdl);
    elseif strcmp(priorpar.type,'star')
        ll_new = compute_log_likelihood_pcn_star(xr_new, datapar, priorpar, fmdl); 
    end
    
    % Compute acceptance probability
    a = min(0,ll_new-ll);

    % Accept or reject
    if log(rand) < a
        xr = xr_new;
        ll = ll_new;
        acc(k) = 1;
    else
        N_reject = N_reject + 1; % do nothing
        acc(k) = 0;
    end

    % Adapt stepsize
    if mod(k-1,Na) == 0
        % Evaluate acceptance rate
        hat_acc(i) = mean(acc(idx:idx+Na));

        % Compute new jump parameter
        zeta = 1/sqrt(i);
        jump_size = exp(log(jump_size) + zeta*(hat_acc(i)-star_acc));
        jump_size = min(jump_size,1);
        disp(['-- jump_size updated to ', num2str(jump_size)])

        % Update counters
        i = i+1;
        idx = idx+Na;
    end

    % Current likelihood
    LL(k) = ll;
    XR(k,:) = reshape(xr,NN,1);
    
    % Print progression every 1/100'th time increment
    if ~mod(k,ceil(N_iter/1000))
        disp(['- - ', num2str(k), '/', num2str(N_iter)]);
        disp(['- - ', 'Rejections: ', num2str(N_reject)])
    end  
    
end
