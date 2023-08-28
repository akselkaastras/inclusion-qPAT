function [LL, N_reject, XR] = ApCNsampler(datapar, samplerpar, priorpar, fmdl, x0)

% Prior parameters
dim = priorpar.dim;

% Update lambda according to prior dimension
lambda = [priorpar.lambda priorpar.lambda];

% Sampler parameters
N_iter = samplerpar.N_iter;
jump_size = samplerpar.jump_size;

% Preburn-in
N_preburnin = 1;
% Burn-in for adaptation
N_burnin = 1e3;

% Initialize and compute initial log likelihood
xr = x0;
if strcmp(priorpar.type,'level')
    ll = compute_log_likelihood_pcn_level(xr, datapar, priorpar, fmdl);
elseif strcmp(priorpar.type,'star')
    ll = compute_log_likelihood_pcn_star(xr, datapar, priorpar, fmdl);
elseif strcmp(priorpar.type,'starDG')
    ll = compute_log_likelihood_pcn_starDG(xr, datapar, priorpar, fmdl); 
elseif strcmp(priorpar.type,'id')
    ll = compute_log_likelihood_pcn_id(xr, datapar, priorpar, fmdl);
end

% Initialize records
NN = prod(dim);
N_reject = 0;          % number of rejections

acc = zeros(N_preburnin,1); % acceptance history
acc(1) = 1;        % initial guess is accepted


% Pre burn-in
for k = 1:N_preburnin
    % Sample preprior
    xi = priorpar.std*randn(dim);
    xi = priorpar.lambdahalf.*xi;

    % Update with preconditioned Crank-Nicolson
    %xr_new = sqrt(1-jump_size^2)*xr + jump_size*xi;
    xr_new = sqrt(1-jump_size^2)*xr + jump_size*xi;

    % Compute log-likelihood
    if strcmp(priorpar.type,'level')
        ll_new = compute_log_likelihood_pcn_level(xr_new, datapar, priorpar, fmdl);
    elseif strcmp(priorpar.type,'star')
        ll_new = compute_log_likelihood_pcn_star(xr_new, datapar, priorpar, fmdl);
    elseif strcmp(priorpar.type,'starDG')
        ll_new = compute_log_likelihood_pcn_starDG(xr_new, datapar, priorpar, fmdl); 
    elseif strcmp(priorpar.type,'id')
        ll_new = compute_log_likelihood_pcn_id(xr_new, datapar, priorpar, fmdl);
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
    
    % Print progression every 1/10'th time increment
    b = ceil(N_preburnin/10);
    if ~mod(k,b)
        disp(['- - Pre-burn-in phase - - ', num2str(k), '/', num2str(N_preburnin)]);
        disp(['- - ', 'Rejections: ', num2str(N_reject)]);
        n = k/b;
        acc_since_last_time = sum(acc(((n-1)*b+1):(n*b)))/b*100;
        disp(['- - ', 'Accepted: ', sprintf('%.1f%%',acc_since_last_time)]);
    end  
end
disp('Pre burn-in complete');
disp('Starting burn-in')

%% Burn-in phase - Now we record the samples
XR_burnin = zeros(N_burnin,NN);
XR_burnin(1,:) = reshape(xr,NN,1);

LL_burnin = zeros(N_burnin,1);  % likelihood
LL_burnin(1) = ll;
acc = zeros(N_burnin,1); % acceptance history
acc(1) = 1;        % initial guess is accepted
N_reject = 0;
for k = 2:N_burnin
    
    % Sample preprior
    xi = priorpar.std*randn(dim);
    xi = priorpar.lambdahalf.*xi;

    % Update with preconditioned Crank-Nicolson
    xr_new = sqrt(1-jump_size^2)*xr + jump_size*xi;

    % Compute log-likelihood
    if strcmp(priorpar.type,'level')
        ll_new = compute_log_likelihood_pcn_level(xr_new, datapar, priorpar, fmdl);
    elseif strcmp(priorpar.type,'star')
        ll_new = compute_log_likelihood_pcn_star(xr_new, datapar, priorpar, fmdl);
    elseif strcmp(priorpar.type,'starDG')
        ll_new = compute_log_likelihood_pcn_starDG(xr_new, datapar, priorpar, fmdl); 
    elseif strcmp(priorpar.type,'id')
        ll_new = compute_log_likelihood_pcn_id(xr_new, datapar, priorpar, fmdl);
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

    % Current likelihood
    LL_burnin(k) = ll;
    XR_burnin(k,:) = reshape(xr,NN,1);
    
    % Print progression every 1/1000'th time increment
    b = ceil(N_burnin/10);
    if ~mod(k,b)
        disp(['- - Burn-in phase - - ', num2str(k), '/', num2str(N_burnin)]);
        disp(['- - ', 'Rejections: ', num2str(N_reject)]);
        n = k/b;
        acc_since_last_time = sum(acc(((n-1)*b+1):(n*b)))/b*100;
        disp(['- - ', 'Accepted: ', sprintf('%.1f%%',acc_since_last_time)]);
    end  
    
end

%% Compute variance of samples based on burn-in
m = mean(XR_burnin,1);
s = sum(XR_burnin.^2,1);
a = var(XR_burnin,1,1);
a = a + 1e-8;
m = reshape(m,dim);
s = reshape(s,dim);
a = reshape(a,dim);
a = min(a,lambda);
a = a./max(a);
%% Continiue sampling with aPCN
disp('Sampling with adapted covariance')
%% Sampling adaptively
XR = zeros(N_iter,NN);
XR(1,:) = reshape(xr,NN,1);

LL = zeros(N_iter,1);  % likelihood
LL(1) = ll;
acc = zeros(N_iter,1); % acceptance history
acc(1) = 1;        % initial guess is accepted
N_reject = 0;

% Adaptation of step-size initialize
Na = floor(0.05*N_iter);
hat_acc = zeros(floor(N_iter/Na),1);
star_acc = 0.3; % target acceptance rate
i = 1;
idx = 1;

for k = 2:N_iter
    
    % Sample preprior
    xi = priorpar.std*randn(dim);
    xi = sqrt(a).*xi;

    % Update with preconditioned Crank-Nicolson
    xr_new = sqrt(1-jump_size^2*a./lambda).*xr + jump_size*xi;

    % Compute log-likelihood
    if strcmp(priorpar.type,'level')
        ll_new = compute_log_likelihood_pcn_level(xr_new, datapar, priorpar, fmdl);
    elseif strcmp(priorpar.type,'star')
        ll_new = compute_log_likelihood_pcn_star(xr_new, datapar, priorpar, fmdl);
    elseif strcmp(priorpar.type,'starDG')
        ll_new = compute_log_likelihood_pcn_starDG(xr_new, datapar, priorpar, fmdl); 
    elseif strcmp(priorpar.type,'id')
        ll_new = compute_log_likelihood_pcn_id(xr_new, datapar, priorpar, fmdl);
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
%     if mod(k-1,Na) == 0
%         % Evaluate acceptance rate
%         hat_acc(i) = mean(acc(idx:idx+Na));
% 
%         % Compute new jump parameter
%         zeta = 1/sqrt(i);
%         jump_size = exp(log(jump_size) + zeta*(hat_acc(i)-star_acc));
%         jump_size = min(jump_size,1);
%         disp(['- - jump_size updated to ', num2str(jump_size)])
% 
%         % Update counters
%         i = i+1;
%         idx = idx+Na;
%     end

    % Current likelihood
    LL(k) = ll;
    XR(k,:) = reshape(xr,NN,1);

    % Update proposal covariance
    n = N_burnin+k-1;
    m = n/(n+1)*m+1/(n+1)*xr;
    s = s + xr.^2;
    a = 1/(n+1)*s - m.^2;
    a = min(a,lambda);
    
    % Print progression every 1/1000'th time increment
    b = ceil(N_iter/1000);
    if ~mod(k,b)
        disp([' - - ', num2str(k), '/', num2str(N_iter)]);
        disp(['- - ', 'Rejections: ', num2str(N_reject)]);
        n = k/b;
        acc_since_last_time = sum(acc(((n-1)*b+1):(n*b)))/b*100;
        disp(['- - ', 'Accepted: ', sprintf('%.1f%%',acc_since_last_time)]);
    end  
    
end
