
function postm = get_tc(xx, ffmat, xgrid, rhoff, lenff)
% Compute tuning curve F at xgrid given xx and ffmat using robust GP regression

ffTYPE = 2;
covfun = covariance_fun(rhoff, lenff, ffTYPE);

% 1. Mean-center ffmat (Explicit GP Prior Mean)
% Prevents weight explosion when regressing uncentered log-rates
mu_ff = mean(ffmat, 1);
ffmat_centered = ffmat - mu_ff;

% 2. Compute covariance matrix with adaptive jitter
[K0, ~] = covfun(xx, xx);
n_samples = size(xx, 1);
jitter = max(1e-4, 1e-3 * mean(diag(K0))); % Adaptive diagonal stabilization
K = K0 + jitter * eye(n_samples);

% 3. Numerically stable solve using Cholesky decomposition
try
    L = chol(K, 'lower');
    alpha = L' \ (L \ ffmat_centered);
catch
    % Fallback to pseudo-inverse if matrix is near-singular
    alpha = pinv(K) * ffmat_centered;
end

% 4. Compute test-train covariance
K1 = covfun(xgrid, xx);

% 5. Predict log-rates and re-add baseline mean
postm = K1 * alpha + mu_ff;

% 6. Physiological upper/lower bound safety clip on log-rates
% Prevents exp(ff) from blowing up during downstream evaluation
max_log_rate = max(ffmat(:)) + 2; % ~7x max observed rate limit
min_log_rate = min(ffmat(:)) - 2;
postm = min(max(postm, min_log_rate), max_log_rate);
end

% function postm = get_tc(xx,ffmat,xgrid,rhoff,lenff)
% % Compute the F of xgrid based on xx and ffmat
% 
% ffTYPE = 2;
% covfun = covariance_fun(rhoff,lenff,ffTYPE); % get the covariance function
% 
% [K0,dK] = covfun(xx,xx);
% K = K0+1e-4*eye(size(xx,1));
% K1 = covfun(xgrid,xx);
% 
% postm = K1*(K\ffmat);
