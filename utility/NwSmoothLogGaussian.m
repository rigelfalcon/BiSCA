
function [yq, L, dbg] = NwSmoothLogGaussian(x, y, h, xq)
% Log-Gaussian kernel smoother with numerical stability

if nargin < 4 || isempty(xq)
    xq = x;
end

[n, dx] = size(x);
[nq, ~] = size(xq);
y=permute(y,[1,3,2]);
% Dimension management
if size(x,2) > 1
    x = reshape(x, [n,1,dx]);
end
if size(xq,2) > 1
    xq = reshape(xq, [nq,1,dx]);
end

h = permute(h, [dx+2:-1:1]);
D = x - permute(xq, [2,1,3]);

% Log-Gaussian kernel computation
log_kn = -0.5 * sum((D./h).^2, 3);

% Numerical stability
max_log = max(log_kn, [], 1);
log_sum = log(sum(exp(log_kn - max_log), 1)) + max_log;
weights = exp(log_kn - log_sum);

% Prediction calculation
yq = permute((sum(weights .* y, 1) ./ sum(weights, 1)),[2,3,1]);

% Debug outputs

if nargout > 1
    dbg.s = sum(weights, 1);
    L = weights;
end
end