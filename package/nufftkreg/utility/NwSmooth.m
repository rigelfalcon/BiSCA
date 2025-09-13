function [yq, L, dbg] = NwSmooth(x, y, h, xq)
    [x, x_mean, x_std] = zscore(x);
    if nargin < 4 || isempty(xq)
        xq = x;
    end
    xq = (xq - x_mean) ./ x_std;

    [n, dx] = size(x);
    [nq, ~] = size(xq);
    if ~isvector(xq)
        xq = reshape(xq, [nq, 1, dx]);
        x = reshape(x, [n, 1, dx]);
    end
    h = permute(h, [dx + 2:-1:1]);

    D = x - permute(xq, [2, 1, 3]);
    Dkn = gaussian_kernel(D, h);
    yq = permute(sum(Dkn .* y, 1) ./ sum(Dkn, 1), [2, 1, 3]);
    dbg.s = sum(Dkn, 1);
    if nargout > 1
        L = Dkn ./ sum(Dkn, 1).'; %
    end
end

function u = epan_kernel(u, b)
    u = u ./ b;
    u = max(0, 3/4 * (1 - sum(u .^ 2, 3)));
end

function u = gaussian_kernel(u, b)
    u = u ./ b;
    u = (1 / sqrt(2 * pi)) * exp(-0.5 * sum((u .^ 2), 3));
end