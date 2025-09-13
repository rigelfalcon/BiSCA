function [yq, L, dbg] = NwSmooth(x, y, h, xq)

    if nargin < 4 || isempty(xq)
        xq = x;
    end

    [n, dx] = size(x);
    [nq, ~] = size(xq);

    if size(x, 2) > 1
        x = reshape(x, [n, 1, dx]);
    end

    if size(xq, 2) > 1
        xq = reshape(xq, [nq, 1, dx]);
    end

    h = permute(h, [dx + 2:-1:1]);

    D = x - permute(xq, [2, 1, 3]);

    Dkn = squeeze(gaussian_kernel(D, h));

    yq = squeeze(sum(Dkn .* y, 1) ./ sum(Dkn, 1));
    yq = yq';
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