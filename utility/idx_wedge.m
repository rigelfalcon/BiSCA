function [idx, Ne] = idx_wedge(N, keeptip)
    % special for the bispectrum
    % the openset of the wedge
    if nargin < 2 || (keeptip)
        keeptip = true;
    end

    idx = triu(true(N), 0);
    idx = flip(triu(flip(idx, 2), 0), 2);

    if isodd(N) && ~keeptip
        idx(floor(N / 2) + 1, floor(N / 2) + 1) = false;
    end

    if nargout > 1
        Ne = sum(idx, 'all');
    else
        Ne = [];
    end

end