function A = wedge2full(x, keeptip)
    % mat2wedge()
    % need get the N from N_elemement
    % wedge2full(mat2wedge(randn(6),true))
    if nargin < 2 || isempty(keeptip)
        keeptip = true;
    end

    assert(isvector(x), 'now only support x to be a vector')

    Ne = length(x);
    if keeptip
        Ne = Ne - 1;
    end

    k1 = -1 + 2 * sqrt(1 + Ne);
    k2 = -1 + sqrt(1 + 4 * Ne);

    if isintegernum(k1)
        k = k1;
    else
        k = k2;
    end

    A = zeros(k);
    idx = idx_wedge(k);
    A(idx) = x;
end