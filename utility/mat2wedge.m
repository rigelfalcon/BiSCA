function [matwedge, idx] = mat2wedge(data, isvec, keeptip)
    if nargin < 3 || (keeptip)
        keeptip = true;
    end
    if nargin < 2 || isempty(isvec)
        isvec = true;
    end

    [r, c, p] = size(data);

    idx = idx_wedge(r, keeptip);

    if isvec
        matwedge = data(idx);
    else
        matwedge = zeros(r);
        matwedge(idx) = data(idx);
    end

end