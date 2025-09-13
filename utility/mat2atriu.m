function [matatriu, idx] = mat2atriu(data, k, isvec)
    % examp: 10*10*100, get upper anti-triangle
    % Yingwang 5/10/2020
    %  k= 0 is the main diagonal
    %  k > 0 is above the main diagonal
    %  k < 0 is below the main diagonal.

    if nargin < 3 || isempty(isvec)
        isvec = true;
    end
    if nargin < 2 || isempty(k)
        k = 0;
    end
    [r, c, p] = size(data);

    if isvec
        idx = atriu(true(r, c), k);
        idx = find(idx);
        idx = idx(:) + r * c * (0:p - 1);
        matatriu = reshape(data(idx), [], p);
    else
        matatriu = zeros(size(data));
        for i = 1:p
            matatriu(:, :, i) = atriu(data(:, :, i), k);
        end
        if nargout > 1
            idx = atriu(true(r, c), k);
            idx = repmat(idx, [1, 1, p]);
        end
    end

end