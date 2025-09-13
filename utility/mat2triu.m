function [mattriu, idx] = mat2triu(data, k, isvec, onlyidx, isfirstslice)
    % examp: 10*10*100, get upper triangle
    % Yingwang 5/10/2020
    %  ind= 0 is the main diagonal
    %  ind > 0 is above the main diagonal
    %  ind < 0 is below the main diagonal.
    if nargin < 5 || isempty(isfirstslice)
        isfirstslice = false;
    end
    if nargin < 4 || isempty(onlyidx)
        onlyidx = false;
    end
    if nargin < 3 || isempty(isvec)
        isvec = true;
    end
    if nargin < 2 || isempty(k)
        k = 0;
    end
    [r, c, p] = size(data);

    if isvec
        if ~isfirstslice
            idx = triu(true(r, c), k);
            idx = find(idx);
            idx = idx(:) + r * c * (0:p - 1);
            if ~onlyidx
                mattriu = reshape(data(idx), [], p);
            end
        else
            idx = triu(true(r, c), k);
            idx = find(idx);
            if ~onlyidx
                data = reshape(data, [r * c, p]);
                mattriu = reshape(data(idx, :), [], p);
            end
        end
    else
        if ~onlyidx
            mattriu = zeros(size(data));
            for i = 1:p
                mattriu(:, :, i) = triu(data(:, :, i), k);
            end
        end
        if nargout > 1
            idx = triu(true(r, c), k);
            idx = repmat(idx, [1, 1, p]);
        end
    end
    if onlyidx
        mattriu = [];
    end

end