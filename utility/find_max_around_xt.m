function [maxVal, maxCoord] = find_max_around_xt(x, y, xt, radius, opt)
    % x: grids, must increase, don't have be equal dist
    % y: n-dimensional array
    % xt: array of coordinates, where each row represents a set of coordinates
    % radius: maximum distance from each coordinate to consider

    if nargin < 5 || isempty(opt)
        opt = [];
    end
    if ~isfield(opt, 'maxopt')
        opt.maxopt = 'abs';
    end
    if ~isfield(opt, 'norm')
        opt.norm = 'l2';
    end

    % Compute the number of dimensions of the array
    ndimsy = ndims(y);
    if ndimsy == 2 && size(y, 2) == 1
        ndimsy = 1;
    end
    xsz = size(x);
    ndc = get_ndcolon(ndimsy);

    % Compute the number of coordinates
    nxt = size(xt, 1);
    for iDim = ndimsy:-1:1
        xi{iDim} = x(ndc{:}, iDim);
        xi{iDim} = unique(xi{iDim}(:), 'rows', 'sorted');
    end

    % Compute the minimum and maximum indices for each dimension
    xt_min = xt - radius;
    xt_max = xt + radius;
    for i = nxt:-1:1
        % Create a cell array of index vectors for each dimension
        idx = cell(ndimsy, 1);
        for iDim = ndimsy:-1:1

            idx_min(:, iDim) = find_closest(xi{iDim}, xt_min(i, iDim));
            idx_max(:, iDim) = find_closest(xi{iDim}, xt_max(i, iDim));
            idx{iDim} = idx_min(:, iDim):idx_max(:, iDim);
        end
        idx = get_ndgrid(idx, 'cellList');
        ind = sub2ind(xsz, idx{:});

        % Use ndgrid to create a grid of indices for each dimension
        x_in_range = [];
        for iDim = ndimsy:-1:1
            xnd_in_range = x(ndc{:}, iDim);
            x_in_range(:, iDim) = xnd_in_range(ind);
        end
        y_in_range = y(ind);

        if strcmpi(opt.norm, 'l1')
            % do nothing for l1
            in_radius_idx = true(length(y_in_range), 1);
        else % default l2
            % Compute the distances between the coordinates and the grid points
            distances = pdist2(x_in_range, xt(i, :));

            % Find the indices of the points within the radius for each coordinate
            in_radius_idx = distances <= radius;
        end

        % Compute the maximum value over all points within the radius for each coordinate

        if strcmpi(opt.maxopt, 'abs')

            y_in_range = y_in_range(in_radius_idx);
            [~, idx] = max(abs(y_in_range), [], 1);
            maxVal_itarget = y_in_range(idx);
            x_in_range = x_in_range(in_radius_idx, :);
            maxCoord_itarget = x_in_range(idx, :);

        else
            [maxVal_itarget, idx] = max(y_in_range(in_radius_idx), [], 1);
            x_in_range = x_in_range(in_radius_idx, :);
            maxCoord_itarget = x_in_range(idx, :);
        end

        if isempty(maxVal_itarget)
            maxVal(i, :) = NaN;
            maxCoord(i, :) = NaN(1, ndimsy);
        else
            maxVal(i, :) = maxVal_itarget;
            maxCoord(i, :) = maxCoord_itarget;
        end
    end

end