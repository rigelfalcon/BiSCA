function grid = get_ndgrid(x, type)
    if nargin < 2 || isempty(type)
        type = 'cell';
    end
    if ismatrix(x) && ~iscell(x)
        N = size(x, 1);
        dx = size(x, 2);
        x = mat2cell(x, N(1), ones(1, dx));
    else
        dx = numel(x);
    end
    grid = cell(dx, 1);
    [grid{1:dx}] = ndgrid(x{:});

    if strcmp(type, 'array')
        grid = permute(grid, dx + 1:-1:1);
        grid = cell2mat2(grid);
    end
    if strcmp(type, 'list') && ~iscell(x)
        grid = permute(grid, dx + 1:-1:1);
        grid = cell2mat(grid);
        grid = reshape(grid, [], dx);
    elseif strcmp(type, 'list') && iscell(x)
        grid = cellfun(@(x) reshape(x, [], 1), grid, 'UniformOutput', false);
        grid = cat(2, grid{:});
    end
    if strcmp(type, 'cellList')
        grid = cellfun(@(x) reshape(x, [], 1), grid, 'UniformOutput', false);
    end
end