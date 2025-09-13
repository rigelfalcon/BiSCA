function [Hlist] = get_hList(M, range, space)

    dh = length(M);
    if nargin < 2 || isempty(range)
        range = repmat([0.01, 1], dh);
    end
    if nargin < 3 || isempty(space)
        space = @linspace;
    end

    Hlist = multispace(range(:, 1), range(:, 2), M, space)';
    Hlist = get_ndgrid(Hlist, 'list');
end