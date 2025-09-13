function [index, res] = find_closest(A, a, n)
    if nargin < 3
        n = 1;
    end
    if n == 1
        index = arrayfun(@(x) minidx(abs(A - x)), a);
        res = A(index) - a;
    else
        index = arrayfun(@(x) sortidx(abs(A - x), n), a, 'UniformOutput', false);
        index = cell2mat(index);
        res = A(index) - a;
    end

    function idx = minidx(x)
        [~, idx] = min(x);
    end

    function idx = sortidx(x, n)
        [~, idx] = sort(x);
        idx = idx(1:n);
    end

end