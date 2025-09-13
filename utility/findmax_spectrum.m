function [smax, fmax] = findmax_spectrum(s, f, range)
    idx = find(f >= range(1) & f <= range(2));
    f = f(idx, :);
    s = s(idx, :);
    [smax, idx_max] = max(s);
    fmax = f(idx_max);

end