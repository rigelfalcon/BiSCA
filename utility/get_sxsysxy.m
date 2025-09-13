function [sxsysxy, f] = get_sxsysxy(f, s, opt, fxy)

    if nargin < 3 || isempty(opt)
        opt = struct();
    end

    if nargin < 4 || isempty(fxy)
        fxy = [];
    end

    opt = set_defaults(opt, 'ExtrapolationMethod', 'none');
    opt = set_defaults(opt, 'doubleside', false);

    if opt.doubleside
        s = [flip(s(f > 0)); s];
        f = [flip(-f(f > 0)); f];
    end

    sxy = get_sxy(s, f, opt, fxy);
    if size(s, 1) == size(f, 1) && size(s, 2) > 1
        s = permute(s, [1, 3, 2]);
    end

    sxsysxy = s .* permute(s, [2, 1, 3]) .* sxy;
end