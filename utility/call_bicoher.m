function [bc, f, psd, w] = call_bicoher(y, opt) %,f,bs,psd
    [N, Ntrial] = size(y);
    opt = set_defaults(opt, 'nfft', N / 2);
    opt = set_defaults(opt, 'wind', []);
    opt = set_defaults(opt, 'nsamp', []);
    opt = set_defaults(opt, 'overlap', []);
    opt = set_defaults(opt, 'frange', [0, 40]);
    opt = set_defaults(opt, 'normalization', 'shahbazi');

    [bc, fraw, psdraw] = bicoher(y, opt.nfft, opt.wind, opt.nsamp, opt.overlap);
    fraw = fraw .* opt.Fs;
    wraw = fraw * 2 * pi;
    frange = opt.frange;

    [fx, fy] = ndgrid(fraw);
    bcf = griddedInterpolant(fx, fy, bc, 'nearest', 'none');

    if ~isempty(frange)
        if ~opt.compact && frange(1) >= 0
            frange = [-frange(2), frange(2)];
        end
        fidx = find(fraw >= frange(1) & fraw <= frange(2));
        f = fraw(fidx);
        w = wraw(fidx);
        psd = psdraw(fidx);
        NFFTraw = opt.nfft;
        NFFT = length(f);
        [fx, fy] = ndgrid(f);
        bc = bcf(fx, fy);
    end
    if opt.compact
        [~, idx] = mat2tril(bc, 0, false);
        bc(setdiff(1:NFFT * NFFT * 1, find(idx))) = NaN;
    end

end