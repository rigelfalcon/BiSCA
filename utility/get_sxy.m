function sxy = get_sxy(s, f, opt, fxy)
    if nargin < 4 || isempty(fxy)
        fxy = [];
    end
    if nargin < 3 || isempty(opt)
        opt = struct();
    end
    if nargin < 2 || isempty(f)
        f = 0:length(s) - 1;
    end
    if ~isfield(opt, 'InterpolationMethod') || isempty(opt.ExtrapolationMethod)
        opt.InterpolationMethod = 'nearest'; %linear
    end
    if ~isfield(opt, 'ExtrapolationMethod') || isempty(opt.ExtrapolationMethod)
        opt.ExtrapolationMethod = 'none';
    end

    f = f(:);
    N = length(f);
    Ncomp = size(s, 2);
    sf = griddedInterpolant(f, conj(s), opt.InterpolationMethod, opt.ExtrapolationMethod); % for bispectrum,use conj(s) as default
    if isempty(fxy)
        [fx, fy] = ndgrid(f);
        sxy = reshape(sf(reshape(fx + fy, N * N, 1)), N, N, Ncomp); %flip(fx,1)+fy %reshape(fx+fy,Nf*Nf,1))
    else
        sxy = reshape(sf(fxy), N, N, Ncomp);
    end

end