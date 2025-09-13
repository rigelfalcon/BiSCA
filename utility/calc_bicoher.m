function [bc, w, f, bs, psd] = calc_bicoher(X, Fs, normalization, nfft, frange, compact, window, overlap)
    [N, Ntrial] = size(X);

    if nargin < 8 || isempty(window)
        window = N;
    end

    if ischar(window) || isstring(window)
        switch lower(window)
            case {'hanning'}
                window = hanning(N);
            case {'none'}
                window = [];
            otherwise
                error('not finish')
        end
    elseif isscalar(window)
        window = hanning(window);
    end
    Nlw = length(window); %number of length of window
    idx = 1:(Nlw * (1 - overlap)):(N - Nlw + 1);
    Ns = length(idx);
    idx = (idx(:) + (0:(Nlw - 1))).';
    idx = idx(:) + N * (0:(Ntrial - 1));
    idx = reshape(idx, [Nlw, Ns * Ntrial]);
    X = X(idx);
    [N, Ntrial] = size(X);

    if nargin < 9 || isempty(overlap)
        overlap = 0;
    end

    if nargin < 7 || isempty(compact)
        compact = false;
    end

    if nargin < 6 || isempty(frange)
        frange = [0, Fs / 2];
    end

    if nargin < 5 || isempty(nfft)
        nfft = N;
    end
    if nargin < 4 || isempty(normalization)
        normalization = 'haubrich';
    end

    if nargin < 2 || isempty(Fs)
        Fs = 1;
    end

    if ~isempty(window)
        X = X .* window;
    end

    [sraw, wraw, fraw] = calc_fft(permute(X, [3, 1, 2]), false, false, Fs, nfft);

    sraw = permute(sraw, [2, 1, 3]) ./ N;
    %% normalize psd
    if strcmpi(normalization, 'shahbazi') %Three
        psdraw = abs(sraw).^3;
        psdraw = (mean(psdraw, 3)).^(1 / 3);
    else
        psdraw = mean(sraw .* conj(sraw), 3);
    end
    %% take spectrum in the range
    if ~isempty(frange)
        if ~compact && frange(1) >= 0
            frange = [-frange(2), frange(2)];
        end
        fidx = find(fraw >= frange(1) & fraw <= frange(2));
        f = fraw(fidx);
        w = wraw(fidx);
        s = sraw(fidx, :, :);
        psd = psdraw(fidx);
        nfft = length(f);
    end
    bs = s .* permute(s, [2, 1, 3]);
    if compact
        [~, idx] = mat2tril(bs, 0, false);
        bs(setdiff(1:nfft * nfft * Ntrial, find(idx))) = NaN;
    end

    sf = griddedInterpolant(fraw, conj(sraw), 'nearest', 'none');
    [fx, fy] = ndgrid(f);
    bsxy = reshape(sf(reshape(fx + fy, nfft * nfft, 1)), nfft, nfft, Ntrial); %flip(fx,1)+fy %reshape(fx+fy,Nf*Nf,1))

    bs = bs .* bsxy;
    bs = mean(bs, 3);

    psdf = griddedInterpolant(fraw, psdraw, 'linear', 'none');
    psdxy = reshape(psdf(reshape(fx + fy, nfft * nfft, 1)), nfft, nfft, 1);

    switch lower(normalization)
        case {'haubrich'}
            % Haubrich normalization which has no upper bound
            % [see Birkelund et al. (2001) and Haubrich (1965)]
            bc = bs ./ (sqrt(psd .* psd.' .* psdxy));
        case {'skewness'}
            % skewness function [Collis et al. (1998)]
            bc = abs(bs).^2 ./ (psd .* psd.' .* psdxy);
        case {'shahbazi'}
            % Shahbazi et al. (2014)
            bc = bs ./ (psd .* psd.' .* psdxy);
        otherwise
            error('no such method')
    end
end