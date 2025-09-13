function [bc, w, f, bs, psd] = calc_bicoher_mtm(X, Fs, NFFT, frange, compact, average, overlap, NW, isseg, normalization)
    [N, Ntrial] = size(X); % ONLY FOR UNIVARIATE
    if nargin < 10 || isempty(normalization)
        normalization = 'haubrich';
    end

    if nargin < 9 || isempty(isseg)
        isseg = true;
    end
    if nargin < 8 || isempty(NW)
        NW = 1;
    end

    if nargin < 7 || isempty(overlap)
        overlap = 0;
    end

    if nargin < 6 || isempty(average)
        average = true;
    end

    if nargin < 5 || isempty(compact)
        compact = false;
    end

    if nargin < 4 || isempty(frange)
        frange = [0, Fs / 2];
    end

    if nargin < 3 || isempty(NFFT)
        NFFT = N;
    end
    if nargin < 2 || isempty(Fs)
        Fs = 1;
    end

    K = round(2 * NW); % PMTM is K<(2NWdt-1) here seems different
    if isseg
        Ns = ceil(N / NFFT); % number of segments in which the EEG signal is wrapped
        Ns = max(1, Ns);
        X = [X; zeros(NFFT * Ns - N, Ntrial)]; % zero padding
        X = reshape(X, NFFT, Ns * Ntrial); % 'resized' EEG data
        X = permute(X, [1, 3, 2]);
        [vk, lambda] = dpss(NFFT, NW); % discrete prolate spheroidal (Slepian) sequences
        X = reshape(X .* vk, [NFFT, K * Ns * Ntrial]);
        M = NFFT;
    else
        Ns = 1;
        X = permute(X, [1, 3, 2]);
        [vk, lambda] = dpss(N, NW);
        X = reshape(X .* vk, [N, K * Ns * Ntrial]);
        M = N;
    end

    %% fft
    [sraw, wraw, fraw] = calc_fft(permute(X, [3, 1, 2]), false, false, Fs, NFFT);
    sraw = permute(sraw, [2, 1, 3]) ./ sqrt(2); %./M;

    %% normalize psd

    if strcmpi(normalization, 'shahbazi') %Three
        U2 = (M * K) ^ 2;
        psdraw = abs(sraw).^3;
        psdraw = (sum(psdraw, 3) ./ U2 ./ (Ns * Ntrial)).^(1 / 3);
    else
        U2 = M * K;
        psdraw = sraw .* conj(sraw);
        psdraw = sum(psdraw, 3) ./ U2 ./ (Ns * Ntrial);
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
        NFFTraw = NFFT;
        NFFT = length(f);
    end
    Ki = get_ndgrid(repmat((1:K).', 1, 3), 'list');
    K3 = size(Ki, 1);

    %% sx and sy
    s = reshape(s, [NFFT, 1, K, Ns * Ntrial]);
    sx = s(:, :, Ki(:, 1), :); %N,1,Ntrial
    sy = permute(s(:, :, Ki(:, 2), :), [2, 1, 3, 4]);
    bs = sx .* sy;
    if compact
        [~, idx] = mat2tril(bs, 0, true);
        bss = NaN(size(bs));
        bss(idx) = bs(idx);
        bs = bss;
    end

    %% s(x+y)
    sraw = reshape(sraw, NFFTraw, 1, K, Ns * Ntrial);
    sraw = sraw(:, :, Ki(:, 3), :);

    sf = griddedInterpolant(fraw, conj(sraw), 'nearest', 'none');
    [fx, fy] = ndgrid(f);
    bsxy = reshape(sf(reshape(fx + fy, NFFT * NFFT, 1)), NFFT, NFFT, K3, Ns * Ntrial); %flip(fx,1)+fy %reshape(fx+fy,Nf*Nf,1))

    bs = bs .* bsxy;
    %% P
    P = sum(vk(:, Ki(:, 1)) .* vk(:, Ki(:, 2)) .* vk(:, Ki(:, 3)), 1);
    P = permute(P, [1, 3, 2]);
    bs = bs .* P;
    %% gamma (Thomson, 1989) or U3 in (Birkelund et al., 2001)

    bs = sum(bs, [3, 4]);
    gamma = M ^ 2 * sum(P(:).^2, 1);
    bs = bs ./ gamma ./ (Ns * Ntrial);

    %% bicoherence
    psdf = griddedInterpolant(fraw, conj(psdraw), 'linear', 'none');
    psdxy = reshape(psdf(reshape(fx + fy, NFFT * NFFT, 1)), NFFT, NFFT, 1);

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