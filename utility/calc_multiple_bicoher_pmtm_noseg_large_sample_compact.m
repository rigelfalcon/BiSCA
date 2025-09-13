function [bc, w, f, M3, M2, denomM3, X, Xraw, varM3, varM2, resM2, resM3] = calc_multiple_bicoher_pmtm_noseg_large_sample_compact(x, Fs, NFFT, frange, compact, NW, normalization, taper, droplast, returnseg, log_var_s)
    x = single(x);

    if true %test_gpu(4e8, [], true)
        x = convert_to_gpuarray(x);
    end

    [N, Ns, Nx] = size(x);
    isgpu = isgpuarray(x);

    if nargin < 11 || isempty(log_var_s)
        log_var_s = false;
    end

    if nargin < 10 || isempty(returnseg)
        returnseg = false;
    end

    if nargin < 9 || isempty(droplast)
        droplast = false;
    end

    if nargin < 8 || isempty(taper)
        taper = 'Slepian'; %Sine
    end

    if nargin < 7 || isempty(normalization)
        normalization = 'haubrich';
    end

    if nargin < 6 || isempty(NW)
        NW = 1; %2*NW = BW/f0 = BW*M*dt
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

    reg = 1e-6;

    K = round(2 * NW);
    x = permute(x, [1, 4, 2, 3]);

    switch taper
        case 'Slepian'
            vk = dpss(N, NW, K);
        case 'Sine'
            vk = sine_taper(N, NW);
        case 'SineScaled'
            vk = sine_taper_scaled(N, NW);
        otherwise
            vk = dpss(N, NW, K);
    end

    if droplast
        vk = vk(:, 1:end - 1);
        K = K - 1;
    end

    if Nx == 1
        xorig = x;
    end

    x = permute(x .* vk, [1, 3, 2, 4]);
    x = reshape(x, [N, Ns * K, Nx]);
    M = N;

    % mixing
    Ki = get_ndgrid(repmat((1:K).', 1, 3), 'list');

    %% fft
    if isgpu
        nxloop = 5000;
    else
        nxloop = 50000;
    end

    idx = gen_loopgroup(Nx, nxloop);
    num_loop = numel(idx);

    if num_loop > 1
        x = gather(x);
    end

    for iloop = num_loop:-1:1
        ix = idx{iloop};
        xloc = x(:, :, ix);
        [Xloc, w, f] = calc_fft(permute(xloc, [3, 1, 2]), false, false, Fs, NFFT);
        X(ix, :, :) = gather(Xloc);
    end

    clear x xloc Xloc; %
    X = permute(X, [2, 3, 1]) ./ N; %./ Vk_scale

    %% M2
    if strcmpi(normalization, 'shahbazi') %Triple
        U2 = M * K;
        M2 = (sum(abs(X), 2) ./ U2 * Ns) .^ (1/3);
    else
        U2 = K / (M);
        M2raw = X .* conj(X) ./ U2; %M2
        M2raw = reshape(M2raw, [], Ns, K, Nx);
        M2raw = sum(M2raw, 3, 'omitnan'); %mean over K taper

        if returnseg
            M2 = M2raw; %M2
        else
            M2 = mean(M2raw, 2); %M2
        end

        meanM2raw = mean(M2raw, 2, 'omitnan');

        % residual spectrum
        if log_var_s
            resM2 = log(M2raw) - log(meanM2raw);
            resM2(isinf(resM2)) = NaN;
        else
            resM2 = M2raw - meanM2raw;
        end

        % variance spectrum
        if returnseg
            varM2 = var(resM2, [], 2, 'omitnan'); % for Ns average not for single segment
        else
            varM2 = var(resM2, [], 2, 'omitnan') ./ Ns;
        end

    end

    % take quantity in the range
    if ~isempty(frange)

        if ~compact && frange(1) >= 0
            frange = [-frange(2), frange(2)];
        end

        fidx = (f >= frange(1) & f <= frange(2));
        f = f(fidx);
        w = w(fidx);
        X = X(fidx, :, :, :);
        M2 = M2(fidx, :, :, :);
        meanM2 = meanM2raw(fidx, :, :, :);
        M2seg = M2raw(fidx, :, :, :);
        varM2 = varM2(fidx, :, :, :);
        resM2 = resM2(fidx, :);
        NFFT = length(f);
    end

    clear M2raw
    Xraw = [];
    [fx, fy] = ndgrid(f); % % fxy=reshape(fx+fy,NFFT*NFFT,1);

    %% P
    P = sum(vk(:, Ki(:, 1)) .* vk(:, Ki(:, 2)) .* vk(:, Ki(:, 3)), 1)';

    P = permute(P, [2, 3, 4, 1]);

    %% gamma (Thomson, 1989) or U3 in (Birkelund et al., 2001)
    U3 = (sum(P(:) .^ 2, 1) ./ M);

    %% M3
    % loop for multiple signal
    if isgpu

        if isunix
            nxloop = 30;
        elseif ispc
            nxloop = 1;
        else
            nxloop = 1;
        end

    else

        if isunix
            nxloop = 2;
        elseif ispc
            nxloop = 1;
        else
            nxloop = 1;
        end

    end

    idx = gen_loopgroup(Nx, nxloop);
    num_loop = numel(idx);

    if isgpu
        [P, Ki, U3] = convert_to_gpuarray(P, Ki, U3);
    end

    X = reshape(X, [NFFT, 1, Ns, K, Nx]);

    if returnseg
        M3 = NaN(NFFT, NFFT, Ns, Nx); %,'single'

        switch lower(normalization)
            case {'hagihira', 'hag'}
                denomM3 = NaN(NFFT, NFFT, Ns, Nx); %,'single'
            otherwise
        end

    else
        M3 = NaN(NFFT, NFFT, Nx); %,'single'

        switch lower(normalization)
            case {'hagihira', 'hag'}
                denomM3 = NaN(NFFT, NFFT, Nx); %,'single'

            otherwise
        end

    end

    if nargout > 9
        resM3 = NaN(NFFT, NFFT, Ns, Nx);
    end

    for iloop = num_loop:-1:1
        ix = idx{iloop};
        nxix = numel(ix);
        Xloc = X(:, :, :, :, ix);

        if isgpu

            if check_matlab_version('23.2')
                Xloc = convert_to_gpuarray(Xloc);
            end

        end

        if iloop == num_loop
            Xf = griddedInterpolant(f, conj(Xloc), 'nearest', 'none'); %nearest 'linear'
        else
            Xf.Values = conj(Xloc);
        end

        Xxy = reshape(Xf(reshape(fx + fy, NFFT * NFFT, 1)), NFFT, NFFT, Ns, K, nxix); %flip(fx,1)+fy

        if isgpu
            Xxy = convert_to_gpuarray(Xxy);
        end

        %% take the expectation and apply scale
        switch lower(normalization)
            case {'hagihira', 'hag'}
                M3loc = P ./ U3 .* Xloc(:, :, :, Ki(:, 1), :) ...
                    .* permute(Xloc(:, :, :, Ki(:, 2), :), [2, 1, 3, 4, 5]) ...
                    .* Xxy(:, :, :, Ki(:, 3), :); %clear Xxy
                %%{
                if returnseg
                    denomM3(:, :, :, ix) = sum(abs(M3loc), [4]);
                else
                    denomM3(:, :, ix) = sum(abs(M3loc), [3, 4]) ./ Ns;
                end

                %%}
                M3loc = sum(M3loc, 4);

            otherwise
                M3loc = sum(P ./ U3 .* Xloc(:, :, :, Ki(:, 1), :) ...
                    .* permute(Xloc(:, :, :, Ki(:, 2), :), [2, 1, 3, 4, 5]) ...
                    .* Xxy(:, :, :, Ki(:, 3), :), 4); %clear Xxy
        end

        M3locmean = mean(M3loc, 3);

        if returnseg
            M3(:, :, :, ix) = gather(squeeze(M3loc));
        else
            M3(:, :, ix) = gather(squeeze(M3locmean));
        end

        if nargout > 9
            M3locres = M3loc - M3locmean; %clear M3loc
            if returnseg
                M3varemp = gather(squeeze(real(sum(M3locres .* conj(M3locres), 3) / Ns))); % for Ns average not for single segment
            else
                M3varemp = gather(squeeze(real(sum(M3locres .* conj(M3locres), 3) / Ns)) ./ Ns); % for Ns average not for single segment
            end
            varM3(:, :, ix) = M3varemp;
            resM3(:, :, :, ix) = M3locres;
            clear M3locres
        else
            resM3 = [];
            varM3 = [];
        end

    end

    clear M3loc Xxy
    %% bicoherence
    if returnseg
        M2d = M2;
    else

        M2d = mean(M2seg, 2);
    end

    %% M2xy

    switch lower(normalization)
        case {'hagihira', 'hag'}

        otherwise
            M2f = griddedInterpolant(f, conj(M2d), 'nearest', 'none');

            if returnseg
                M2xy = reshape(M2f(reshape(fx + fy, NFFT * NFFT, 1)), NFFT, NFFT, size(M2d, 2), Nx); %Ns->size(M2d,2)
            else
                M2xy = reshape(M2f(reshape(fx + fy, NFFT * NFFT, 1)), NFFT, NFFT, Nx);
            end

    end

    switch lower(normalization)
        case {'haubrich'}
            % Haubrich normalization which has no upper bound
            % [see Birkelund et al. (2001) and Haubrich (1965)]
            if returnseg
                denomM3 = sqrt(permute(M2d, [1, 3, 2, 4]) .* permute(M2d, [3, 1, 2, 4]) .* M2xy);
            else
                denomM3 = sqrt(permute(M2d, [1, 2, 4, 3]) .* permute(M2d, [2, 1, 4, 3]) .* M2xy);
            end

            bc = M3 ./ (denomM3 + reg);
        case {'skewness'}
            % skewness function [Collis et al. (1998)] square root of haubrich
            denomM3 = (sqrt(permute(M2d, [1, 3, 2, 4]) .* permute(M2d, [3, 1, 2, 4]) .* M2xy));
            bc = sqrt(abs(M3) .^ 2 ./ (denomM3 + reg));
        case {'shahbazi'}
            % Shahbazi et al. (2014)
            bc = M3 ./ (M2d .* permute(M2d, [2, 1, 3, 4]) .* M2xy + reg);
            denomM3 = [];
        case {'hagihira', 'hag'}
            % hagihira 2001
            bc = M3 ./ (denomM3 + reg);
        otherwise
            error('no such method')
    end

    if isgpu
        [bc, w, f, M3, M2, denomM3, X, Xraw, varM3, varM2, resM2, resM3] = convert_to_double(bc, w, f, M3, M2, denomM3, X, Xraw, varM3, varM2, resM2, resM3);
    end

end

function compare_hosa_with_theoretical(M2MThat, M3MThat, M2MTvaremp, M3MTvaremp, W2, W3, bchat, f, x, Nx, Ns, N, NW, K)
    meanM2MThat = M2MThat;
    meanM2MThat(1) = NaN;
    meanM2MThat = mean(abs(meanM2MThat), [1], 'omitnan');

    meanM3MThat = M3MThat;
    meanM3MThat(:, 1) = NaN;
    meanM3MThat(1, :) = NaN;
    meanM3MThat = abs(mean((meanM3MThat), [1, 2], 'omitnan'));

    meanM2MTvaremp = mean(M2MTvaremp(:), 'omitnan');
    meanM3MTvaremp = mean(M3MTvaremp(:), 'omitnan');

    M2MThat(1) = 0;
    M2xM2y = (M2MThat) .^ 2;
    M2MTvarhat = (1 / (Ns)) * M2xM2y;

    [M2xM2yM2xy] = get_sxsysxy(f, M2MThat);
    M3MTvarhat = (1 / Ns) * (4 / (3 * (2 * NW) ^ 2)) * M2xM2yM2xy;

    meanbcmthat = bchat;
    meanbcmthat(:, 1, :, :) = NaN;
    meanbcmthat(1, :, :, :) = NaN;
    meanbcmthat = abs(mean(meanbcmthat, [1, 2], 'omitnan'));

    x = reshape(x, [], Nx);
    gamma = cumulants(x, 3);
    gamma = abs(gamma);
    gamma_bc = abs(gamma(3) ./ (gamma(2) .^ (3/2)));

    meanM2MThat = mean(meanM2MThat) * N; %P2; %(N); %/sqrt(pi)
    meanM3MThat = mean(meanM3MThat) * (N ^ 2); %P3; %(N ^ 2); %/pi
    meanbcmthat = mean(meanbcmthat) * ((N ^ 2) / (N ^ (3/2))); %((N^(3/2))/(N^2)); %((N^2)/(N^(3/2))); %/pi
    disp(newline);

    %
    disp(['M2MThat/gamma(2): ', num2str(meanM2MThat ./ gamma(2))])
    disp(['M3MThat/gamma(3): ', num2str(meanM3MThat ./ (gamma(3)))])
    disp(['bcmthat/gamma_bc: ', num2str(meanbcmthat ./ (gamma_bc))])
    disp(['bcmthat: ', num2str(meanbcmthat)])
    disp(['gamma_bc: ', num2str(gamma_bc)])

    idx_use = M2MTvarhat ~= 0;
    disp(['M2MTvaremp: ', num2str(meanM2MTvaremp)])
    disp(['M2MTvarhat: ', num2str(mean(M2MTvarhat(:), 'omitnan'))])
    disp(['M2MTvaremp/M2MTvarhat: ', num2str(mean((M2MTvaremp(idx_use)) ./ (M2MTvarhat(idx_use)), 'omitnan'))])
    disp(newline);

    idx_use = M3MTvarhat ~= 0;
    disp(['M3MTvaremp: ', num2str((meanM3MTvaremp))])
    disp(['M3MTvarhat: ', num2str((mean(M3MTvarhat(:), 'omitnan')))])
    disp(['M3MTvaremp/M3MTvarhat: ', num2str((mean(M3MTvaremp(idx_use) ./ (M3MTvarhat(idx_use)), 'omitnan')))])

    disp(newline);

end