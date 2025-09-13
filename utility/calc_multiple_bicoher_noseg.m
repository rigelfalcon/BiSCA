function [bc, w, f, M3, M2, denomM3, X, Xraw, varM3, varM2, resM2, resM3] = calc_multiple_bicoher_noseg(x, Fs, NFFT, frange, compact, normalization, returnseg, intrinsic_mean_s)

    [N, Ns, Nx] = size(x);

    if nargin < 6 || isempty(compact)
        compact = false;
    end

    if nargin < 5 || isempty(frange)
        frange = [0, Fs / 2];
    end

    if nargin < 4 || isempty(NFFT)
        NFFT = N;
    end

    if nargin < 3 || isempty(normalization)
        normalization = 'haubrich';
    end

    if nargin < 2 || isempty(Fs)
        Fs = 1;
    end

    %%
    [X, w, f] = calc_fft(permute(x, [3, 1, 2]), false, false, Fs, NFFT);

    X = permute(X, [2, 3, 4, 1]) ./ N;
    %% normalize psd
    if strcmpi(normalization, 'shahbazi') %Triple
        M2 = abs(X) .^ 3;
        M2 = (mean(M2, 2)) .^ (1/3);
    else
        M2raw = X .* conj(X); %M2

        if intrinsic_mean_s %Geometric mean fits the distribution of spectrum
            M2raw = log(M2raw + eps);
            meanM2raw = mean(M2raw, 2, 'omitnan');

            if returnseg
                M2 = exp(M2raw);
            else
                M2 = exp(meanM2raw);
            end

        else
            meanM2raw = mean(M2raw, 2, 'omitnan');

            if returnseg
                M2 = M2raw; %M2
            else
                M2 = meanM2raw; %M2
            end

        end

        resM2 = M2raw - meanM2raw;

        if returnseg
            varM2 = var(resM2, [], 2, 'omitnan'); % for Ns average not for single segment
        else
            varM2 = var(resM2, [], 2, 'omitnan') ./ Ns;
        end

    end

    %% take spectrum in the range
    if ~isempty(frange)

        if ~compact && frange(1) >= 0
            frange = [-frange(2), frange(2)];
        end

        fidx = find(f >= frange(1) & f <= frange(2));
        f = f(fidx);
        w = w(fidx);
        X = X(fidx, :, :, :);
        M2 = M2(fidx, :, :, :);
        varM2 = varM2(fidx, :, :, :);
        NFFT = length(f);
    end

    Xraw = [];
    M3 = permute(X, [1, 3, 2, 4]) .* permute(X, [3, 1, 2, 4]);

    if compact
        [~, idxM3] = mat2tril(M3, 1, true, true, true);
        M3 = reshape(M3, NFFT * NFFT, Ns, Nx);
        M3(~idxM3, :, :) = NaN;
        M3 = reshape(M3, NFFT, NFFT, Ns, Nx);
    end

    Xf = griddedInterpolant(f, conj(X), 'nearest', 'none');
    [fx, fy] = ndgrid(f);
    Xxy = reshape(Xf(reshape(fx + fy, NFFT * NFFT, 1)), NFFT, NFFT, Ns, Nx);

    M3temp = M3 .* Xxy;

    if returnseg
        M3 = M3temp;
    else
        M3 = mean(M3temp, 3);
    end

    denomM3 = mean(abs(M3temp), [3, 4]);
    resM3 = M3temp - mean(M3temp, 3);
    varM3 = mean(abs(resM3) .^ 2, 3);
    clear M3temp;
    varM3 = permute(varM3, [1, 2, 4, 3]);

    M2f = griddedInterpolant(f, conj(M2), 'nearest', 'none');
    M2xy = reshape(M2f(reshape(fx + fy, NFFT * NFFT, 1)), NFFT, NFFT, size(M2, 2), Nx);

    switch lower(normalization)
        case {'haubrich'}
            % Haubrich normalization which has no upper bound
            % [see Birkelund et al. (2001) and Haubrich (1965)]
            denomM3 = (sqrt(permute(M2, [1, 3, 2, 4]) .* permute(M2, [3, 1, 2, 4]) .* M2xy));
            bc = M3 ./ denomM3;
        case {'skewness'}
            % skewness function [Collis et al. (1998)]
            denomM3 = ((permute(M2, [1, 3, 2]) .* permute(M2, [3, 1, 2]) .* M2xy));
            bc = abs(M3) .^ 2 ./ denomM3;
        case {'shahbazi'}
            % Shahbazi et al. (2014)
            bc = M3 ./ (M2 .* permute(M2, [2, 1, 3]) .* M2xy);
            denomM3 = [];
        case {'hagihira', 'hag'}
            % hagihira 2001
            bc = M3 ./ denomM3;
        otherwise
            error('no such method')
    end

    M2 = squeeze(M2);

    compare_hosa_with_theoretical(M2, M3, varM2, varM3, bc, f, x, Nx, Ns, N);
end

function compare_hosa_with_theoretical(M2hat, M3hat, varM2emp, varM3emp, bc, f, x, Nx, Ns, N)
    meanM2hat = M2hat;
    meanM2hat(1) = NaN;
    meanM2hat = mean(abs(meanM2hat), [1], 'omitnan');

    meanM3hat = M3hat;
    meanM3hat(:, 1) = NaN;
    meanM3hat(1, :) = NaN;
    meanM3hat = abs(mean((meanM3hat), [1, 2], 'omitnan'));

    meanvarM2emp = mean(varM2emp(:), 'omitnan');
    meanvarM3emp = mean(varM3emp(:), 'omitnan');

    M2hat(1) = 0;
    varM2hat = mean(M2hat, 2) .^ 2;

    %

    [M2xM2yM2xy] = get_sxsysxy(f, M2hat);
    varM3hat = N * mean(M2xM2yM2xy, 3);

    bchat = bc;
    bchat(:, 1, :, :) = NaN;
    bchat(1, :, :, :) = NaN;
    bcreal = abs(mean(bchat, [1, 2], 'omitnan'));

    x = reshape(x, [], Nx);
    gamma = cumulants(x, 3);
    gamma = abs(gamma);
    gamma_bc = abs(gamma(3) ./ (gamma(2) .^ (3/2)));

    meanM2hat = mean(meanM2hat) * (N); %; %/sqrt(pi)
    meanM3hat = mean(meanM3hat) * (N ^ 2); %; %/pi
    bcreal = mean(bcreal) * ((N ^ 2) / (N ^ (3/2))); %; %((N^(3/2))/(N^2)); %((N^2)/(N^(3/2))); %/pi

    %
    disp(['M2/gamma(2): ', num2str(meanM2hat ./ gamma(2))])
    disp(['M3/gamma(3): ', num2str(meanM3hat ./ (gamma(3)))])
    disp(['bc/gamma_bc: ', num2str(bcreal ./ (gamma_bc))])
    disp(['bcreal: ', num2str(bcreal)])
    disp(['gamma_bc: ', num2str(gamma_bc)])
    disp(newline);

    idx_use = varM2hat ~= 0;
    disp(['varM2emp: ', num2str(meanvarM2emp)])
    disp(['varM2hat: ', num2str(mean(varM2hat(:), 'omitnan'))])
    disp(['varM2emp/varM2hat: ', num2str(mean((varM2emp(idx_use)) ./ (varM2hat(idx_use)), 'omitnan'))])
    disp(newline);

    idx_use = varM3hat ~= 0;
    disp(['varM3emp: ', num2str((meanvarM3emp))])
    disp(['varM3hat: ', num2str((mean(varM3hat(:), 'omitnan')))])
    disp(['M3MTvaremp/M3MTvarhat: ', num2str(10 * log10(mean(varM3emp(idx_use) ./ (varM3hat(idx_use)), 'omitnan')))])

end