function ar_noise = ar_noise(Nt, Nvar, type_inno)
    % ADD_AR_NOISE Add AR(1) noise to a signal
    %   YT_NOISY = ADD_AR_NOISE(YT, FS, P, NOISE_SCALE) adds AR(1) noise to
    %   the input signal YT with sampling frequency FS, using an AR model of
    %   order P and a noise scaling factor NOISE_SCALE.

    % Generate AR(1) noise with the same length as yt
    if nargin < 3 || isempty(type_inno)
        type_inno = 'gaussian';
    end

    ar_coeff = 0.85; % 0.99 AR coefficient %
    ar_noise = sqrt(1 - ar_coeff ^ 2) .* randn(Nt, Nvar); % white noise

    for i = 2:Nt
        if isCharString(type_inno)
            switch type_inno
                case 'gaussian'
                    ar_noise(i, :) = ar_coeff .* ar_noise(i - 1, :) + sqrt(1 - ar_coeff ^ 2) .* randn(1, Nvar);
                case 'pearson3'
                    ar_noise(i, :) = ar_coeff .* ar_noise(i - 1, :) + sqrt(1 - ar_coeff ^ 2) .* gen_pearson3_rnd(0, 1, 6, [1, Nvar]);
            end
        elseif isnumeric(type_inno)
            ar_noise(i, :) = ar_coeff .* ar_noise(i - 1, :) + sqrt(1 - ar_coeff ^ 2) .* type_inno(i);
        elseif isa(type_inno, 'prob.KernelDistribution')
            ar_noise(i, :) = ar_coeff .* ar_noise(i - 1, :) + sqrt(1 - ar_coeff ^ 2) .* random(pd, 1, 1);
        end
    end

    % Normalize the AR noise
    ar_noise = (ar_noise - mean(ar_noise, 1)) ./ std(ar_noise, 1);

end