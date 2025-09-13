function taperinfo = update_taperinfo(taperinfo)
    N = taperinfo.N;
    nfft = taperinfo.nfft;
    N = 2 * taperinfo.N;
    nfft = 2 * taperinfo.nfft;
    Fs = taperinfo.Fs;
    K = round(2 * taperinfo.NW);

    if ~isfield(taperinfo, 'vk') || N ~= size(taperinfo.vk, 1)

        switch taperinfo.taper
            case {'Slepian', 'DPSS'}
                taperinfo.vk = dpss(N, taperinfo.NW, K);
            case 'Sine'
                taperinfo.vk = sine_taper(N, taperinfo.NW);
            case 'SineScaled'
                taperinfo.vk = sine_taper_scaled(N, taperinfo.NW);
            otherwise
                taperinfo.vk = dpss(N, taperinfo.NW, K);
        end

        if taperinfo.droplast
            taperinfo.vk = taperinfo.vk(:, 1:end - 1);
        end

        taperinfo.Vk = fftshift(fft(ifftshift(taperinfo.vk))); %
        taperinfo.Vk_e = sum(taperinfo.Vk, [1, 2]) ./ size(taperinfo.vk, 2);
        taperinfo.vk = permute(taperinfo.vk, [1, 3, 4, 2]);
        taperinfo.vk = padarray(taperinfo.vk, [nfft - N, 0], 0, 'post');
    end

    if ~isfield(taperinfo, 'f') || nfft ~= size(taperinfo.f, 1)

        if rem(nfft, 2) == 0
            f = [-nfft / 2:(nfft / 2 - 1)]' / nfft * Fs;
        else
            f = [- (nfft - 1) / 2:(nfft - 1) / 2]' / (nfft - 1) * Fs;
        end

        taperinfo.f = f;
    end

end