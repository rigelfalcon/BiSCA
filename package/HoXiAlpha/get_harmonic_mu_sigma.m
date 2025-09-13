function [mu, sigma, order] = get_harmonic_mu_sigma(mu, sigma, ks, kh, no_zero_peak)
    order = [1 ./ ((ks:-1:1) + 1), (1:kh)];

    mu = mu(1) .* order;

    if ~no_zero_peak
        mu = [0, mu];
    end

    if no_zero_peak
        sigma = [sigma(1) * ones(1, ks + kh)];
    else
        sigma = [sigma(1) * ones(1, 1 + ks + kh)];
    end

end