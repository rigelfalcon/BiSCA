function X = gen_pearson3_rnd(mu, sigma, gamma, sz)

    % Set the desired mean, variance, and skewness

    % Calculate the Pearson parameters
    [alpha, beta, xi] = get_pearson3_para(mu, sigma, gamma);

    %Generate a random sample using the Pearson 3 distribution
    X = pearson3_rnd(alpha, beta, xi, gamma, sz);

end