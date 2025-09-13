function [x] = pearson3_rnd(alpha, beta, xi, Gamma, sz)

    %% Generate a random array with size sz with values between 0 and 1
    p = rand(sz);

    %% If Gamma = 0, the Normal distribution is used
    if Gamma == 0
        norm_inv = @(P, mu, sigma) sigma .* (-sqrt(2) .* erfcinv(2 .* P)) + mu;
        x = norm_inv(p, alpha, beta);
        return
    end

    %% Compute the inverse cdf of the matrix p to generate the random sample
    if Gamma > 0 %Case of positive skewness
        X = gammaincinv(p, alpha .* ones(size(p)));
        x = (X .* beta) + xi;
    else %Case of negative skewness
        X = gammaincinv(1 - p, alpha .* ones(size(p)));
        x = (-X .* beta) + xi;
    end
end