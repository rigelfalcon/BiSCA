function [alpha, beta, xi] = get_pearson3_para(mu, Sigma, Gamma)
sigma = sqrt(Sigma);
if Gamma~=0
    % Compute the parameters alpha, beta, and xi of the Pearson Type III distribution
    alpha = 4 / Gamma^2;
    beta = 0.5 * sigma * abs(Gamma);
    xi = mu - 2 * sigma / Gamma;
else
    alpha=mu;
    beta=sigma;
    xi=[];
end