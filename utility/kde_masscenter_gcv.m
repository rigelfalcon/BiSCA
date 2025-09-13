function [fq, selected_k, h_i, cv_scores] = kde_masscenter_gcv(x, kList, xq)
    % KDE_MASSCENTER_CV - Adaptive bandwidth KDE with LOO cross-validation
    % Uses Likelihood Cross-Validation (LCV) for bandwidth selection

    % Handle input defaults
    if nargin < 3 || isempty(xq), xq = x; end

    [n, dx] = size(x);
    max_k = max(kList);

    % Precompute neighbor distances for all k
    [~, dist] = knnsearch(x, x, 'K', max_k + 1); % Include self
    h_i_all = sqrt(dist(:, 2:end).^2 ./ dx); % Exclude self-distance

    % Compute LOO likelihood scores
    cv_scores = zeros(length(kList), 1);

    for i = 1:length(kList)
        k = kList(i);
        h_i = h_i_all(:, k);

        % Compute LOO density estimates
        f_loo = zeros(n, 1);
        for j = 1:n
            % Compute density at x(j) using all points except j
            mask = true(n, 1);
            mask(j) = false;
            f_loo(j) = kde_masscenter(x(mask, :), x(j, :), h_i(mask));
        end

        % Compute cross-validation score (avoid log(0))
        epsilon = 1e-10;
        cv_scores(i) = -mean(log(max(f_loo, epsilon)));
    end

    % Select optimal k
    [~, idx] = min(cv_scores);
    selected_k = kList(idx);
    h_i = h_i_all(:, selected_k);

    % Final density estimation
    fq = kde_masscenter(x, xq, h_i);
end

function fq = kde_masscenter(x, xq, h_i)
    % Improved KDE implementation with proper normalization
    [n, dx] = size(x);
    nq = size(xq, 1);

    % Normalization constant (Gaussian kernel)
    kernel_const = (2 * pi)^(-dx / 2);

    % Compute density estimates
    fq = zeros(nq, 1);
    for i = 1:nq
        % Scaled distances
        D = (x - xq(i, :)) ./ h_i;
        % Kernel evaluations
        kvals = kernel_const * exp(-0.5 * sum(D.^2, 2)) ./ (h_i.^dx);
        % Density estimate
        fq(i) = mean(kvals);
    end
end