function stat = calc_bc_stat(bc, type_stat, Ns, p_threshold, top_level)
    % Validate inputs
    if nargin < 5
        top_level = [];
    end
    validateattributes(bc, {'numeric'}, {'2d'});
    validateattributes(p_threshold, {'numeric'}, {'scalar', '>', 0, '<', 1});
    validateattributes(Ns, {'numeric'}, {'scalar', 'positive'});
    if strcmp(type_stat, 'top') && (isempty(top_level) || top_level <= 0 || top_level >= 100)
        error('top_level must be between 0 and 100 for ''top'' statistic.');
    end

    % Reshape and square bicoherence values
    bcvec = abs(reshape(bc, [], size(bc, 3))) .^ 2;

    % Compute test statistic and degrees of freedom
    [T, df, P] = compute_statistic(bcvec, type_stat, top_level);

    % Determine null distribution, compute p-value and threshold
    stat.T = T;
    stat.df = df;
    stat.Ns = Ns;
    stat.P = P;
    stat.p_threshold = p_threshold;

    if strcmp(type_stat, 'iqr')
        % For the 'iqr' test, rescale the squared bicoherence
        scaling_factor = 2 * Ns;
        scaled_bcvec = bcvec * scaling_factor;
        valid = scaled_bcvec(~isnan(scaled_bcvec));
        T_scaled = iqr(valid);

        % Compute sample mean from scaled values and derive noncentrality
        sample_mean = median(valid);
        lambda0 = max(sample_mean - 2, 0);

        q = safe_ncx2inv([0.25, 0.75], 2, lambda0, 1e4);
        q1 = q(1);
        q3 = q(2);
        theoretical_iqr = q3 - q1;
        f_q1 = safe_ncx2pdf(q1, 2, lambda0, 1e4);
        f_q3 = safe_ncx2pdf(q3, 2, lambda0, 1e4);

        n = numel(valid);
        var_iqr = (3 / (16 * n * f_q1 ^ 2) + 3 / (16 * n * f_q3 ^ 2) - 1 / (8 * n * f_q1 * f_q3));
        se_iqr = sqrt(var_iqr);

        z_score = abs(T_scaled - theoretical_iqr) / se_iqr;
        p_value = 2 * (1 - normcdf(z_score));
        p_value = min(p_value, 1);

        stat.th = theoretical_iqr / scaling_factor;
        stat.bc_th = sqrt(stat.th / scaling_factor);
        stat.bc_T = sqrt(T_scaled / scaling_factor);
        stat.p = p_value;
        stat.zscore = z_score;
        pd = [];

    else
        % Scale bicoherence values to estimate non-centrality (for applicable cases)
        scaled_bcvec = 2 * Ns * bcvec;
        valid_scaled = scaled_bcvec(~isnan(scaled_bcvec));
        [~, idx_outlier] = rmoutliers(valid_scaled, "quartiles", "ThresholdFactor", 1);
        sample_mean = median(valid_scaled(~idx_outlier));
        lambda0 = max(sample_mean - 2, 0); % Estimate non-centrality

        switch type_stat
            case 'sum'
                scaled_T = T * 2 * Ns;
                scaling_factor = 2 * Ns;
                p = 1 - safe_ncx2cdf(scaled_T, 2 * P, 0, 1e4); %P*lambda0
                th = safe_ncx2inv(1 - p_threshold, 2 * P, 0, 1e4);
                pd = struct('type', 'noncentralchisquare', 'df', 2 * P, 'lambda', P * lambda0);

            case {'mean', 'robustmean', 'logmean', 'robustlogmean', 'mean_median'}
                scaled_T = T * 2 * Ns * P;
                scaling_factor = 2 * Ns * P;
                p = 1 - safe_ncx2cdf(scaled_T, 2 * P, 0, 1e4); %P*lambda0
                th = safe_ncx2inv(1 - p_threshold, 2 * P, 0, 1e4);
                pd = struct('type', 'noncentralchisquare', 'df', 2 * P, 'lambda', P * lambda0);

            case 'median'
                scaled_T = T * 2 * Ns;

                population_median = safe_ncx2inv(0.5, 2, 0, 1e4);
                f_median = safe_ncx2pdf(population_median, 2, 0, 1e4);
                asymptotic_variance_of_median = 1 / (4 * P * f_median ^ 2);
                pd = makedist('Normal', 'mu', population_median, 'sigma', sqrt(asymptotic_variance_of_median));
                scaling_factor = 2 * Ns;
                p = 1 - pd.cdf(scaled_T);
                th = pd.icdf(1 - p_threshold);
            case 'emedian'
                scaled_T = T * 2 * Ns; % å°†ä¸­ä½æ•°ç»Ÿè?¡é‡ç¼©æ”¾è‡³Ï‡Â?(2)å°ºåº¦
                pd = makeChiSqDist(df); % åˆ›å»ºè‡?ç”±åº¦ä¸?2çš„å¡æ–¹åˆ†å¸?
                scaling_factor = 2 * Ns;
                p = 1 - pd.cdf(scaled_T);
                th = pd.icdf(1 - p_threshold);
            case 'max'
                b_P = safe_ncx2inv(1 - 1 / P, 2, lambda0, 1e4);
                f_b_P = safe_ncx2pdf(b_P, 2, lambda0, 1e4);
                a_P = 1 / (P * f_b_P);
                mu = b_P;
                sigma = a_P;
                pd = makedist('GeneralizedExtremeValue', 'k', 0, 'sigma', sigma, 'mu', mu);
                scaling_factor = 2 * Ns;
                scaled_T = T * scaling_factor;
                p = 1 - pd.cdf(scaled_T);
                th = pd.icdf(1 - p_threshold);

            case 'cmax'
                [~, idx_outlier] = rmoutliers(valid_scaled, "quartiles", "ThresholdFactor", 1);
                sample_mean = median(valid_scaled(~idx_outlier));
                lambda0 = max(sample_mean - 2, 0);
                theoretical_var = 4 + 4 * lambda0;
                scale_factor = sqrt(theoretical_var / 4);
                mu = scale_factor * (2 * log(P) - (log(log(P)) + log(4 * pi)) / 2) + lambda0;
                sigma = scale_factor * 2;
                pd = makedist('GeneralizedExtremeValue', 'k', 0, 'sigma', sigma, 'mu', mu);
                scaling_factor = 2 * Ns;
                scaled_T = T * scaling_factor;
                p = 1 - pd.cdf(scaled_T);
                th = pd.icdf(1 - p_threshold);

            case 'tmax'
                [~, idx_outlier] = rmoutliers(valid_scaled, "quartiles", "ThresholdFactor", 1);
                s2 = var(valid_scaled(~idx_outlier));
                P = sum(~idx_outlier);
                scale_factor = sqrt(s2 / 4);
                mu = scale_factor * (2 * log(P) - (log(log(P)) + log(4 * pi)) / 2);
                sigma = scale_factor * 2;
                pd = makedist('GeneralizedExtremeValue', 'k', 0, 'sigma', sigma, 'mu', mu);
                scaling_factor = 2 * Ns;
                scaled_T = T * scaling_factor;
                p = 1 - pd.cdf(scaled_T);
                th = pd.icdf(1 - p_threshold);

            case 'emax'
                scaled_T = T * 2 * Ns;
                p = 1 - safe_ncx2cdf(scaled_T, 2, 0, 1e4);
                th = safe_ncx2inv(1 - p_threshold, 2, 0, 1e4);
                pd = struct('type', 'noncentralchisquare', 'df', 2, 'lambda', lambda0);
                scaling_factor = 2 * Ns;

            case 'top'
                scaled_T = T * 2 * Ns;
                p = 1 - safe_ncx2cdf(scaled_T, 2, lambda0, 1e4);
                th = safe_ncx2inv(1 - p_threshold, 2, lambda0, 1e4);
                pd = struct('type', 'noncentralchisquare', 'df', 2, 'lambda', lambda0);
                scaling_factor = 2 * Ns;

            otherwise
                error('Unsupported statistic for parametric test.');
        end
        stat.T = scaled_T;
        stat.th = th;
        stat.bc_th = sqrt(stat.th / scaling_factor);
        stat.bc_T = sqrt(scaled_T / scaling_factor);
        stat.p = p;
        if isfield(pd, 'type') && strcmp(pd.type, 'noncentralchisquare')
            mean_val = pd.df + pd.lambda;
            std_val = sqrt(2 * pd.df + 4 * pd.lambda);
            stat.zscore = (scaled_T - mean_val) / std_val;
        else
            stat.zscore = (scaled_T - pd.mean()) / pd.std();
        end
    end

    stat.nlogp = -log10(stat.p);
    stat.h = stat.p < p_threshold;
    if strcmp(type_stat, 'iqr')
        mask = (abs(bc) .^ 2) < stat.th;
    else
        mask = (abs(bc)) > stat.bc_th;
    end
    bc(isnan(bc)) = 0;
    stat.bcth = sparse(double(bc .* mask));

    if strcmp(type_stat, 'max')
        stat.distribution = 'gumbel';
    else
        stat.distribution = 'gamma';
    end
    stat.pd = pd;
end

function [T, df, P] = compute_statistic(bcvec, type_stat, top_level)
    P = sum(~isnan(bcvec(:, 1)));

    switch type_stat
        case 'emax'
            valid = bcvec(~isnan(bcvec));
            T = max(valid, [], 1, 'omitnan');
            df = 2;
        case {'max', 'tmax'}
            valid = bcvec(~isnan(bcvec));
            T = sort(valid, 'descend');
            T = mean(T(1:2));
            df = 2;
        case 'top'
            T = prctile(bcvec, top_level, 1);
            df = 2;
        case 'mean'
            T = mean(bcvec, 1, 'omitnan');
            df = 2 * P;
        case {'mean_median'}
            T = median(bcvec, 1, 'omitnan');
            df = 2 * P;
        case {'median'}
            [~, idx_outlier] = rmoutliers(bcvec, "quartiles", "ThresholdFactor", 1); %3
            T = median(bcvec(~idx_outlier), 1, 'omitnan');
            P = sum(~idx_outlier);
            df = 2 * P;
        case {'emedian'}
            T = median(bcvec, 1, 'omitnan');
            df = 2;
        case 'robustmean'
            valid = bcvec(~isnan(bcvec));
            T = trimmean(valid, 10);
            df = 2 * P;
        case 'logmean'
            valid = bcvec(~isnan(bcvec));
            T = 10.^mean(log10(valid), 'omitnan');
            df = 2 * P;
        case 'robustlogmean'
            valid = log10(bcvec(~isnan(bcvec)));
            valid = valid(~isinf(valid));
            T = 10.^trimmean(valid, 10);
            df = 2 * P;
        case 'sum'
            T = sum(bcvec, 1, 'omitnan');
            df = 2 * P;
        case 'iqr'
            valid = bcvec(~isnan(bcvec));
            T = iqr(valid);
            df = [];
        otherwise
            error('Unsupported statistic.');
    end
end

function x = safe_ncx2inv(p, df, lambda0, lambda0_threshold)
    % SAFE_NCX2INV Wrapper for ncx2inv, using normal approximation for large lambda0.
    %
    % Inputs:
    %   p                 - Probability for inverse CDF (0 <= p <= 1).
    %   df                - Degrees of freedom for non-central chi-square.
    %   lambda0           - Non-centrality parameter (nonnegative).
    %   lambda0_threshold - Threshold for normal approximation (default: 1e4).
    %
    % Output:
    %   x                 - Quantile (inverse CDF) for probability p.

    if nargin < 4
        lambda0_threshold = 1e4;
    end

    validateattributes(p, {'numeric'}, {'vector', '>=', 0, '<=', 1});
    validateattributes(df, {'numeric'}, {'scalar', 'positive'});
    validateattributes(lambda0, {'numeric'}, {'scalar', 'nonnegative'});
    validateattributes(lambda0_threshold, {'numeric'}, {'scalar', 'positive'});

    if lambda0 > lambda0_threshold
        mu = df + lambda0;
        sigma = sqrt(2 * (df + 2 * lambda0));
        x = norminv(p, mu, sigma);
    else
        x = ncx2inv(p, df, lambda0);
    end
end

function f = safe_ncx2pdf(x, df, lambda0, lambda0_threshold)
    % SAFE_NCX2PDF Wrapper for ncx2pdf, using normal approximation for large lambda0.
    %
    % Inputs:
    %   x                 - Value(s) at which to evaluate the PDF.
    %   df                - Degrees of freedom for non-central chi-square.
    %   lambda0           - Non-centrality parameter (nonnegative).
    %   lambda0_threshold - Threshold for normal approximation (default: 1e4).
    %
    % Output:
    %   f                 - PDF value at x.

    if nargin < 4
        lambda0_threshold = 1e4;
    end

    validateattributes(x, {'numeric'}, {'vector', 'nonnegative'});
    validateattributes(df, {'numeric'}, {'scalar', 'positive'});
    validateattributes(lambda0, {'numeric'}, {'scalar', 'nonnegative'});
    validateattributes(lambda0_threshold, {'numeric'}, {'scalar', 'positive'});

    if lambda0 > lambda0_threshold
        mu = df + lambda0;
        sigma = sqrt(2 * (df + 2 * lambda0));
        f = normpdf(x, mu, sigma);
    else
        f = ncx2pdf(x, df, lambda0);
    end
end

function p = safe_ncx2cdf(x, df, lambda0, lambda0_threshold)
    % SAFE_NCX2CDF Wrapper for ncx2cdf, using normal approximation for large lambda0.
    %
    % Inputs:
    %   x                 - Value(s) at which to evaluate the CDF.
    %   df                - Degrees of freedom for non-central chi-square.
    %   lambda0           - Non-centrality parameter (nonnegative).
    %   lambda0_threshold - Threshold for normal approximation (default: 1e4).
    %
    % Output:
    %   p                 - CDF value at x.

    if nargin < 4
        lambda0_threshold = 1e4;
    end

    validateattributes(x, {'numeric'}, {'vector', 'nonnegative'});
    validateattributes(df, {'numeric'}, {'scalar', 'positive'});
    validateattributes(lambda0, {'numeric'}, {'scalar', 'nonnegative'});
    validateattributes(lambda0_threshold, {'numeric'}, {'scalar', 'positive'});

    if lambda0 > lambda0_threshold
        mu = df + lambda0;
        sigma = sqrt(2 * (df + 2 * lambda0));
        p = normcdf(x, mu, sigma);
    else
        p = ncx2cdf(x, df, lambda0);
    end
end