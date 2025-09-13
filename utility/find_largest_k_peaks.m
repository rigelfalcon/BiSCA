function [mu, h, baseline, spectrum_detrended] = find_largest_k_peaks(frequency, spectrum, k, f_min, f_max, min_distance, h_smooth, opt)
    frequency = frequency(:);
    spectrum = spectrum(:);

    if length(frequency) ~= length(spectrum)
        error('Frequency and spectrum vectors must be the same length.');
    end
    if k <= 0 || floor(k) ~= k
        error('k must be a positive integer.');
    end
    if nargin < 4 || isempty(f_min)
        f_min = min(frequency);
    end
    if nargin < 5 || isempty(f_max)
        f_max = max(frequency);
    end
    if nargin < 6 || isempty(min_distance)
        min_distance = (max(frequency) - min(frequency)) / length(frequency);
    end
    if nargin < 7 || isempty(h_smooth)
        h_smooth = 0.5;
    end
    if nargin < 8 || isempty(opt)
        opt = struct();
    end

    if isempty(opt) || ~isfield(opt, 'prominence_ratio') || isempty(opt.prominence_ratio)
        opt.isdefault_prominence_ratio = true;
    else
        opt.isdefault_prominence_ratio = false;
    end
    opt = set_defaults(opt, 'prominence_ratio', 0.05);
    try
        [baseline] = NwSmooth(frequency, spectrum, h_smooth, frequency);

    catch ME
        warning('Smoothing spline fitting failed: %s. Using zero baseline.', ME.message);
        baseline = zeros(size(spectrum));
    end

    spectrum_detrended = spectrum - baseline;

    detection_indices = (frequency >= f_min) & (frequency <= f_max);
    if ~any(detection_indices)
        error('No data points found within the specified frequency range [f_min, f_max].');
    end
    freq_detection = frequency(detection_indices);
    spectrum_detection = spectrum_detrended(detection_indices);

    min_prominence = max(spectrum_detection) * 0.05;

    [peak_heights, peak_locs] = findpeaks(spectrum_detection, freq_detection, ...
        'MinPeakProminence', min_prominence, ...
        'MinPeakDistance', min_distance, 'NPeaks', k);
    if length(peak_heights) < k && opt.isdefault_prominence_ratio
        [peak_heights, peak_locs] = findpeaks(spectrum_detection, freq_detection, ...
            'MinPeakDistance', min_distance, 'NPeaks', k);
    end
    if isempty(peak_heights)
        warning('No peaks found in the specified frequency range.');
        mu = [];
        h = [];
        return;
    end

    [sorted_heights, sorted_mu] = deal(peak_heights, peak_locs);

    num_peaks_found = length(sorted_heights);
    if k > num_peaks_found
        warning('Requested k=%d exceeds the number of available peaks (%d). Returning all peaks.', k, num_peaks_found);
        k = num_peaks_found;
    end

    mu = sorted_mu(1:k);
    h = sorted_heights(1:k);
end

function [mu, h] = findLargestKPeaksWeightedPolyfit(frequency, spectrum, k, min_distance, p_order)
    n = length(frequency);
    weights = ones(n, 1);
    end_fraction = 0.05; % 5% at each end
    num_end = max(1, ceil(n * end_fraction));
    weights(1:num_end) = 10;
    weights(end - num_end + 1:end) = 10;

    X = polyval_matrix(frequency, p_order);
    W = diag(weights);
    coeffs = (X' * W * X) \ (X' * W * spectrum);

    baseline = polyval(coeffs, frequency);
    spectrum_detrended = spectrum - baseline;

    min_prominence = max(spectrum_detrended) * 0.05;

    [peak_heights, peak_locs] = findpeaks(spectrum_detrended, frequency, ...
        'MinPeakProminence', min_prominence, ...
        'MinPeakDistance', min_distance);

    if isempty(peak_heights)
        warning('No peaks found after polynomial baseline correction.');
        mu = [];
        h = [];
        return;
    end

    [sorted_heights, sort_order] = sort(peak_heights, 'descend');
    sorted_mu = peak_locs(sort_order);

    num_peaks_found = length(sorted_heights);
    if k > num_peaks_found
        warning('Requested k=%d exceeds the number of available peaks (%d). Returning all peaks.', k, num_peaks_found);
        k = num_peaks_found;
    end

    mu = sorted_mu(1:k);
    h = sorted_heights(1:k);
end

function X = polyval_matrix(x, degree)
    X = ones(length(x), degree + 1);
    for d = 1:degree
        X(:, d + 1) = x.^d;
    end
end

function [yq, L, dbg] = NwSmooth(x, y, h, xq)
    if nargin < 4 || isempty(xq)
        xq = x;
    end

    [n, dx] = size(x);
    [nq, ~] = size(xq);

    if size(x, 2) > 1
        x = reshape(x, [n, 1, dx]);
    end

    if size(xq, 2) > 1
        xq = reshape(xq, [nq, 1, dx]);
    end

    h = permute(h, [dx + 2:-1:1]);

    D = x - permute(xq, [2, 1, 3]);

    Dkn = squeeze(gaussian_kernel(D, h));

    yq = squeeze(sum(Dkn .* y, 1) ./ sum(Dkn, 1));
    yq = yq';
    dbg.s = sum(Dkn, 1);

    if nargout > 1
        L = Dkn ./ sum(Dkn, 1).';
    end
end

function u = epan_kernel(u, b)
    u = u ./ b;
    u = max(0, 3/4 * (1 - sum(u .^ 2, 3)));
end

function u = gaussian_kernel(u, b)
    u = u ./ b;
    u = (1 / sqrt(2 * pi)) * exp(-0.5 * sum((u .^ 2), 3));
end